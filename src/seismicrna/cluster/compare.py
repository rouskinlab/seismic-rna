import re
from functools import cached_property
from itertools import permutations
from logging import getLogger
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .em import EMRun
from .names import LOG_EXP_NAME, LOG_OBS_NAME
from .uniq import UniqReads
from ..core.header import NUM_CLUSTS_NAME
from ..core.mu import calc_rmsd, calc_nrmsd, calc_pearson
from ..core.report import NOCONV

logger = getLogger(__name__)

EXP_COUNT_PRECISION = 3  # Number of digits to round expected log counts


def get_common_k(runs: list[EMRun]):
    """ Find the number of clusters (k) from among EM clustering runs.
    If there are multiple ks, then raise a ValueError. """
    ks: list[int] = sorted({run.k for run in runs})
    if len(ks) != 1:
        raise ValueError(f"Expected 1 unique number of clusters, but got {ks}")
    return ks[0]


def sort_runs(runs: list[EMRun]):
    """ Sort the runs of EM clustering by decreasing likelihood so that
    the run with the best (largest) likelihood comes first. """
    # Verify that every run has the same k; otherwise, the likelihood is
    # not directly comparable between runs.
    if runs:
        get_common_k(runs)
    return sorted(runs, key=lambda run: run.log_like, reverse=True)


class EMRunsK(object):
    """ One or more EM runs with the same number of clusters. """

    def __init__(self,
                 runs: list[EMRun],
                 max_pearson_run: float,
                 min_nrmsd_run: float):
        if not runs:
            raise ValueError(f"{self} got no EM runs")
        # Sort the runs from largest to smallest likelihood.
        runs = sort_runs(runs)
        # Flag runs that fail to meet the filters.
        self.max_pearson_run = max_pearson_run
        self.min_nrmsd_run = min_nrmsd_run
        # To filter invalid runs, need to use "not" with the opposite of
        # the desired inequality because runs with just one cluster will
        # produce NaN values, which should always compare as True here.
        self.run_valid = np.array([
            not (run.calc_max_pearson() > max_pearson_run
                 or run.calc_min_nrmsd() < min_nrmsd_run)
            for run in runs
        ])
        if self.n_runs_valid == 0:
            logger.warning(f"{self} got no EM runs where all pairs of clusters "
                           f"had Pearson correlation ≤ {self.max_pearson_run} "
                           f"and NRMSD ≥ {self.min_nrmsd_run}")
        # Select the best run.
        self.best = runs[self.best_index]
        # Number of runs.
        self.n_runs_total = len(runs)
        # Number of clusters (K).
        self.k = get_common_k(runs)
        # Number of iterations until convergenge for each run.
        self.converged = np.array([run.iter if run.converged else NOCONV
                                   for run in runs])
        # Log-likelihood of each run.
        self.log_likes = np.array([run.log_like for run in runs])
        # BIC of each run.
        self.bics = np.array([run.bic for run in runs])
        # Minimum NRMSD between any two clusters in each run.
        self.min_nrmsds = np.array([run.calc_min_nrmsd() for run in runs])
        # Maximum correlation between any two clusters in each run.
        self.max_pearsons = np.array([run.calc_max_pearson() for run in runs])
        # NRMSD between each run and the best run.
        self.nrmsds_vs_best = np.array([calc_rms_nrmsd(run, self.best)
                                        for run in runs])
        # Correlation between each run and the best run.
        self.pearsons_vs_best = np.array([calc_mean_pearson(run, self.best)
                                          for run in runs])
        self._passing = None

    @cached_property
    def n_runs_valid(self):
        """ Number of valid runs. """
        return int(np.count_nonzero(self.run_valid))

    def get_valid_index(self, i: int | list[int] | np.ndarray):
        """ Index(es) of valid run number(s) `i`. """
        return np.flatnonzero(self.run_valid)[i]

    @cached_property
    def best_index(self) -> int:
        """ Index of the best valid run. """
        if self.n_runs_valid > 0:
            # The best run is the valid run with the largest likelihood.
            return self.get_valid_index(0)
        # If no runs are valid, then use the best invalid run.
        return 0

    @cached_property
    def subopt_indexes(self):
        """ Indexes of the valid suboptimal runs. """
        return self.get_valid_index(np.arange(1, self.n_runs_valid))

    @cached_property
    def loglike_vs_best(self):
        """ Log likelihood difference between the best and second-best
        runs. """
        if self.n_runs_valid < 2:
            return np.nan
        index1, index2 = self.get_valid_index([0, 1])
        return float(self.log_likes[index1] - self.log_likes[index2])

    @cached_property
    def pearson_vs_best(self):
        """ Maximum Pearson correlation between the best run and any
        other run. """
        if self.n_runs_valid < 2:
            return np.nan
        return float(np.max(self.pearsons_vs_best[self.subopt_indexes]))

    @cached_property
    def nrmsd_vs_best(self):
        """ Minimum NRMSD between the best run and any other run. """
        if self.n_runs_valid < 2:
            return np.nan
        return float(np.min(self.nrmsds_vs_best[self.subopt_indexes]))

    @property
    def best_bic(self):
        """ BIC of the best run. """
        return self.best.bic

    @property
    def passing(self):
        """ Whether this number of clusters passes the filters. """
        if not isinstance(self._passing, bool):
            raise ValueError(f"{self}.valid has not been set")
        return self._passing

    def set_passing(self,
                    max_loglike_vs_best: float,
                    min_pearson_vs_best: float,
                    max_nrmsd_vs_best: float):
        """ Set whether this number of clusters passes the filters. """
        # Use "not" followed by the opposite of the desired inequality
        # so that if any attribute is NaN, the inequality will evaluate
        # to False, and its "not" will be True, as desired.
        self._passing = not (self.n_runs_valid == 0
                             or self.loglike_vs_best > max_loglike_vs_best
                             or self.pearson_vs_best < min_pearson_vs_best
                             or self.nrmsd_vs_best > max_nrmsd_vs_best)
        return self.passing


def find_best_k(ks: Iterable[EMRunsK]):
    """ Find the best number of clusters. """
    # Select only the numbers of clusters that pass the filters.
    ks = [runs for runs in ks if runs.passing]
    if not ks:
        logger.warning("No numbers of clusters pass the filters")
        return 0
    # Of the remaining numbers of clusters, find the number that gives
    # the smallest BIC.
    ks = sorted(ks, key=lambda k: k.best_bic)
    # Return that number of clusters.
    return ks[0].k


def format_exp_count_col(k: int):
    return f"{LOG_EXP_NAME}, {NUM_CLUSTS_NAME} {k}"


def parse_exp_count_col(col: str):
    if not (match := re.match(f"^{LOG_EXP_NAME}, {NUM_CLUSTS_NAME} ([0-9]+)$", col)):
        raise ValueError(f"Invalid expected count column: {repr(col)}")
    k, = match.groups()
    return int(k)


def get_log_exp_obs_counts(uniq_reads: UniqReads, ks: list[EMRunsK]):
    """ Get the expected and observed log counts of each bit vector. """
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_reads.counts_per_uniq),
                        index=uniq_reads.get_uniq_names())
    # For each number of clusters, compute the expected log counts.
    log_exp = [(format_exp_count_col(runs.k),
                pd.Series(runs.best.logn_exp, index=log_obs.index))
               for runs in ks]
    # Assemble all log counts into one DataFrame.
    log_counts = pd.DataFrame.from_dict({LOG_OBS_NAME: log_obs} | dict(log_exp))
    # Round the log counts.
    return log_counts.round(EXP_COUNT_PRECISION)


def _compare_groups(func: Callable, mus1: np.ndarray, mus2: np.ndarray):
    """ Compare two groups of clusters using a comparison function and
    return a matrix of the results. """
    _, n1 = mus1.shape
    _, n2 = mus2.shape
    return np.array([[func(mus1[:, cluster1], mus2[:, cluster2])
                      for cluster2 in range(n2)]
                     for cluster1 in range(n1)]).reshape((n1, n2))


def calc_rmsd_groups(mus1: np.ndarray, mus2: np.ndarray):
    """ Calculate the RMSD of each pair of clusters in two groups. """
    return _compare_groups(calc_rmsd, mus1, mus2)


def calc_nrmsd_groups(mus1: np.ndarray, mus2: np.ndarray):
    """ Calculate the NRMSD of each pair of clusters in two groups. """
    return _compare_groups(calc_nrmsd, mus1, mus2)


def calc_pearson_groups(mus1: np.ndarray, mus2: np.ndarray):
    """ Calculate the Pearson correlation of each pair of clusters in
    two groups. """
    return _compare_groups(calc_pearson, mus1, mus2)


def assign_clusterings(mus1: np.ndarray, mus2: np.ndarray):
    """ Optimally assign clusters from two groups to each other. """
    # Make a cost matrix for assigning clusters using the RMSD.
    costs = np.square(calc_rmsd_groups(mus1, mus2))
    n, m = costs.shape
    if n != m:
        raise ValueError(f"Got different numbers of clusters in groups 1 ({n}) "
                         f"and 2 ({m})")
    # Find the assignment of clusters that gives the minimum cost using
    # the naive approach of checking every possible pairwise assignment.
    # While other algorithms (e.g. the Jonker-Volgenant algorithm) solve
    # the assignment problem in O(n³) time, and this naive approach runs
    # in O(n!) time, the latter is simpler and still sufficiently fast
    # when n is no more than about 6, which is almost always true.
    ns = np.arange(n)
    best_assignment = ns
    min_cost = None
    for cols in permutations(ns):
        assignment = np.array(cols, dtype=int)
        cost = np.sum(costs[ns, assignment])
        if min_cost is None or cost < min_cost:
            min_cost = cost
            best_assignment = assignment
    return best_assignment


def calc_rms_nrmsd(run1: EMRun, run2: EMRun):
    """ Compute the root-mean-square NRMSD between the clusters. """
    nrmsds = calc_nrmsd_groups(run1.p_mut, run2.p_mut)
    assignment = assign_clusterings(run1.p_mut, run2.p_mut)
    return float(np.sqrt(np.mean([np.square(nrmsds[row, col])
                                  for row, col in enumerate(assignment)])))


def calc_mean_pearson(run1: EMRun, run2: EMRun):
    """ Compute the mean Pearson correlation between the clusters. """
    correlations = calc_pearson_groups(run1.p_mut, run2.p_mut)
    assignment = assign_clusterings(run1.p_mut, run2.p_mut)
    return float(np.mean([correlations[row, col]
                          for row, col in enumerate(assignment)]))

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
