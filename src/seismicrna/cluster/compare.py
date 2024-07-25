import re
from itertools import permutations
from logging import getLogger
from typing import Callable

import numpy as np
import pandas as pd

from .em import EMRun
from .names import LOG_EXP_NAME, LOG_OBS_NAME
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


def filter_runs(runs: list[EMRun],
                max_pearson_run: float,
                min_nrmsd_run: float):
    """ Filter EM runs. """
    if not runs:
        logger.warning(f"Got no runs to filter")
    # Remove each number of clusters where any pair of clusters has a
    # Pearson correlation greater than the limit.
    if runs and max_pearson_run < 1.:
        runs = [run for run in runs
                if not run.calc_max_pearson() > max_pearson_run]
        if not runs:
            logger.warning(f"For all runs, the maximum Pearson correlation "
                           f"is > {max_pearson_run}")
    # Remove each number of clusters where any pair of clusters has an
    # NRMSD less than the limit.
    if runs and min_nrmsd_run > 0.:
        runs = [run for run in runs
                if not run.calc_min_nrmsd() < min_nrmsd_run]
        if not runs:
            logger.warning(f"For all runs, the minimum NRMSD "
                           f"is < {min_nrmsd_run}")
    return runs


class EMRunsK(object):
    """ One or more EM runs with the same number of clusters. """

    def __init__(self,
                 runs: list[EMRun],
                 max_pearson_run: float,
                 min_nrmsd_run: float):
        if not runs:
            raise ValueError("Got no clustering runs")
        runs = sort_runs(runs)
        # Number of clusters (k).
        self.k = get_common_k(runs)
        # Number of runs.
        self.n_runs_total = len(runs)
        # Number of iterations until convergenge for each run.
        self.converged = [run.iter if run.converged else NOCONV
                          for run in runs]
        # List of log-likelihoods for each run.
        self.log_likes = [run.log_like for run in runs]
        # Root-mean-square deviations between each run and run 0.
        self.nrmsds_vs_best = [calc_rms_nrmsd(run, runs[0])
                               for run in runs]
        # Correlations between each run and run 0.
        self.pearsons_vs_best = [calc_mean_pearson(run, runs[0])
                                 for run in runs]
        # Minimum NRMSD between any two clusters
        self.min_nrmsds = [run.calc_min_nrmsd() for run in runs]
        # Maximum Pearson correlation between any two clusters
        self.max_pearsons = [run.calc_max_pearson() for run in runs]
        # Remove runs for which any pair of clusters has an invalid
        # Pearson correlation or NRMSD.
        runs = filter_runs(runs,
                           max_pearson_run=max_pearson_run,
                           min_nrmsd_run=min_nrmsd_run)
        self.n_runs_filtered = len(runs)
        # Keep the remaining run with the best (largest) likelihood.
        self.best = runs[0] if runs else None

    @property
    def loglike_vs_best(self):
        """ Log likelihood difference between the best and second-best
        runs. """
        if self.n_runs_total > 1:
            return self.log_likes[0] - self.log_likes[1]
        # No difference if only 1 run.
        return 0.

    @property
    def pearson_vs_best(self):
        """ Maximum Pearson correlation between the best run and any
        other run. """
        return max(self.pearsons_vs_best[1:], default=1.)

    @property
    def nrmsd_vs_best(self):
        """ Minimum NRMSD between the best run and any other run. """
        return min(self.nrmsds_vs_best[1:], default=0.)

    @property
    def best_bic(self):
        """ BIC of the best run. """
        if self.best is None:
            # Infinitely bad (large) if no best run.
            return np.inf
        return self.best.bic


def filter_ks(ks: list[EMRunsK],
              max_loglike_vs_best: float,
              min_pearson_vs_best: float,
              max_nrmsd_vs_best: float):
    """ Filter numbers of clusters. """
    if not ks:
        logger.warning(f"Got no numbers of clusters to filter")
    else:
        # Remove each number of clusters with no runs.
        ks = [k for k in ks if k.n_runs_filtered > 0]
        if not ks:
            logger.warning(f"No number of clusters had at least one run")
    if ks:
        # Remove each number of clusters where the log likelihood gap
        # is too large.
        ks = [k for k in ks if k.loglike_vs_best <= max_loglike_vs_best]
        if not ks:
            logger.warning(f"No number of clusters had a log likelihood "
                           f"vs. the best run that was ≤ {max_loglike_vs_best}")
    if ks:
        # Remove each number of clusters where the maximum Pearson
        # correlation with the best run is too small.
        ks = [k for k in ks if k.pearson_vs_best >= min_pearson_vs_best]
        if not ks:
            logger.warning(f"No number of clusters had a Pearson correlation "
                           f"vs. the best run that was ≥ {min_pearson_vs_best}")
    if ks:
        # Remove each number of clusters where the minimum NRMSD with
        # the best run is too large.
        ks = [k for k in ks if k.nrmsd_vs_best <= max_nrmsd_vs_best]
        if not ks:
            logger.warning(f"No number of clusters had an NRMSD "
                           f"vs. the best run that was ≤ {max_nrmsd_vs_best}")
    return ks


def get_common_best_run_attr(ks: list[EMRunsK], attr: str):
    """ Get an attribute of the best clustering run from every k, and
    confirm that `key(attribute)` is identical for every k. """
    # Start by getting the attribute from the first k.
    value = getattr(ks[0].best, attr)
    # Verify that the best run from every other k has an equal value
    # of that attribute.
    if any(getattr(k.best, attr) != value for k in ks[1:]):
        raise ValueError(f"Found more than 1 value for attribute {repr(attr)} "
                         f"among k values {ks}")
    return value


def find_best_k(ks: list[EMRunsK],
                max_loglike_vs_best: float,
                min_pearson_vs_best: float,
                max_nrmsd_vs_best: float):
    """ Find the best number of clusters. """
    if not ks:
        raise ValueError("Got no groups of EM runs with any number of clusters")
    # Remove numbers of clusters that do not pass the filters.
    ks = filter_ks(ks,
                   max_loglike_vs_best=max_loglike_vs_best,
                   min_pearson_vs_best=min_pearson_vs_best,
                   max_nrmsd_vs_best=max_nrmsd_vs_best)
    if not ks:
        raise ValueError("No groups of EM runs with any numbers of clusters "
                         "passed the filters")
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


def get_log_exp_obs_counts(ks: list[EMRunsK]):
    """ Get the expected and observed log counts of each bit vector. """
    # Retrieve the unique bit vectors from the clusters.
    uniq_reads = get_common_best_run_attr(ks, "uniq_reads")
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_reads.counts_per_uniq),
                        index=uniq_reads.get_uniq_names())
    # For each number of clusters, compute the expected log counts.
    log_exp = ((format_exp_count_col(runs.k),
                pd.Series(runs.best.logn_exp, index=log_obs.index))
               for runs in ks)
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
