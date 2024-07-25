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
    get_common_k(runs)
    return sorted(runs, key=lambda run: run.log_like, reverse=True)


class EMRunsK(object):
    """ One or more EM runs with the same number of clusters. """

    def __init__(self, runs: list[EMRun], max_pearson: float, min_nrmsd: float):
        if not runs:
            raise ValueError("Got no clustering runs")
        runs = sort_runs(runs)
        # Number of clusters (k).
        self.k = get_common_k(runs)
        # Number of runs.
        self.n_runs = len(runs)
        # Number of iterations until convergenge for each run.
        self.converged = [run.iter if run.converged else NOCONV
                          for run in runs]
        # List of log-likelihoods for each run.
        self.log_likes = [run.log_like for run in runs]
        # Root-mean-square deviations between each run and run 0.
        self.nrmsd_vs_0 = [calc_rms_nrmsd(run, runs[0]) for run in runs]
        # Correlations between each run and run 0.
        self.pearson_vs_0 = [calc_mean_pearson(run, runs[0]) for run in runs]
        # Minimum NRMSD between any two clusters
        self.min_nrmsds = [run.calc_min_nrmsd() for run in runs]
        # Maximum Pearson correlation between any two clusters
        self.max_pearsons = [run.calc_max_pearson() for run in runs]
        # Remove runs for which any pair of clusters has an invalid
        # Pearson correlation or NRMSD.
        runs = filter_runs(runs, max_pearson=max_pearson, min_nrmsd=min_nrmsd)
        # Keep the remaining run with the best (largest) likelihood.
        self.best = runs[0]

    @property
    def bic(self):
        return self.best.bic

    def calc_max_pearson(self):
        return self.best.calc_max_pearson()

    def calc_min_nrmsd(self):
        return self.best.calc_min_nrmsd()


def filter_runs(runs: list[EMRun] | list[EMRunsK],
                max_pearson: float,
                min_nrmsd: float):
    # Remove each number of clusters where any pair of clusters has a
    # Pearson correlation greater than the limit.
    if max_pearson < 1.:
        runs_use = [run for run in runs
                    if not run.calc_max_pearson() > max_pearson]
        if runs_use:
            runs = runs_use
        else:
            logger.warning(f"For all runs {runs}, the maximum Pearson "
                           f"correlation is > {max_pearson}")
    # Remove each number of clusters where any pair of clusters has an
    # NRMSD less than the limit.
    if min_nrmsd > 0.:
        runs_use = [run for run in runs
                    if not run.calc_min_nrmsd() < min_nrmsd]
        if runs_use:
            runs = runs_use
        else:
            logger.warning(f"For all runs {runs}, the minimum NRMSD "
                           f"is < {min_nrmsd}")
    return runs


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


def find_best_k(ks: list[EMRunsK], max_pearson: float, min_nrmsd: float):
    """ Find the best number of clusters. """
    if not ks:
        raise ValueError("Got no groups of EM runs with any number of clusters")
    # Remove runs for which any pair of clusters has an invalid Pearson
    # correlation or NRMSD.
    ks = filter_runs(ks, max_pearson=max_pearson, min_nrmsd=min_nrmsd)
    # Of the remaining numbers of clusters, find the number that gives
    # the smallest BIC.
    ks = sorted(ks, key=lambda k: k.bic)
    best = ks[0]
    return best.k


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
