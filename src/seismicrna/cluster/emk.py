from itertools import permutations
from typing import Callable, Iterable

import numpy as np

from .em import EMRun
from ..core.logs import logger
from ..core.mu import calc_rmsd, calc_nrmsd, calc_pearson

NOCONV = 0


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
                 min_nrmsd_run: float,
                 max_jackpot_quotient: float,
                 max_loglike_vs_best: float,
                 min_pearson_vs_best: float,
                 max_nrmsd_vs_best: float):
        if not runs:
            raise ValueError(f"{self} got no EM runs")
        # Sort the runs from largest to smallest likelihood.
        runs = sort_runs(runs)
        # Flag runs that fail to meet the filters.
        self.max_pearson_run = max_pearson_run
        self.min_nrmsd_run = min_nrmsd_run
        self.max_jackpot_quotient = max_jackpot_quotient
        # Set the criteria for whether this number of clusters passes.
        self.max_loglike_vs_best = max_loglike_vs_best
        self.min_pearson_vs_best = min_pearson_vs_best
        self.max_nrmsd_vs_best = max_nrmsd_vs_best
        # Check whether each run shows signs of being overclustered.
        # To select only the valid runs, use "not" with the opposite of
        # the desired inequality because runs with just one cluster will
        # produce NaN values, which should always compare as True here.
        self.run_not_overclustered = np.array(
            [not (run.max_pearson > max_pearson_run
                  or run.min_nrmsd < min_nrmsd_run)
             for run in runs]
        )
        # Check whether each run shows signs of being underclustered.
        self.run_not_underclustered = np.array(
            [not (run.jackpot_quotient > max_jackpot_quotient)
             for run in runs]
        )
        # Select the best run.
        self.best = runs[self.best_index()]
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
        # Jackpotting quotient of each run.
        self.jackpot_quotients = np.array([run.jackpot_quotient
                                           for run in runs])
        # Minimum NRMSD between any two clusters in each run.
        self.min_nrmsds = np.array([run.min_nrmsd for run in runs])
        # Maximum correlation between any two clusters in each run.
        self.max_pearsons = np.array([run.max_pearson for run in runs])
        # NRMSD between each run and the best run.
        self.nrmsds_vs_best = np.array([calc_rms_nrmsd(run, self.best)
                                        for run in runs])
        # Correlation between each run and the best run.
        self.pearsons_vs_best = np.array([calc_mean_pearson(run, self.best)
                                          for run in runs])

    def run_passing(self, allow_underclustered: bool = False):
        """ Whether each run passed the filters. """
        if allow_underclustered:
            return self.run_not_overclustered
        return self.run_not_overclustered & self.run_not_underclustered

    def n_runs_passing(self, **kwargs):
        """ Number of runs passing the filters. """
        return int(np.count_nonzero(self.run_passing(**kwargs)))

    def get_valid_index(self, i: int | list[int] | np.ndarray, **kwargs):
        """ Index(es) of valid run number(s) `i`. """
        return np.flatnonzero(self.run_passing(**kwargs))[i]

    def best_index(self, **kwargs) -> int:
        """ Index of the best valid run. """
        try:
            # The best run is the valid run with the largest likelihood.
            return self.get_valid_index(0, **kwargs)
        except IndexError:
            # If no runs are valid, then use the best invalid run.
            logger.warning(f"{self} got no EM runs that passed all filters")
            return 0

    def subopt_indexes(self, **kwargs):
        """ Indexes of the valid suboptimal runs. """
        return self.get_valid_index(np.arange(1, self.n_runs_passing(**kwargs)),
                                    **kwargs)

    def loglike_vs_best(self, **kwargs):
        """ Log likelihood difference between the best and second-best
        runs. """
        try:
            index1, index2 = self.get_valid_index([0, 1], **kwargs)
            return float(self.log_likes[index1] - self.log_likes[index2])
        except IndexError:
            return np.nan

    def pearson_vs_best(self, **kwargs):
        """ Maximum Pearson correlation between the best run and any
        other run. """
        try:
            return float(np.max(
                self.pearsons_vs_best[self.subopt_indexes(**kwargs)]
            ))
        except ValueError:
            return np.nan

    def nrmsd_vs_best(self, **kwargs):
        """ Minimum NRMSD between the best run and any other run. """
        try:
            return float(np.min(
                self.nrmsds_vs_best[self.subopt_indexes(**kwargs)]
            ))
        except ValueError:
            return np.nan

    def enough_runs_passing(self, **kwargs):
        """ Whether enough runs passed. """
        return self.n_runs_passing(**kwargs) >= min(self.n_runs_total, 2)

    def passing(self, **kwargs):
        """ Whether this number of clusters passes the filters. """
        # Use "not" followed by the opposite of the desired inequality
        # so that if any attribute is NaN, the inequality will evaluate
        # to False, and its "not" will be True, as desired.
        return self.enough_runs_passing(**kwargs) and not (
                self.loglike_vs_best(**kwargs) > self.max_loglike_vs_best
                or self.pearson_vs_best(**kwargs) < self.min_pearson_vs_best
                or self.nrmsd_vs_best(**kwargs) > self.max_nrmsd_vs_best
        )

    def summarize(self, **kwargs):
        """ Summarize the results of the runs. """
        lines = [f"EM runs for K={self.k}",
                 "\nPARAMETERS\n"]
        for attr in ["max_pearson_run",
                     "min_nrmsd_run",
                     "max_jackpot_quotient",
                     "max_loglike_vs_best",
                     "min_pearson_vs_best",
                     "max_nrmsd_vs_best"]:
            lines.append(f"{attr} = {getattr(self, attr)}")
        lines.append("\nRUNS\n")
        for attr in ["n_runs_total",
                     "converged",
                     "log_likes",
                     "bics",
                     "jackpot_quotients",
                     "min_nrmsds",
                     "max_pearsons",
                     "nrmsds_vs_best",
                     "pearsons_vs_best",
                     "run_not_overclustered",
                     "run_not_underclustered"]:
            lines.append(f"{attr} = {getattr(self, attr)}")
        lines.append("\nPASSING\n")
        for attr in ["run_passing",
                     "n_runs_passing",
                     "best_index",
                     "loglike_vs_best",
                     "pearson_vs_best",
                     "nrmsd_vs_best",
                     "enough_runs_passing",
                     "passing"]:
            func = getattr(self, attr)
            lines.append(f"{attr} = {func(**kwargs)}")
        return "\n".join(lines)


def find_best_k(ks: Iterable[EMRunsK], **kwargs):
    """ Find the best number of clusters. """
    # Sort the runs by increasing numbers of clusters.
    ks = sorted(ks, key=lambda runs: runs.k)
    if not ks:
        logger.warning("No numbers of clusters exist")
        return 0
    # Select only the numbers of clusters that pass the filters.
    # For the largest number of clusters, underclustering can be allowed
    # to permit this number to be identified as the best number so far;
    # for all other numbers, use all filters.
    ks = [runs for runs in ks[:-1] if runs.passing()] + (
        [ks[-1]] if ks[-1].passing(**kwargs) else []
    )
    if not ks:
        logger.warning("No numbers of clusters pass the filters")
        return 0
    # Of the remaining numbers of clusters, find the number that gives
    # the smallest BIC.
    ks = sorted(ks, key=lambda runs: runs.best.bic)
    # Return that number of clusters.
    return ks[0].k


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
    mus1 = run1.mus.values
    mus2 = run2.mus.values
    nrmsds = calc_nrmsd_groups(mus1, mus2)
    assignment = assign_clusterings(mus1, mus2)
    return float(np.sqrt(np.mean([np.square(nrmsds[row, col])
                                  for row, col in enumerate(assignment)])))


def calc_mean_pearson(run1: EMRun, run2: EMRun):
    """ Compute the mean Pearson correlation between the clusters. """
    mus1 = run1.mus.values
    mus2 = run2.mus.values
    correlations = calc_pearson_groups(mus1, mus2)
    assignment = assign_clusterings(mus1, mus2)
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
