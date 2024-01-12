"""
Cluster Comparison Module
========================================================================
Auth: Matty

Collect and compare the results from independent runs of EM clustering.
"""

from __future__ import annotations

from itertools import permutations
from typing import Callable

import numpy as np
import pandas as pd

from .em import EmClustering
from .names import EXP_NAME, OBS_NAME
from ..core.header import ORDER_NAME
from ..core.mu import calc_nrmsd, calc_pearson

EXP_COUNT_PRECISION = 3  # Number of digits to round expected log counts
NON_CONVERGED = -1  # Number indicating a run did not converge


def get_common_order(runs: list[EmClustering]):
    """ Find the order of the clustering (the number of clusters) from a
    list of EM clustering runs of the same order. If multiple orders are
    found, then raise a ValueError. """
    orders = sorted(set(run.order for run in runs))
    if len(orders) != 1:
        raise ValueError(f"Expected 1 unique order, but got {orders}")
    return orders[0]


def sort_replicate_runs(runs: list[EmClustering]):
    """ Sort the runs of EM clustering by decreasing likelihood so that
    the run with the best (largest) likelihood comes first. """
    # Verify that every run has the same order.
    get_common_order(runs)
    return sorted(runs, key=lambda run: run.log_like, reverse=True)


class RunOrderResults(object):
    """ Results of clustering runs of the same order. """

    def __init__(self, runs: list[EmClustering]):
        if not runs:
            raise ValueError("Got no clustering runs")
        runs = sort_replicate_runs(runs)
        # Order of the clustering (i.e. number of clusters).
        self.order = get_common_order(runs)
        # Number of runs.
        self.n_runs = len(runs)
        # Number of iterations until convergenge for each run.
        self.converged = [run.iter if run.converged else NON_CONVERGED
                          for run in runs]
        # List of log likelihoods for each run.
        self.log_likes = [run.log_like for run in runs]
        # Root-mean-square deviations between adjacent runs.
        self.rmsds = [calc_rms_nrmsd(runs[i], runs[i + 1])
                      for i in range(len(runs) - 1)]
        # Correlations between adjacent runs.
        self.meanr = [calc_mean_pearson(runs[i], runs[i + 1])
                      for i in range(len(runs) - 1)]
        # Run with the best (smallest) BIC.
        self.best = runs[0]


def orders_list_to_dict(orders: list[RunOrderResults]):
    """ Return a dict of the RunOrderResults keyed by their orders. """
    runs_per_order = dict()
    for runs in orders:
        if runs.order in runs_per_order:
            raise ValueError(f"Duplicate order: {runs.order}")
        runs_per_order[runs.order] = runs
    for order in range(1, max(runs_per_order) + 1):
        if order not in runs_per_order:
            raise ValueError(f"Missing order: {order}")
    return runs_per_order


def get_common_best_run_attr(orders: list[RunOrderResults], attr: str):
    """ Get an attribute of the best clustering run from every order,
    and confirm that `key(attribute)` is identical for all orders. """
    # Start by getting the attribute from order 1, which always exists.
    value = getattr(orders_list_to_dict(orders)[1].best, attr)
    # Verify that the best run from every other order has an equal value
    # of that attribute.
    if any(getattr(order.best, attr) != value
           for order in orders if order.order != 1):
        raise ValueError(f"Found more than 1 value for attribute {repr(attr)} "
                         f"among EM clustering runs {orders}")
    return value


def find_best_order(orders: list[RunOrderResults]) -> int:
    """ Find the number of clusters with the best (smallest) BIC. """
    return sorted(orders, key=lambda runs: runs.best.bic)[0].order


def format_exp_count_col(order: int):
    return f"{EXP_NAME}, {ORDER_NAME} {order}"


def get_log_exp_obs_counts(orders: list[RunOrderResults]):
    """ Get the expected and observed log counts of each bit vector. """
    # Retrieve the unique bit vectors from the clusters.
    uniq_reads = get_common_best_run_attr(orders, "uniq_reads")
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_reads.counts_per_uniq),
                        index=uniq_reads.get_uniq_names())
    # For each order of clustering, compute the expected log counts.
    log_exp = ((format_exp_count_col(runs.order),
                pd.Series(runs.best.logn_exp, index=log_obs.index))
               for runs in orders)
    # Assemble all log counts into one DataFrame.
    log_counts = pd.DataFrame.from_dict({OBS_NAME: log_obs, **dict(log_exp)})
    # Sort the data by expected count at order 1, then round it.
    return log_counts.sort_values(by=[format_exp_count_col(order=1)],
                                  ascending=False).round(EXP_COUNT_PRECISION)


def _compare_groups(func: Callable, mus1: np.ndarray, mus2: np.ndarray):
    """ Compare two groups of clusters using a comparison function and
    return a matrix of the results. """
    _, n1 = mus1.shape
    _, n2 = mus2.shape
    return np.array([[func(mus1[:, cluster1], mus2[:, cluster2])
                      for cluster2 in range(n2)]
                     for cluster1 in range(n1)])


def calc_nrmsd_groups(mus1: np.ndarray, mus2: np.ndarray):
    """ Calculate the NRMSD of each pair of clusters in two groups. """
    return _compare_groups(calc_nrmsd, mus1, mus2)


def calc_pearson_groups(mus1: np.ndarray, mus2: np.ndarray):
    """ Calculate the Pearson correlation of each pair of clusters in
    two groups. """
    return _compare_groups(calc_pearson, mus1, mus2)


def assign_clusterings(mus1: np.ndarray, mus2: np.ndarray):
    """ Optimally assign clusters from two groups to each other. """
    # Make a cost matrix for assigning clusters using the MSD.
    costs = np.square(calc_nrmsd_groups(mus1, mus2))
    n, m = costs.shape
    if n != m:
        raise ValueError(f"Got different numbers of clusters in groups 1 ({n}) "
                         f"and 2 ({m})")
    # Find the assignment of clusters that gives the minimum cost using
    # the naive approach of checking every possible pairwise assignment.
    # While other algorithms (e.g. the Jonker-Volgenant algorithm) solve
    # the assignment problem in O(n³) time, and this naive approach runs
    # in O(n!) time, the latter is simpler and still sufficiently fast
    # when n is no more than 5 or 6, which is almost always true.
    best_assignment = list()
    min_cost = None
    for columns in permutations(range(n)):
        assignment = [(i, columns[i]) for i in range(n)]
        cost = sum(costs[i, j] for i, j in assignment)
        if min_cost is None or cost < min_cost:
            min_cost = float(cost)
            best_assignment = assignment
    return best_assignment


def calc_rms_nrmsd(run1: EmClustering, run2: EmClustering):
    """ Compute the root-mean-square NRMSD between the clusters. """
    costs = np.square(calc_nrmsd_groups(run1.mus, run2.mus))
    assignment = assign_clusterings(run1.mus, run2.mus)
    return float(np.sqrt(np.mean([costs[i, j] for i, j in assignment])))


def calc_mean_pearson(run1: EmClustering, run2: EmClustering):
    """ Compute the mean Pearson correlation between the clusters. """
    correlations = calc_pearson_groups(run1.mus, run2.mus)
    assignment = assign_clusterings(run1.mus, run2.mus)
    return float(np.mean([correlations[i, j] for i, j in assignment]))

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
