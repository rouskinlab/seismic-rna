"""
Cluster Comparison Module
========================================================================
Auth: Matty

Collect and compare the results from independent runs of EM clustering.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from .em import EmClustering
from .names import EXP_NAME, OBS_NAME
from ..core.header import ORDER_NAME

EXP_COUNT_PRECISION = 3  # Number of digits to round expected log counts


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
        self.converged = [run.iter if run.converged else 0 for run in runs]
        # List of log likelihoods for each run.
        self.log_likes = [run.log_like for run in runs]
        # Log likelihood mean and standard deviation.
        self.log_like_mean = np.mean(self.log_likes)
        self.log_like_std = np.std(self.log_likes)
        # Variation of information.
        self.var_info = calc_mean_var_info(runs)
        # Run with the best (smallest) BIC.
        self.best = runs[0]


def get_common_best_run_attr(ord_runs: dict[int, RunOrderResults], attr: str):
    """ Get an attribute of the best clustering run from every order,
    and confirm that `key(attribute)` is identical for all orders. """
    # Start by getting the attribute from order 1, which always exists.
    value = getattr(ord_runs[1].best, attr)
    # Verify that the best run from every other order has an equal value
    # of that attribute.
    if any(getattr(runs.best, attr) != value
           for order, runs in ord_runs.items() if order != 1):
        raise ValueError(f"Found more than 1 value for attribute '{attr}' "
                         f"among EM clustering runs {ord_runs}")
    return value


def find_best_order(ord_runs: dict[int, RunOrderResults]) -> int:
    """ Find the number of clusters with the best (smallest) BIC. """
    return sorted(ord_runs.items(), key=lambda runs: runs[1].best.bic)[0][0]


def format_exp_count_col(order: int):
    return f"{EXP_NAME}, {ORDER_NAME} {order}"


def get_log_exp_obs_counts(ord_runs: dict[int, RunOrderResults]):
    """ Get the expected and observed log counts of each bit vector. """
    # Retrieve the unique bit vectors from the clusters.
    uniq_reads = get_common_best_run_attr(ord_runs, "uniq_reads")
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_reads.counts_per_uniq),
                        index=uniq_reads.get_uniq_names())
    # For each order of clustering, compute the expected log counts.
    log_exp = ((format_exp_count_col(order),
                pd.Series(runs.best.logn_exp, index=log_obs.index))
               for order, runs in ord_runs.items())
    # Assemble all log counts into one DataFrame.
    log_counts = pd.DataFrame.from_dict({OBS_NAME: log_obs, **dict(log_exp)})
    # Sort the data by expected count at order 1, then round it.
    return log_counts.sort_values(by=[format_exp_count_col(order=1)],
                                  ascending=False).round(EXP_COUNT_PRECISION)


def calc_information_entropy(x: np.ndarray, validate: bool = True):
    """
    Calculate the information entropy of a set of observations, x.
    """


def calc_var_info_pqr(x: np.ndarray, y: np.ndarray, xy: np.ndarray,
                      validate: bool = True):
    """
    Calculate the variation of information for two partitions, X and Y,
    of the same set A. For more details and the source of the formula,
    see https://en.wikipedia.org/wiki/Variation_of_information.

    Parameters
    ----------
    x: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of X. All must be in (0, 1] and sum to 1.
    y: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of Y. All must be in (0, 1] and sum to 1.
    xy: ndarray
        A 2-dimensional array (len(p) x len(q)) of the fraction of the
        total elements that are in each pair of partitions from X and Y.
        All must be in (0, 1] and sum to 1.
    validate: bool = True
        Whether to validate the values of p, q, and r before calculating
        the variation of information.

    Returns
    -------
    float
        Variation of information
    """
    # Validate dimensions (always done even if validate is False because
    # the dimensions must be correct or else the entire calculation will
    # fail, and this step is much faster than those that require math).
    if x.ndim != 1:
        raise ValueError(f"Dimension of x must be 1, but got {x.ndim}")
    if x.size == 0:
        raise ValueError(f"x contained no elements")
    if y.ndim != 1:
        raise ValueError(f"Dimension of y must be 1, but got {y.ndim}")
    if y.size == 0:
        raise ValueError(f"y contained no elements")
    if xy.ndim != 2:
        raise ValueError(f"Dimension of xy must be 2, but got {xy.ndim}")
    if xy.shape != (x.size, y.size):
        raise ValueError(f"r must have dimensions of {x.size, y.size},"
                         f"but got {xy.shape}")
    if validate:
        # Validate bounds.
        if np.any(x <= 0.) or np.any(x > 1.):
            raise ValueError(f"All x must be in (0, 1], but got {x}")
        if np.any(y <= 0.) or np.any(y > 1.):
            raise ValueError(f"All y must be in (0, 1], but got {y}")
        if np.any(xy <= 0.) or np.any(xy > 1.):
            raise ValueError(f"All xy must be in (0, 1], but got {xy}")
        # Validate sums.
        if not np.isclose(x.sum(), 1.):
            raise ValueError(f"x must sum to 1, but got {x.sum()}")
        if not np.isclose(y.sum(), 1.):
            raise ValueError(f"y must sum to 1, but got {y.sum()}")
        if not np.isclose(xy.sum(), 1.):
            raise ValueError(f"xy must sum to 1, but got {xy.sum()}")
    # Compute the mutual information of X and Y.
    log_pq_grid = np.log(x)[:, np.newaxis] + np.log(y)[np.newaxis, :]
    return float(np.sum(xy * (log_pq_grid - 2. * np.log(xy))))


def calc_var_info_pair(resps1: pd.DataFrame, resps2: pd.DataFrame):
    """ Calculate the variation of information for a pair of clustering
    results. """
    # Find the number of reads and clusters.
    if not resps1.index.equals(resps2.index):
        raise ValueError(f"Read names of resps1 {resps1.index} "
                         f"and resps2 {resps2.index} differ")
    if (n_reads := resps1.index.size) == 0:
        # There is zero variation if there are zero reads.
        return 0.
    # Verify that resps1 and resps2 have the same reads (indexes).
    if not resps1.index.equals(resps2.index):
        raise ValueError("Read names of resps1 and clusts2 differ")
    # Verify that the probability that each read belongs to any
    # cluster equals 1.
    if not np.allclose(resps1.sum(axis=1), 1.):
        raise ValueError("Probabilities of reads in resps1 did not sum to 1")
    if not np.allclose(resps2.sum(axis=1), 1.):
        raise ValueError("Probabilities of reads in resps2 did not sum to 1")
    # For each run, compute the log proportion of reads in each cluster.
    logp1 = np.log(resps1.values.mean(axis=0))
    logp2 = np.log(resps2.values.mean(axis=0))
    # For each pair of clusters, compute the proportion of reads that
    # would be expected if the clusterings were independent.
    logp12_exp = logp1[:, np.newaxis] + logp2[np.newaxis, :]
    # For each pair of clusters, compute the proportion of reads that
    # was actually observed in both clusters.
    p12_obs = np.array([[np.vdot(resps1[c1], resps2[c2]) / n_reads
                         for c2 in resps2.columns]
                        for c1 in resps1.columns])
    logp12_obs = np.log(p12_obs)
    # Compute the variation of information between X and Y.
    return np.sum(p12_obs * (logp12_exp - 2 * logp12_obs))


def calc_mean_var_info(runs: list[EmClustering]):
    """ Calculate the expected variation of information among ≥ 2 runs
    of EM clustering. """
    # List every pair of EM runs.
    pairs = list(combinations(range(len(runs)), 2))
    if not pairs:
        # Variation of information defaults to 0 if no pairs exist.
        return 0.
    # FIXME: implement variation of information
    return 0.
    # Compute and cache the responsibilities of each run.
    resps = [run.get_resps() for run in runs]
    # Compute the variation of information for each pair of EM runs.
    vinfo = {(i, j): calc_var_info_pair(resps[i], resps[j]) for i, j in pairs}
    # Return the mean variation of information among pairs of EM runs.
    return sum(vinfo.values()) / len(vinfo)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
