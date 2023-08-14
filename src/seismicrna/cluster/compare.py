"""
Cluster Comparison Module
========================================================================
Auth: Matty

Collect and compare the results from independent runs of EM clustering.
"""

from __future__ import annotations

from itertools import combinations, islice
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .em import EmClustering
from .names import EXP_NAME, OBS_NAME, ORD_NAME


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
    try:
        get_common_order(runs)
    except ValueError:
        raise ValueError("Cannot sort replicate runs of multiple orders")
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
    value = ord_runs[1].best.__getattribute__(attr)
    # Verify that the best run from every other order has an equal value
    # of that attribute.
    if any(runs.best.__getattribute__(attr) != value
           for order, runs in ord_runs.items() if order != 1):
        raise ValueError(f"Found more than 1 value for attribute '{attr}' "
                         f"among EM clustering runs {ord_runs}")
    return value


def find_best_order(ord_runs: dict[int, RunOrderResults]) -> int:
    """ Find the number of clusters with the best (smallest) BIC. """
    return sorted(ord_runs.items(), key=lambda runs: runs[1].best.bic)[0][0]


def format_exp_count_col(order: int):
    return f"{EXP_NAME}, {ORD_NAME} {order}"


def get_log_exp_obs_counts(ord_runs: dict[int, RunOrderResults]):
    """ Get the expected and observed log counts of each bit vector. """
    # Retrieve the unique bit vectors from the clusters.
    uniq_muts = get_common_best_run_attr(ord_runs, "muts")
    # Compute the observed log counts of the bit vectors.
    log_obs = pd.Series(np.log(uniq_muts.counts),
                        index=uniq_muts.get_uniq_names())
    # For each order of clustering, compute the expected log counts.
    log_exp = ((format_exp_count_col(order),
                pd.Series(runs.best.log_exp_counts, index=log_obs.index))
               for order, runs in ord_runs.items())
    # Assemble all log counts into one DataFrame.
    log_counts = pd.DataFrame.from_dict({OBS_NAME: log_obs, **dict(log_exp)})
    # Sort the data by expected count at order 1, then round it.
    return log_counts.sort_values(by=[format_exp_count_col(order=1)],
                                  ascending=False).round(EXP_COUNT_PRECISION)


def calc_var_info_pqr(p: np.ndarray, q: np.ndarray, r: np.ndarray,
                      validate: bool = True):
    """
    Calculate the variation of information for two partitions, X and Y,
    of the same set A. For more details and the source of the formula,
    see https://en.wikipedia.org/wiki/Variation_of_information.

    Parameters
    ----------
    p: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of X. All must be in (0, 1] and sum to 1.
    q: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of Y. All must be in (0, 1] and sum to 1.
    r: ndarray
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
    if p.ndim != 1:
        raise ValueError(f"Dimension of p must be 1, but got {p.ndim}")
    if p.size == 0:
        raise ValueError(f"p contained no elements")
    if q.ndim != 1:
        raise ValueError(f"Dimension of q must be 1, but got {q.ndim}")
    if q.size == 0:
        raise ValueError(f"q contained no elements")
    if r.ndim != 2:
        raise ValueError(f"Dimension of r must be 2, but got {r.ndim}")
    if r.shape != (p.size, q.size):
        raise ValueError("r must have dimensions of p.size x q.size "
                         f"{p.size, q.size}, but got {r.shape}")
    if validate:
        # Validate bounds.
        if np.any(p <= 0.) or np.any(p > 1.):
            raise ValueError(f"All values in p must be in (0, 1], but got {p}")
        if np.any(q <= 0.) or np.any(q > 1.):
            raise ValueError(f"All values in q must be in (0, 1], but got {q}")
        if np.any(r <= 0.) or np.any(r > 1.):
            raise ValueError(f"All values in r must be in (0, 1], but got {r}")
        # Validate sums.
        if not np.isclose(p.sum(), 1.):
            raise ValueError(f"p must sum to 1, but got {p.sum()}")
        if not np.isclose(q.sum(), 1.):
            raise ValueError(f"q must sum to 1, but got {q.sum()}")
        if not np.isclose(r.sum(), 1.):
            raise ValueError(f"r must sum to 1, but got {r.sum()}")
    # Compute the variation of information.
    log_pq_grid = np.log(p)[:, np.newaxis] + np.log(q)[np.newaxis, :]
    return float(np.sum(r * (log_pq_grid - 2. * np.log(r))))


def calc_var_info_run_pair(clusts1: pd.DataFrame, clusts2: pd.DataFrame,
                           validate: bool = True):
    """ Calculate the variation of information between two EM runs. """
    # Find the number of reads and clusters.
    n_reads, n_clusts = clusts1.shape
    if clusts2.shape != (n_reads, n_clusts):
        raise ValueError(f"Dimensions of clusts1 {clusts1.shape} "
                         f"and clusts2 {clusts2.shape} differ")
    if n_clusts <= 0:
        raise ValueError(f"Number of clusters must be ≥ 1, but got {n_clusts}")
    if n_reads == 0:
        # There is zero variation if there are zero reads.
        return 0.
    if validate:
        # Verify that run1 and run2 have the same reads (indexes) and
        # clusters (columns).
        if not clusts1.index.equals(clusts2.index):
            raise ValueError("Read names of clusts1 and clusts2 differ")
        if not clusts1.columns.equals(clusts2.columns):
            raise ValueError("Clusters of clusts1 and clusts2 differ")
        # Verify that the probability that each read belongs to any
        # cluster equals 1.
        if not np.allclose(clusts1.sum(axis=1), 1.):
            raise ValueError(
                "Probabilities of reads in clusts1 did not all sum to 1")
        if not np.allclose(clusts2.sum(axis=1), 1.):
            raise ValueError(
                "Probabilities of reads in clusts2 did not all sum to 1")
    # For each run, compute the proportion of reads in each cluster.
    props1 = clusts1.mean(axis=0).values
    props2 = clusts2.mean(axis=0).values
    # For each pair of clusters, compute the proportion of reads in both
    # clusters.
    props12 = np.array([[np.vdot(clusts1[c1], clusts2[c2]) / n_reads
                         for c2 in clusts2.columns]
                        for c1 in clusts1.columns])
    # Compute the variation of information.
    return calc_var_info_pqr(props1, props2, props12, validate=validate)


def calc_mean_var_info(runs: list[EmClustering]):
    """ Calculate the expected variation of information among ≥ 2 runs
    of EM clustering. """
    # List every pair of EM runs.
    pairs = list(combinations(range(len(runs)), 2))
    if not pairs:
        # Variation of information defaults to 0 if no pairs exist.
        return 0.
    # Compute and cache the responsibilities of each run.
    resps = [run.output_resps() for run in runs]
    # Find the mean variation of information among pairs of EM runs.
    sum_var_info = sum(calc_var_info_run_pair(resps[i1], resps[i1])
                       for i1, i2 in pairs)
    return sum_var_info / len(pairs)


def iter_all_likelihoods(p: np.ndarray):
    """ Yield every bit vector likelihood from greatest to least. """
    if p.ndim != 1:
        raise ValueError(f"p must have 1 dimension, but got {p.ndim}")
    # If no bits remain, then just return the base value.
    if p.size == 0:
        return iter([0.])
    if np.min(p) < 0. or np.max(p) > 1.:
        raise ValueError(f"p must all be in [0, 1], but got {p}")
    # Remove any bits that are certain to be 0 or 1.
    p = p[np.logical_not(np.logical_or(p == 0., p == 1.))]
    # If no bits remain, then just return the base value.
    if p.size == 0:
        return iter([0.])
    # Compute the log probabilities (p) and anti-probabilities (q).
    log_p = np.log(p)
    log_q = np.log(1. - p)
    # Sum the larger of the two as the "base" weight.
    log_base = float(np.sum(np.maximum(log_p, log_q)))
    # Compute the difference in logs and sort from least to greatest.
    log_diffs = np.abs(log_p - log_q)
    log_diffs.sort()

    class LikelihoodSeries(object):

        def __init__(self, weight: float, inner: Callable[[], Iterable]):
            self._weight = weight
            self._inner = inner
            self._values = list()

        def _cache(self, value: float):
            """ Store a value in the values list. """
            self._values.append(value)
            return value

        def iter(self):
            """ Iterate over the values. """
            # Check if the values were already computed.
            if self._values:
                # If the values were already computed, then yield them.
                yield from self._values
            else:
                # Otherwise, make iterators for the inner series and the
                # weighted series.
                inner = iter(self._inner())
                weighted = iter(ival - self._weight for ival in self._inner())
                # Get the first element from each iterator. Both should
                # have > 0 items and should not raise StopIteration.
                ival = next(inner)
                wval = next(weighted)
                # Yield elements until both iterators are exhausted.
                while True:
                    # Yield the larger value and advance its iterator.
                    if ival is not None and ival >= wval:
                        # Yield the inner value.
                        yield self._cache(ival)
                        # Get the next inner value, or None if empty.
                        ival = next(inner, None)
                    else:
                        # Yield the weighted value.
                        yield self._cache(wval)
                        try:
                            # Get the next weighted value, if any.
                            wval = next(weighted)
                        except StopIteration:
                            # There are no more weighted values.
                            break
                # Confirm that there are no more inner values.
                try:
                    next(inner)
                except StopIteration:
                    # There are no more inner values: success.
                    return
                # There were more inner values: failure.
                raise ValueError(
                    "Weighted values were exhausted before inner values")

    # Make series of log likelihoods.
    series = LikelihoodSeries(log_diffs[0], lambda: [0.])
    for log_diff in log_diffs[1:]:
        series = LikelihoodSeries(log_diff, series.iter)

    # Return an iterator of log likelihood values adjusted by the base.
    return iter(value + log_base for value in series.iter())


def max_n_likelihoods(p: np.ndarray, n: int | None = None):
    """ Return the top n bit vector likelihoods. """
    return np.array(list(islice(iter_all_likelihoods(p), n)))


def cum_log_exp_obs_counts():
    """ Get the cumulative expected and observed log counts per bit
    vector. """
