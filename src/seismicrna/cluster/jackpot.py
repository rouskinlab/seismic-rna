from typing import Iterable

import numpy as np
from numba import jit

from .marginal import calc_marginal
from ..core.array import find_dims, get_length
from ..core.logs import logger
from ..core.random import stochastic_round
from ..core.unbias import (CLUSTERS,
                           POSITIONS,
                           READS,
                           calc_p_noclose_given_clust,
                           calc_p_nomut_window,
                           calc_p_noclose_given_ends,
                           calc_p_mut_given_span_noclose,
                           calc_p_clust_given_noclose,
                           calc_p_ends_given_clust_noclose,
                           calc_p_clust_given_ends_noclose)

SUM_EXP_PRECISION = 3
UNIQUE = "unique reads"

rng = np.random.default_rng()


def linearize_ends_matrix(p_ends: np.ndarray):
    """ Linearize an N x N matrix of end coordinate probabilities/counts
    into a length-N list of 5'/3' coordinates and their values. """
    end5s, end3s = np.nonzero(np.triu(p_ends))
    return end5s, end3s, p_ends[end5s, end3s]


@jit()
def _assign_clusters(clusters: np.ndarray,
                     p_clust_given_read: np.ndarray,
                     n_reads_per_clust: np.ndarray,
                     read_order: np.ndarray,
                     read_randomness: np.ndarray):
    """ Assign one cluster to each read.

    Parameters
    ----------
    clusters: np.ndarray
        1D (reads) integer array into which the cluster assigned to each
        read will be written.
    p_clust_given_read: np.ndarray
        2D (reads x clusters) float array of the probability that each
        read belongs to each cluster.
    n_reads_per_clust: np.ndarray
        1D (clusters) integer array of the number of reads to assign to
        each cluster.
    read_order: np.ndarray
        1D (reads) integer array of the random order in which to assign
        reads to clusters.
    read_randomness: np.ndarray
        1D (reads) float array of a random number (0 - 1) for each read,
        to be used for choosing the cluster.

    Additionally, this function makes the following assumptions; if any
    assumption is violated, this function can produce incorrect results
    or attempt to access invalid memory addresses, potentially causing
    a segmentation fault.
    - There must be at least 1 cluster.
    - Every element of `p_clust_given_read` is ≥ 0 and ≤ 1.
    - Every row of `p_clust_given_read` sums to 1.
    - All elements of `n_reads_per_clust` are non-negative integers.
    - `n_reads_per_clust` sums to the total number of reads (i.e. the
      lengths of `clusters`, `p_clust_given_read`, `read_order`, and
      `read_randomness`).
    - `read_order` contains every integer from 0 to the number of reads
      minus 1 exactly once (though the order can be random).
    - Every element of `read_randomness` is ≥ 0 and < 1.
    """
    n_clust, = n_reads_per_clust.shape
    # Iterate through the reads in random order to eliminate bias that
    # can arise from the order in which reads are assigned.
    for i in read_order:
        # Calculate the probability of selecting each cluster for this
        # read.
        p_clust_i = p_clust_given_read[i] * n_reads_per_clust
        p_clust_i /= p_clust_i.sum()
        # Select the cluster based on the probabilities.
        p_clust_i_sum = 0.
        for k in range(n_clust):
            p_clust_i_sum += p_clust_i[k]
            if read_randomness[i] < p_clust_i_sum:
                # Assign the read to the cluster.
                clusters[i] = k
                # Decrement the number of reads remaining in that cluster.
                n_reads_per_clust[k] -= 1
                break
    return clusters


def _sim_clusters(p_clust_given_read: np.ndarray):
    """ Simulate a cluster for each read. """
    n_reads, n_clust = p_clust_given_read.shape
    if p_clust_given_read.size > 0 and p_clust_given_read.min() <= 0.:
        raise ValueError("All p_clust_given_read must be > 0, but got "
                         f"{np.count_nonzero(p_clust_given_read <= 0.)} "
                         "probabilities < 0")
    if not np.allclose(p_clust_given_read.sum(axis=1), 1.):
        p_clust_given_read_sum = p_clust_given_read.sum(axis=1)
        not_close = np.logical_not(np.isclose(p_clust_given_read_sum, 1.))
        raise ValueError("All rows of p_clust_given_read must sum to 1, "
                         f"but got {p_clust_given_read_sum[not_close]}")
    # Initially assign all reads to cluster 0.
    clusters = np.zeros(n_reads, dtype=int)
    if n_clust == 1:
        # There is only one cluster, so all reads must be in it.
        return clusters
    # Choose the number of reads for each cluster, ensuring that the sum
    # equals the total number of reads.
    p_clust_given_read_sum = p_clust_given_read.sum(axis=0)
    while np.sum(n_reads_per_clust := stochastic_round(
            p_clust_given_read_sum,
            preserve_sum=True
    )) != n_reads:
        pass
    # Choose the cluster for each read.
    _assign_clusters(clusters,
                     p_clust_given_read,
                     n_reads_per_clust,
                     rng.permutation(n_reads),
                     rng.random(n_reads))
    return clusters


@jit()
def _calc_covered(covered: np.ndarray,
                  end5s: np.ndarray,
                  end3s: np.ndarray):
    """ Mark the bases that each read covers. """
    for i, (end5, end3) in enumerate(zip(end5s, end3s)):
        covered[i, end5: end3 + 1] = True


def _sim_muts(end5s: np.ndarray,
              end3s: np.ndarray,
              clusts: np.ndarray,
              p_mut_given_span_noclose: np.ndarray,
              min_mut_gap: int):
    """ Simulate mutations and write them into reads.

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end of each read.
    end3s: np.ndarray
        1D (reads) array of the 3' end of each read.
    clusts: np.ndarray
        1D (reads) array of the cluster to which each read belongs.
    p_mut_given_span_noclose: np.ndarray
        2D (positions x clusters) array of the probability that a base
        at each position is mutated given the read covers the position
        and no two mutations are too close.
    min_mut_gap: int
        Minimum number of positions between two mutations.

    Returns
    -------
    np.ndarray
        2D (reads x positions) boolean array of whether each base in
        each read is mutated.
    """
    dims = find_dims([(READS,),
                      (READS,),
                      (POSITIONS, CLUSTERS)],
                     [end5s,
                      end3s,
                      p_mut_given_span_noclose],
                     ["end5s",
                      "end3s",
                      "p_mut_given_span_noclose"],
                     nonzero=[CLUSTERS])
    n_reads = dims[READS]
    n_pos = dims[POSITIONS]
    n_clust = dims[CLUSTERS]
    # Initialize an empty matrix for the mutations.
    muts = np.zeros((n_reads, n_pos), dtype=bool)
    # Assign mutations to each cluster separately.
    for k in range(n_clust):
        # List all reads in cluster k.
        i_k = np.flatnonzero(clusts == k)
        # Determine which positions are covered by each read.
        covered = np.zeros((i_k.size, n_pos), dtype=bool)
        _calc_covered(covered, end5s[i_k], end3s[i_k])
        # Calculate the coverage and number of mutations per position.
        coverage = np.count_nonzero(covered, axis=0)
        num_muts = stochastic_round(coverage * p_mut_given_span_noclose[:, k])
        # Start filling in mutations at positions in order of increasing
        # number of non-mutated bases.
        muts_k = np.zeros_like(covered)
        for j in np.argsort(coverage - num_muts):
            # Determine which reads can be mutated at this position,
            # meaning they cover the position and do not have another
            # mutation within min_mut_gap positions.
            j_lo = max(j - min_mut_gap, 0)
            j_up = min(j + min_mut_gap + 1, n_pos)
            mutable = np.flatnonzero(np.logical_and(
                covered[:, j],
                np.count_nonzero(muts_k[:, j_lo: j_up], axis=1) == 0
            ))
            # Mutate a random subset of those reads.
            mutated = rng.choice(mutable,
                                 num_muts[j],
                                 replace=False,
                                 shuffle=False)
            muts_k[mutated, j] = True
        # After all positions have been finalized, copy them to muts.
        muts[i_k] = muts_k
    return muts


def _sim_reads(end5s: np.ndarray,
               end3s: np.ndarray,
               p_clust_given_ends_noclose: np.ndarray,
               p_mut_given_span_noclose: np.ndarray,
               min_mut_gap: int):
    find_dims([(READS,),
               (READS,),
               (POSITIONS, POSITIONS, CLUSTERS,),
               (POSITIONS, CLUSTERS)],
              [end5s,
               end3s,
               p_clust_given_ends_noclose,
               p_mut_given_span_noclose],
              ["end5s",
               "end3s",
               "p_clust_given_ends_noclose",
               "p_mut_given_span_noclose"],
              nonzero=[CLUSTERS])
    # Simulate the clusters and the mutations in each read.
    clusts = _sim_clusters(p_clust_given_ends_noclose[end5s, end3s])
    muts = _sim_muts(end5s,
                     end3s,
                     clusts,
                     p_mut_given_span_noclose,
                     min_mut_gap)
    # Merge the mutation data and 5'/3' ends into one array of reads.
    reads = np.hstack([muts,
                       end5s[:, np.newaxis],
                       end3s[:, np.newaxis]],
                      dtype=int)
    return reads, clusts


def _calc_resps(n_clusts: int,
                clusts: np.ndarray,
                uniq_inverse: np.ndarray,
                uniq_counts: np.ndarray):
    """ Calculate the cluster membership of each unique read.
    """
    dims = find_dims([(READS,), (READS,), (UNIQUE,)],
                     [clusts, uniq_inverse, uniq_counts],
                     ["clusts", "uniq_inverse", "uniq_counts"])
    resps = np.zeros((dims[UNIQUE], n_clusts))
    np.add.at(resps, (uniq_inverse, clusts), 1.)
    # Normalize the memberships so each row sums to 1.
    return resps / resps.sum(axis=1)[:, np.newaxis]


def _calc_obs_exp(reads: np.ndarray,
                  clusts: np.ndarray,
                  p_mut: np.ndarray,
                  p_ends: np.ndarray,
                  p_clust: np.ndarray,
                  min_mut_gap: int,
                  unmasked: np.ndarray):
    dims = find_dims([(READS, POSITIONS), (READS,)],
                     [reads, clusts],
                     ["reads", "clusts"])
    n_reads = dims[READS]
    # Count each unique read.
    uniq_reads, uniq_inverse, num_obs = np.unique(reads,
                                                  axis=0,
                                                  return_inverse=True,
                                                  return_counts=True)
    uniq_end5s = uniq_reads[:, -2]
    uniq_end3s = uniq_reads[:, -1]
    n_uniq_reads_observed, _ = uniq_reads.shape
    # Calculate the properties of the simulated reads.
    muts_per_pos = [np.flatnonzero(uniq_reads[:, j]) for j in unmasked]
    log_exp = np.log(n_reads) + calc_marginal(p_mut,
                                              p_ends,
                                              p_clust,
                                              uniq_end5s,
                                              uniq_end3s,
                                              unmasked,
                                              muts_per_pos,
                                              min_mut_gap)
    return num_obs, log_exp


def sim_obs_exp(end5s: np.ndarray,
                end3s: np.ndarray,
                p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                min_mut_gap: int,
                unmasked: np.ndarray):
    """ Simulate observed and expected counts. """
    find_dims([(READS,),
               (READS,),
               (POSITIONS, CLUSTERS),
               (POSITIONS, POSITIONS),
               (CLUSTERS,)],
              [end5s, end3s, p_mut, p_ends, p_clust],
              ["end5s", "end3s", "p_mut", "p_ends", "p_clust"])
    # Calculate the parameters of reads with no two mutations too close.
    p_nomut_window = calc_p_nomut_window(p_mut, min_mut_gap)
    p_noclose_given_ends = calc_p_noclose_given_ends(p_mut, p_nomut_window)
    p_mut_given_span_noclose = calc_p_mut_given_span_noclose(
        p_mut,
        p_ends,
        p_noclose_given_ends,
        p_nomut_window
    )
    p_noclose_given_clust = calc_p_noclose_given_clust(
        p_ends, p_noclose_given_ends
    )
    p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
        p_ends, p_noclose_given_ends
    )
    p_clust_given_noclose = calc_p_clust_given_noclose(
        p_clust, p_noclose_given_clust
    )
    p_clust_given_ends_noclose = calc_p_clust_given_ends_noclose(
        p_ends_given_clust_noclose, p_clust_given_noclose
    )
    while True:
        # Simulate the reads and the clusters to which they belong.
        reads, clusts = _sim_reads(end5s,
                                   end3s,
                                   p_clust_given_ends_noclose,
                                   p_mut_given_span_noclose,
                                   min_mut_gap)
        yield _calc_obs_exp(reads,
                            clusts,
                            p_mut,
                            p_ends,
                            p_clust,
                            min_mut_gap,
                            unmasked)


def calc_semi_g_anomaly(num_obs: int | np.ndarray,
                        log_exp: float | np.ndarray):
    """ Calculate each read's semi-G-anomaly, i.e. half its contribution
    to the G-test statistic. """
    return num_obs * (np.log(num_obs) - log_exp)


def calc_jackpot_score(semi_g_anomalies: np.ndarray, n_reads: int):
    """ Calculate the jackpotting score.

    The jackpotting score is defined as the average log of the ratio of
    observed to expected counts for a read, weighted by the count of
    each read:

    JS = sum{O_i * log(O_i / E_i)} / sum{O_i}
    where i is each unique read.

    This formula is based on that of the G-test statistic, which is
    identical except that the latter is multiplied by 2 and not divided
    by the number of observations.
    """
    if n_reads == 0:
        if semi_g_anomalies.size > 0:
            raise ValueError("If n_reads is 0, then semi_g_anomalies must have "
                             f"length 0, but got {semi_g_anomalies.size}")
        return 0.
    return semi_g_anomalies.sum() / n_reads


def calc_jackpot_score_ci(jackpot_scores: Iterable[float],
                          confidence_level: float):
    """ Calculate the confidence interval of the mean of an array of
    jackpotting scores. """
    # Ensure g_stats is a NumPy array of floats with no NaN/inf values.
    if not isinstance(jackpot_scores, (np.ndarray, list)):
        jackpot_scores = list(jackpot_scores)
    jackpot_scores = np.asarray_chkfinite(jackpot_scores, dtype=float)
    # Calculate the confidence interval of the mean, assuming that the
    # jackpotting scores are distributed normally.
    n = get_length(jackpot_scores, "jackpot_scores")
    if n > 1:
        # Import SciPy here instead of at the top of the module because
        # the latter would take so long as to slow down the start-up.
        from scipy.stats import t
        t_lo, t_up = t.interval(confidence_level, n - 1)
        mean = jackpot_scores.mean()
        std_err = np.sqrt(jackpot_scores.var() / n)
        ci_lo = mean + std_err * t_lo
        ci_up = mean + std_err * t_up
    else:
        ci_lo = np.nan
        ci_up = np.nan
    return ci_lo, ci_up


def calc_jackpot_quotient(real_jackpot_score: float,
                          null_jackpot_score: float):
    """ Calculate the jackpotting quotient.

    The jackpotting quotient indicates how much more overrepresented the
    average read is in the real dataset compared to the null dataset.

    Since the jackpotting score is the expected log-ratio of a read's
    observed to expected count, raising e to the power of it yields the
    observed-to-expected ratio, which measures jackpotting.

    JS = average{log(O / E)}
    exp(JS) = exp(average{log(O / E)})

    Then, the jackpotting quotient is the jackpotting score of the real
    dataset divided by that of the null dataset:

    JQ = exp(JS_real) / exp(JS_null)
       = exp(JS_real - JS_null)

    Parameters
    ----------
    real_jackpot_score: float
        Jackpotting score of the real dataset.
    null_jackpot_score: float
        Jackpotting score of the null model.

    Returns
    -------
    float
        Jackpotting quotient
    """
    return np.exp(real_jackpot_score - null_jackpot_score)


def bootstrap_jackpot_scores(uniq_end5s: np.ndarray,
                             uniq_end3s: np.ndarray,
                             counts_per_uniq: np.ndarray,
                             p_mut: np.ndarray,
                             p_ends: np.ndarray,
                             p_clust: np.ndarray,
                             min_mut_gap: int,
                             unmasked: np.ndarray,
                             real_jackpot_score: float,
                             confidence_level: float,
                             max_jackpot_quotient: float):
    """ Bootstrap jackpotting scores from the null model.

    Parameters
    ----------
    uniq_end5s: np.ndarray
        1D (unique reads) array of the 5' ends (0-indexed) of the unique
        reads.
    uniq_end3s: np.ndarray
        1D (unique reads) array of the 3' ends (0-indexed) of the unique
        reads.
    counts_per_uniq: np.ndarray
        1D (unique reads) array of the number of times each unique read
        was observed.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    p_ends: np.ndarray
        2D (positions x positions) array of the probability distribution
        of 5'/3' end coordinates.
    p_clust: np.ndarray
        1D (clusters) probability that a read comes from each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    unmasked: np.ndarray
        Positions (0-indexed) of p_mut that are not masked.
    real_jackpot_score: float
        Jackpotting score of the real dataset.
    confidence_level: float
        Confidence level for computing a confidence interval of the
        jackpotting quotient.
    max_jackpot_quotient: float
        Stop bootstrapping once the confidence interval lies entirely
        above or entirely below this threshold.

    Returns
    -------
    np.ndarray
        Array of normalized G-test statistics for the null model.
    """
    find_dims([(UNIQUE,),
               (UNIQUE,),
               (UNIQUE,),
               (POSITIONS, CLUSTERS),
               (POSITIONS, POSITIONS),
               (CLUSTERS,)],
              [uniq_end5s,
               uniq_end3s,
               counts_per_uniq,
               p_mut,
               p_ends,
               p_clust],
              ["uniq_end5s",
               "uniq_end3s",
               "counts_per_uniq",
               "p_mut",
               "p_ends",
               "p_clust"])
    n_reads = counts_per_uniq.sum()
    logger.routine(f"Began boostrapping null jackpotting scores for a dataset "
                   f"with {n_reads} reads and a real jackpotting score of "
                   f"{real_jackpot_score}")
    # Simulate observed and expected read counts.
    end5s = np.repeat(uniq_end5s, counts_per_uniq)
    end3s = np.repeat(uniq_end3s, counts_per_uniq)
    null_jackpotting_scores = list()
    for num_obs, log_exp in sim_obs_exp(end5s,
                                        end3s,
                                        p_mut,
                                        p_ends,
                                        p_clust,
                                        min_mut_gap,
                                        unmasked):
        # Calculate this null model's jackpotting score.
        null_g_anomalies = calc_semi_g_anomaly(num_obs, log_exp)
        null_jackpotting_score = calc_jackpot_score(null_g_anomalies,
                                                    n_reads)
        null_jackpotting_scores.append(null_jackpotting_score)
        logger.detail(f"Null jackpotting score: {null_jackpotting_score}")
        # Calculate a confidence interval for the mean jackpotting score
        # of the null models simulated so far.
        js_ci_lo, js_ci_up = calc_jackpot_score_ci(null_jackpotting_scores,
                                                   confidence_level)
        # Determine a confidence interval for the jackpotting quotient.
        # Since the jackpotting quotient depends inversely on the null
        # jackpotting score, the lower and upper bounds are swapped.
        jq_ci_lo = calc_jackpot_quotient(real_jackpot_score, js_ci_up)
        jq_ci_up = calc_jackpot_quotient(real_jackpot_score, js_ci_lo)
        if not np.isnan(jq_ci_lo) and not np.isnan(jq_ci_up):
            logger.detail(f"{confidence_level * 100.} % confidence interval "
                          f"for jackpotting quotient: {jq_ci_lo} - {jq_ci_up}")
        # Stop when the confidence interval lies entirely below or above
        # max_jackpot_quotient, so it's clear whether the jackpotting
        # quotient is less or greater than max_jackpot_quotient.
        # Avoid "if not jq_ci_lo <= max_jackpot_quotient <= jq_ci_up"
        # because this expression will evaluate to True after the first
        # iteration, when jq_ci_lo and jq_ci_up will both be NaN.
        if jq_ci_lo > max_jackpot_quotient or max_jackpot_quotient > jq_ci_up:
            logger.routine("Ended boostrapping null jackpotting scores")
            return np.array(null_jackpotting_scores)
