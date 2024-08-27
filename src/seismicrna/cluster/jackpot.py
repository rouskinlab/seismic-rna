from logging import getLogger
from typing import Iterable

import numpy as np
from numba import jit

from .marginal import calc_marginal
from ..core.array import find_dims, get_length
from ..core.unbias import (CLUSTERS,
                           POSITIONS,
                           READS,
                           calc_p_noclose_given_clust,
                           calc_p_noclose_given_ends_auto,
                           calc_p_clust_given_noclose,
                           calc_p_noclose,
                           calc_p_ends_given_clust_noclose,
                           calc_p_clust_given_ends_noclose,
                           calc_params,
                           calc_params_observed)

SUM_EXP_PRECISION = 3
UNIQUE = "unique reads"

rng = np.random.default_rng()

logger = getLogger(__name__)


def linearize_ends_matrix(p_ends: np.ndarray):
    """ Linearize an N x N matrix of end coordinate probabilities/counts
    into a length-N list of 5'/3' coordinates and their values. """
    end5s, end3s = np.nonzero(np.triu(p_ends))
    return end5s, end3s, p_ends[end5s, end3s]


@jit()
def _rand_muts(p_mut_read: np.ndarray):
    """ Simulate random mutations for one read. """
    return np.flatnonzero(np.random.random(p_mut_read.size) < p_mut_read)


@jit()
def _sim_muts(reads: np.ndarray,
              clusts: np.ndarray,
              p_mut: np.ndarray,
              min_mut_gap: int):
    """ Simulate mutations and write them into reads.

    Parameters
    ----------
    reads: np.ndarray
        2D (reads x positions + 2) array into which the simulated reads
        will be written.
    clusts: np.ndarray
        1D (reads) array of the cluster to which each read belongs.
    p_mut: np.ndarray
        2D (positions x clusters) array of the probability that each
        position is mutated for each cluster.
    min_mut_gap: int
        Minimum number of positions between two mutations.
    """
    for read, k in zip(reads, clusts):
        # Determine the 5'/3' ends.
        end5 = read[-2]
        end3 = read[-1]
        # Choose which positions are mutated.
        p_mut_read = p_mut[end5: end3 + 1, k]
        muts = _rand_muts(p_mut_read)
        if min_mut_gap > 0:
            # Continue choosing which positions are mutated until no two
            # mutations are too close.
            while muts.size > 1 and np.diff(muts).min() <= min_mut_gap:
                muts = _rand_muts(p_mut_read)
        # Write the mutated positions into the array of reads.
        read[muts + end5] = 1


def _sim_clusters(end5s: np.ndarray,
                  end3s: np.ndarray,
                  p_clust_given_ends_noclose: np.ndarray):
    """ Simulate a cluster for each read. """
    dims = find_dims([(READS,),
                      (READS,),
                      (POSITIONS, POSITIONS, CLUSTERS,)],
                     [end5s, end3s, p_clust_given_ends_noclose],
                     ["end5s", "end3s", "p_clust_given_ends_noclose"])
    # For each read, calculate the cumulative probability that the read
    # belongs to each cluster.
    cum_p_clust_per_read = np.cumsum(p_clust_given_ends_noclose[end5s, end3s],
                                     axis=1)
    # For each read, generate a random number from 0 to 1 and pick the
    # cluster in whose window of probability that random number lies.
    return np.count_nonzero(rng.random((dims[READS], 1)) > cum_p_clust_per_read,
                            axis=1)


def _sim_reads(end5s: np.ndarray,
               end3s: np.ndarray,
               p_clust_given_ends_noclose: np.ndarray,
               p_mut: np.ndarray,
               min_mut_gap: int):
    dims = find_dims([(READS,),
                      (READS,),
                      (POSITIONS, POSITIONS, CLUSTERS,),
                      (POSITIONS, CLUSTERS)],
                     [end5s, end3s, p_clust_given_ends_noclose, p_mut],
                     ["end5s", "end3s", "p_clust_given_ends_noclose", "p_mut"])
    # Assign each read to exactly one cluster.
    clusts = _sim_clusters(end5s, end3s, p_clust_given_ends_noclose)
    # Initialize an empty matrix for the reads.
    reads = np.zeros((dims[READS], dims[POSITIONS] + 2), dtype=int)
    # Write the 5'/3' ends for the reads into the last two columns.
    reads[:, -2] = end5s
    reads[:, -1] = end3s
    # Generate random mutations and write them into the matrix of reads.
    _sim_muts(reads, clusts, p_mut, min_mut_gap)
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
                  n_clusts: int,
                  min_mut_gap: int,
                  unmasked: np.ndarray):
    dims = find_dims([(READS, POSITIONS), (READS,)],
                     [reads, clusts],
                     ["reads", "clusts"])
    n_pos = dims[POSITIONS] - 2
    n_reads, _ = reads.shape
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
    resps = _calc_resps(n_clusts, clusts, uniq_inverse, num_obs)
    (p_mut_observed,
     p_ends_observed,
     p_clust_observed) = calc_params_observed(n_pos,
                                              unmasked,
                                              muts_per_pos,
                                              uniq_end5s,
                                              uniq_end3s,
                                              num_obs,
                                              resps)
    # Correct the bias in the simulated parameters.
    p_mut_sim, p_ends_sim, p_clust_sim = calc_params(p_mut_observed,
                                                     p_ends_observed,
                                                     p_clust_observed,
                                                     min_mut_gap)
    # Calculate the expected number of each read in the simulated data
    # from its own parameters (p_mut_sim, p_ends_sim, p_clust_sim),
    # since the real data are also
    log_exp = np.log(n_reads) + calc_marginal(p_mut_sim,
                                              p_ends_sim,
                                              p_clust_sim,
                                              uniq_reads[:, -2],
                                              uniq_reads[:, -1],
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
    dims = find_dims([(READS,),
                      (READS,),
                      (POSITIONS, CLUSTERS),
                      (POSITIONS, POSITIONS),
                      (CLUSTERS,)],
                     [end5s, end3s, p_mut, p_ends, p_clust],
                     ["end5s", "end3s", "p_mut", "p_ends", "p_clust"])
    n_clusts = dims[CLUSTERS]
    # Calculate the parameters of reads with no two mutations too close.
    p_noclose_given_ends = calc_p_noclose_given_ends_auto(
        p_mut, min_mut_gap
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
    # Inflate the number of reads to simulate because a proportion of
    # them will have two mutations too close.
    p_noclose = calc_p_noclose(p_clust, p_noclose_given_clust)
    if p_noclose == 0.:
        raise ValueError("It is impossible for any read to be observed")
    while True:
        # Simulate the reads and the clusters to which they belong.
        reads, clusts = _sim_reads(end5s,
                                   end3s,
                                   p_clust_given_ends_noclose,
                                   p_mut,
                                   min_mut_gap)
        yield _calc_obs_exp(reads, clusts, n_clusts, min_mut_gap, unmasked)


def calc_semi_g_anomaly(num_obs: float | np.ndarray,
                        log_exp: float | np.ndarray):
    """ Calculate each item's semi-G-anomaly, i.e. half its contribution
    to the G-test statistic. """
    return num_obs * (np.log(num_obs) - log_exp)


def calc_jackpot_score(semi_g_anomalies: np.ndarray, n_reads: int):
    """ Calculate the jackpotting score.

    The jackpotting score is defined as the average log of the ratio of
    observed to expected counts for a read, weighted by the count of
    each read:

    JS = sum{O_i * log(O_i / E_i)} / sum{O_i}
    where i is each unique read.

    This formula was based on that of the G-test statistic, which is
    identical except that the latter is multiplied by 2.
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
    if not isinstance(jackpot_scores, np.ndarray):
        if not isinstance(jackpot_scores, list):
            jackpot_scores = list(jackpot_scores)
        jackpot_scores = np.array(jackpot_scores)
    jackpot_scores = np.asarray_chkfinite(jackpot_scores, dtype=float)
    # Calculate the confidence interval of the mean, assuming that the
    # jackpotting scores are distributed normally.
    n = get_length(jackpot_scores, "jackpot_scores")
    if n > 1:
        mean = jackpot_scores.mean()
        std_dev = jackpot_scores.std()
        std_err = std_dev / np.sqrt(n)
        # Import SciPy here instead of at the top of the module because
        # the latter would take so long as to slow down the start-up.
        from scipy.stats import t
        t_lo, t_up = t.interval(confidence_level, n - 1)
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
    logger.info(f"Began boostrapping null jackpotting scores for a dataset "
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
        logger.debug(f"Null jackpotting score: {null_jackpotting_score}")
        # Calculate a confidence interval for the mean jackpotting score
        # of the null models simulated so far.
        njs_ci_lo, njs_ci_up = calc_jackpot_score_ci(null_jackpotting_scores,
                                                     confidence_level)
        # Determine a confidence interval for the jackpotting quotient.
        # Since the jackpotting quotient depends inversely on the null
        # jackpotting score, the lower and upper bounds are swapped.
        jq_ci_lo = calc_jackpot_quotient(real_jackpot_score, njs_ci_up)
        jq_ci_up = calc_jackpot_quotient(real_jackpot_score, njs_ci_lo)
        if not np.isnan(jq_ci_lo) and not np.isnan(jq_ci_up):
            logger.debug(f"{confidence_level * 100.} % confidence interval "
                         f"for jackpotting quotient: {jq_ci_lo} - {jq_ci_up}")
        # Stop when the confidence interval lies entirely below or above
        # max_jackpot_quotient, so it's clear whether the jackpotting
        # quotient is less or greater than max_jackpot_quotient.
        # Avoid "if not jq_ci_lo <= max_jackpot_quotient <= jq_ci_up"
        # because this expression will evaluate to True after the first
        # iteration, when jq_ci_lo and jq_ci_up will both be NaN.
        if jq_ci_lo > max_jackpot_quotient or max_jackpot_quotient > jq_ci_up:
            logger.info("Ended boostrapping null jackpotting scores")
            return np.array(null_jackpotting_scores)
