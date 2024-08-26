from logging import getLogger
from typing import Iterable

import numpy as np

from .marginal import calc_marginal
from ..core.array import (ensure_same_length,
                          find_dims,
                          get_length)
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


# @jit()
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
        muts = np.flatnonzero(rng.random(p_mut_read.size) < p_mut_read)
        if min_mut_gap > 0:
            # Continue choosing which positions are mutated until no two
            # mutations are too close.
            while muts.size > 1 and np.diff(muts).min() <= min_mut_gap:
                muts = np.flatnonzero(rng.random(p_mut_read.size) < p_mut_read)
        # Write the mutated positions into the array of reads.
        read[muts + end5] = 1


def _find_reads_no_close(muts: np.ndarray, min_mut_gap: int):
    """ Find reads with no two mutations closer than `min_mut_gap`. """
    n_reads, n_pos = muts.shape
    # It is impossible for two mutations to be separated by more than
    # (n_pos - 1) positions.
    min_mut_gap = min(min_mut_gap, n_pos - 1)
    if min_mut_gap <= 0:
        # It is impossible for two mutations to be too close.
        return np.ones(n_reads, dtype=bool)
    # Using a sliding window of (min_mut_gap + 1) positions, label which
    # reads have no more than 1 mutation within any window.
    window_size = min_mut_gap + 1
    muts_cumsum = np.hstack([np.zeros((n_reads, window_size), dtype=int),
                             np.cumsum(muts, axis=1)])
    window_sum = muts_cumsum[:, window_size:] - muts_cumsum[:, :-window_size]
    return window_sum.max(axis=1) < 2


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


def calc_half_g_anomaly(num_obs: float | np.ndarray,
                        log_exp: float | np.ndarray):
    """ Calculate half of an item's G anomaly, i.e. its contribution to
    the overall G-test statistic. """
    return num_obs * (np.log(num_obs) - log_exp)


def calc_jackpot_g_stat(num_obs: np.ndarray,
                        log_exp: np.ndarray,
                        min_exp: float = 0.):
    """ Calculate the G-test statistic for jackpotting.

    Parameters
    ----------
    num_obs: np.ndarray
        Number of times each read was observed (1D, int).
    log_exp: np.ndarray
        Log of the number of times each read was expected to be observed
        (1D, float, same shape as `num_obs`).
    min_exp: float
        Aggregate all reads whose expected count is less than `min_exp`
        into a single category. Doing so corrects two violations of the
        assumptions of G-tests: that all reads' expected counts are at
        least approximately 5 (the rarest reads' expected counts are
        very close to 0) and that the observed and expected counts have
        equal sums (the sum of the expected counts can be less because
        it omits the expected counts of reads that were not observed).

    Returns
    -------
    tuple[float, int]
        - G-test statistic, normalized by total number of reads observed
        - Number of degrees of freedom
    """
    if min_exp < 0.:
        raise ValueError(f"min_exp must be ≥ 0, but got {min_exp}")
    num_uniq = ensure_same_length(num_obs, log_exp, "observed", "expected")
    if num_uniq > 0:
        min_obs = num_obs.min()
        if min_obs < 0:
            raise ValueError(
                f"All num_obs must be ≥ 0, but got {num_obs[num_obs < 0]}"
            )
    # Total number of reads.
    total_obs = round(num_obs.sum())
    # Expected count of each unique read.
    num_exp = np.exp(log_exp)
    # Total expected count of all reads that were observed 0 times.
    total_exp = round(num_exp.sum(), SUM_EXP_PRECISION)
    obs_0_exp = total_obs - total_exp
    # Degree of freedom for all reads that were observed 0 times.
    if obs_0_exp > 0.:
        df_obs_0 = 1
    elif obs_0_exp == 0.:
        df_obs_0 = 0
    else:
        raise ValueError(f"Total observed reads ({total_obs}) is less than "
                         f"total expected reads ({total_exp})")
    if min_exp > 0.:
        # Calculate the G-test statistic of reads whose expected counts
        # are at least min_exp.
        is_at_least_min_exp = num_exp >= min_exp
        at_least_min_exp = np.flatnonzero(is_at_least_min_exp)
        g_stat = 2. * calc_half_g_anomaly(num_obs[at_least_min_exp],
                                          log_exp[at_least_min_exp]).sum()
        # Degrees of freedom for those reads.
        df_at_least_min_exp = at_least_min_exp.size - 1
        # Add reads with expected counts less than min_exp to the G-test.
        # Assume that every read that was observed 0 times has an expected
        # count less than min_exp; this assumption should be valid as long
        # as min_exp is at least 5 (the usual minimum for a G- or χ²-test).
        # Making this assumption alleviates the need to compute the expected
        # count of every unseen read, which is computationally infeasible.
        # Thus, simply add obs_0_exp to the total of the expected count of
        # reads with expected counts less than num_exp.
        less_than_min_exp = np.flatnonzero(np.logical_not(is_at_least_min_exp))
        num_obs_less_than_min_exp = round(num_obs[less_than_min_exp].sum())
        if num_obs_less_than_min_exp > 0:
            num_exp_less_than_min_exp = (num_exp[less_than_min_exp].sum()
                                         + obs_0_exp)
            if num_exp_less_than_min_exp == 0.:
                raise ValueError("Expected number of reads with expected "
                                 f"counts < {min_exp} cannot be 0 if the "
                                 "observed number of such reads is "
                                 f"{num_obs_less_than_min_exp}")
            g_stat += 2. * calc_half_g_anomaly(
                num_obs_less_than_min_exp,
                np.log(num_exp_less_than_min_exp)
            )
            # Degree of freedom for all reads with expected counts less
            # than min_exp.
            df_less_than_min_exp = 1
        else:
            df_less_than_min_exp = 0
    else:
        # Calculate the G-test statistic of all reads.
        g_stat = 2. * calc_half_g_anomaly(num_obs, log_exp).sum()
        # Degrees of freedom for those reads.
        df_at_least_min_exp = num_exp.size - 1
        # Degree of freedom for all reads with expected counts less
        # than min_exp.
        df_less_than_min_exp = 0
    # Calculate total degrees of freedom.
    df = max(df_at_least_min_exp + df_less_than_min_exp + df_obs_0, 0)
    if df == 0 and g_stat != 0.:
        raise ValueError("G-test statistic must be 0 if the degrees of "
                         f"freedom equals 0, but got {g_stat}")
    return g_stat, df


def normalize_g_stat(g_stat: float, df: float):
    """ Normalize the G-test statistic.

    The normalized G-test statistic is defined as the expected value of
    the log of the ratio of observed to expected counts for a read
    (or, in other words, the weighted average of those log ratios):

    normG = sum{O_i * log(O_i / E_i)} / sum{O_i}
    where i is each unique read.

    As sum{O_i * log(O_i / E_i)} equals half of the G-test statistic (G)
    and sum{O_i} equals the total number of reads observed (n), normG
    simply equals G / (2n).
    """
    return g_stat / df


def calc_norm_g_stat_ci(norm_g_stats: Iterable[float],
                        confidence_level: float):
    """ Calculate the confidence interval of the mean of an array of
    normalized G-test statistics. """
    # Ensure g_stats is a NumPy array of floats with no NaN/inf values.
    if not isinstance(norm_g_stats, np.ndarray):
        if not isinstance(norm_g_stats, list):
            norm_g_stats = list(norm_g_stats)
        norm_g_stats = np.array(norm_g_stats)
    norm_g_stats = np.asarray_chkfinite(norm_g_stats,
                                        dtype=float)
    # Calculate the confidence interval of the mean, assuming that the
    # G-test statistics are distributed normally.
    n = get_length(norm_g_stats, "mean_obs_exp_log_ratios")
    if n > 1:
        mean = norm_g_stats.mean()
        std_dev = norm_g_stats.std()
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


def calc_jackpot_index(real_norm_g_stat: float,
                       null_norm_g_stat: float):
    """ Calculate the jackpotting index.

    The jackpotting index indicates by what factor the average read is
    overrepresented in the real dataset compared to the null dataset.

    Since the normalized G-stat is the expected log-ratio of a read's
    observed to expected count, raising e to the power of it yields the
    observed-to-expected ratio, which measures jackpotting:

    Jratio = O / E
           = exp{log(O / E)}
           = exp(normG)

    Then, the jackpotting index is the quotient of the jackpotting ratio
    of the real dataset versus the null dataset:

    Jindex = Jratio_real / Jratio_null
           = exp(normG_real) / exp(normG_null)
           = exp(normG_real - normG_null)

    Parameters
    ----------
    real_norm_g_stat: float
        Normalized G-test statistic of the real dataset.
    null_norm_g_stat: float
        Normalized G-test statistic(s) based on the null model.

    Returns
    -------
    float
        Jackpotting index
    """
    return np.exp(real_norm_g_stat - null_norm_g_stat)


def bootstrap_norm_g_stats(uniq_end5s: np.ndarray,
                           uniq_end3s: np.ndarray,
                           counts_per_uniq: np.ndarray,
                           p_mut: np.ndarray,
                           p_ends: np.ndarray,
                           p_clust: np.ndarray,
                           min_mut_gap: int,
                           unmasked: np.ndarray,
                           real_norm_g_stat: float,
                           confidence_level: float,
                           max_jackpot_index: float):
    """ Bootstrap normalized G-test statistics from the null model.

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
    real_norm_g_stat: float
        Normalized G-test statistic of the real dataset.
    confidence_level: float
        Level for computing a confidence interval of the mean normalized
        G-test statistic and of the jackpotting index.
    max_jackpot_index: float
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
    logger.info(f"Began boostrapping null normalized G-test statistics for a "
                f"dataset with {n_reads} reads and a real normalized G-test "
                f"statistic of {real_norm_g_stat}")
    # Simulate observed and expected read counts.
    end5s = np.repeat(uniq_end5s, counts_per_uniq)
    end3s = np.repeat(uniq_end3s, counts_per_uniq)
    null_norm_g_stats = list()
    for num_obs, log_exp in sim_obs_exp(end5s,
                                        end3s,
                                        p_mut,
                                        p_ends,
                                        p_clust,
                                        min_mut_gap,
                                        unmasked):
        # Calculate and normalize this null model G-test statistic.
        null_g_stat, _ = calc_jackpot_g_stat(num_obs, log_exp)
        null_norm_g_stat = normalize_g_stat(null_g_stat, n_reads)
        null_norm_g_stats.append(null_norm_g_stat)
        logger.debug(f"Null normalized G-test statistic: {null_norm_g_stat}")
        # Determine a confidence interval for the mean of all normalized
        # G-test statistics simulated so far.
        g_ci_lo, g_ci_up = calc_norm_g_stat_ci(null_norm_g_stats,
                                               confidence_level)
        # Determine a confidence interval for the jackpotting index.
        # Since the jackpotting index is inversely related to the null
        # G-test statistic, the lower and upper bounds are swapped.
        ji_ci_lo = calc_jackpot_index(real_norm_g_stat, g_ci_up)
        ji_ci_up = calc_jackpot_index(real_norm_g_stat, g_ci_lo)
        if not np.isnan(ji_ci_lo) and not np.isnan(ji_ci_up):
            logger.debug(f"{confidence_level * 100.} % confidence interval "
                         f"for jackpotting index: {ji_ci_lo} - {ji_ci_up}")
        # Stop when the confidence interval lies entirely below or above
        # max_jackpot_index, so it's unambiguous whether the jackpotting
        # index is less or greater than max_jackpot_index.
        # Do not use "if not ji_ci_lo <= max_jackpot_index <= ji_ci_up"
        # because this expression will evaluate to True after the first
        # iteration, when ji_ci_lo and ji_ci_up will both be NaN.
        if ji_ci_lo > max_jackpot_index or max_jackpot_index > ji_ci_up:
            logger.info("Ended boostrapping null normalized G-test statistics")
            return np.array(null_norm_g_stats)
