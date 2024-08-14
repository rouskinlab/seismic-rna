from typing import Iterable

import numpy as np
from numba import jit, int64

from .marginal import calc_marginal
from ..core.array import get_length, ensure_same_length, find_dims
from ..core.mu import calc_p_noclose, calc_p_noclose_given_ends

SUM_EXP_PRECISION = 3
POSITIONS = "positions"
CLUSTERS = "clusters"


def linearize_ends_matrix(p_ends: np.ndarray):
    """ Linearize an N x N matrix of end coordinate probabilities/counts
    into a length-N list of 5'/3' coordinates and their values. """
    end5s, end3s = np.nonzero(np.triu(p_ends))
    return end5s, end3s, p_ends[end5s, end3s]


@jit()
def _sim_reads(p_mut: np.ndarray,
               end5s: np.ndarray,
               end3s: np.ndarray,
               end_choices: np.ndarray):
    """ Simulate reads.

    Parameters
    ----------
    end_choices: np.ndarray
        Number of reads to simulate
    p_mut: np.ndarray
        Probability that each position is mutated (1D, float)
    end5s: np.ndarray
        5' end coordinates: 0-indexed (1D, int)
    end3s: np.ndarray
        3' end coordinates: 0-indexed (1D, int)
    """
    n_reads, = end_choices.shape
    # Initialize reads to zeros; reserve the last two columns for the
    # 5' and 3' end coordinates.
    reads = np.zeros((n_reads, p_mut.size + 2), dtype=int64)
    for i, end in enumerate(end_choices):
        # Choose a pair of end coordinates.
        end5 = end5s[end]
        end3 = end3s[end]
        reads[i, -2] = end5
        reads[i, -1] = end3
        # Choose which positions are mutated.
        n_pos = end3 - end5 + 1
        muts = end5 + np.flatnonzero(np.random.random(n_pos)
                                     < p_mut[end5: end3 + 1])
        reads[i, muts] = 1
    return reads


def _filter_reads(reads: np.ndarray, min_mut_gap: int):
    """ Remove reads with two mutations closer than `min_mut_gap`. """
    n_reads, n_pos = reads.shape
    # It is impossible for two mutations to be separated by more than
    # (n_pos - 1) positions.
    min_mut_gap = min(min_mut_gap, n_pos - 1)
    # Initially, select all reads.
    selected = np.ones(n_reads, dtype=bool)
    if min_mut_gap <= 0:
        # No two mutations can be too close.
        return selected
    # Using a sliding window of (min_mut_gap + 1) positions, flag reads
    # with ≥ 2 mutations within any window.
    offset = min_mut_gap + 1
    for i in range(n_pos - min_mut_gap):
        j = i + offset
        selected &= np.count_nonzero(reads[:, i: j], axis=1) < 2
    return selected


def sim_obs_exp(p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                n_reads: int,
                min_mut_gap: int,
                unmasked: np.ndarray):
    """ Simulate observed and expected counts. """
    rng = np.random.default_rng()
    dims = find_dims([(POSITIONS, CLUSTERS),
                      (POSITIONS, POSITIONS),
                      (CLUSTERS,)],
                     [p_mut, p_ends, p_clust],
                     ["p_mut", "p_ends", "p_clust"])
    # Linearize the end coordinates and their probabilities.
    end5s, end3s, p_ends_linear = linearize_ends_matrix(p_ends)
    n_clust = dims[CLUSTERS]
    while True:
        # Choose how many reads to put in each cluster.
        if n_clust > 1:
            n_reads_per_clust = np.unique(rng.choice(np.arange(n_clust),
                                                     n_reads,
                                                     p=p_clust),
                                          return_counts=True)[1]
        elif n_clust == 1:
            n_reads_per_clust = np.array([n_reads])
        else:
            raise ValueError(f"Must have ≥ 1 cluster(s), but got {n_clust}")
        if len(n_reads_per_clust) != n_clust:
            # This should be impossible, but checking just in case.
            raise ValueError("len(n_reads_per_clust) != n_clust")
        # Simulate reads for each cluster.
        reads = [_sim_reads(p_mut[:, k],
                            end5s,
                            end3s,
                            rng.choice(p_ends_linear.size, n, p=p_ends_linear))
                 for k, n in enumerate(n_reads_per_clust)]
        if len(reads) == 1:
            reads = reads[0]
        else:
            reads = np.vstack(reads)
        if min_mut_gap > 0:
            # Remove reads with two mutations too close.
            reads = reads[_filter_reads(reads[:, :-2], min_mut_gap)]
        # Count the number of times each unique read was observed.
        uniq_reads, num_obs = np.unique(reads, axis=0, return_counts=True)
        # List the reads with mutations at each unmasked position.
        muts_per_pos = [np.flatnonzero(uniq_reads[:, j]) for j in unmasked]
        log_exp = np.log(n_reads) + calc_marginal(p_mut,
                                                  p_ends,
                                                  p_clust,
                                                  uniq_reads[:, -2],
                                                  uniq_reads[:, -1],
                                                  unmasked,
                                                  muts_per_pos,
                                                  min_mut_gap,
                                                  calc_resps=False)
        yield num_obs, log_exp


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
        g_stat = 2. * np.sum(num_obs[at_least_min_exp]
                             * (np.log(num_obs[at_least_min_exp])
                                - log_exp[at_least_min_exp]))
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
            g_stat += 2. * (num_obs_less_than_min_exp
                            * np.log(num_obs_less_than_min_exp
                                     / num_exp_less_than_min_exp))
            # Degree of freedom for all reads with expected counts less
            # than min_exp.
            df_less_than_min_exp = 1
        else:
            df_less_than_min_exp = 0
    else:
        # Calculate the G-test statistic of all reads.
        g_stat = 2. * np.sum(num_obs * (np.log(num_obs) - log_exp))
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


def normalize_g_stat(g_stat: float, n_reads: int):
    """ Normalize the G-test statistic. """
    if n_reads == 0:
        if g_stat != 0.:
            raise ValueError(
                f"If n_reads is 0, then g_stat must be 0, but got {g_stat}"
            )
        return 0.
    return g_stat / n_reads


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
                       null_norm_g_stat: float | list[float] | np.ndarray):
    """ Calculate the jackpotting index.

    Parameters
    ----------
    real_norm_g_stat: float
        Normalized G-test statistic of the real dataset.
    null_norm_g_stat: float | list[float] | np.ndarray
        Normalized G-test statistic(s) based on the null model.

    Returns
    -------
    float
        Jackpotting index (log-ratio of the real to the null).
    """
    # If null_g_stat is array-like, then find its mean.
    null_norm_g_stat = np.atleast_1d(null_norm_g_stat)
    if get_length(null_norm_g_stat, "null_g_stat") > 0:
        null_norm_g_stat = null_norm_g_stat.mean()
    else:
        null_norm_g_stat = np.nan
    return real_norm_g_stat - null_norm_g_stat


def bootstrap_norm_g_stats(p_mut: np.ndarray,
                           p_ends: np.ndarray,
                           p_clust: np.ndarray,
                           n_reads: int,
                           min_mut_gap: int,
                           unmasked: np.ndarray,
                           real_norm_g_stat: float,
                           confidence_level: float,
                           max_jackpot_index: float):
    """ Bootstrap the mean observed/expected log-ratios.

    Parameters
    ----------
    p_mut: np.ndarray
        (positions x clusters) Probability that each position is mutated
    p_ends: np.ndarray
        (positions x positions) Probability that a read has each pair of
        5' and 3' end coordinates
    p_clust: np.ndarray
        (clusters) Probability that a read comes from each cluster
    n_reads: int
        Number of reads in the real dataset.
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
    """
    # Inflate the number of reads because some will drop out due to the
    # observer bias.
    p_noclose_given_ends = calc_p_noclose_given_ends(p_mut, min_mut_gap)
    p_noclose = np.vdot(calc_p_noclose(p_ends, p_noclose_given_ends), p_clust)
    if p_noclose == 0.:
        raise ValueError("Observing a read is impossible given the mutation "
                         f"rates {p_mut} and minimum gap of {min_mut_gap}")
    n_reads_sim = round(n_reads / p_noclose)
    # Simulate observed and expected read counts.
    null_norm_g_stats = list()
    for num_obs, log_exp in sim_obs_exp(p_mut,
                                        p_ends,
                                        p_clust,
                                        n_reads_sim,
                                        min_mut_gap,
                                        unmasked):
        # Calculate and normalize this null model G-test statistic.
        null_g_stat, _ = calc_jackpot_g_stat(num_obs, log_exp)
        null_norm_g_stats.append(normalize_g_stat(null_g_stat, num_obs.sum()))
        # Determine a confidence interval for the mean of all normalized
        # G-test statistics simulated so far.
        g_ci_lo, g_ci_up = calc_norm_g_stat_ci(null_norm_g_stats,
                                               confidence_level)
        # Determine a confidence interval for the jackpotting index.
        # Since the jackpotting index is inversely related to the null
        # G-test statistic, the lower and upper bounds are swapped.
        i_ci_lo = calc_jackpot_index(real_norm_g_stat, g_ci_up)
        i_ci_up = calc_jackpot_index(real_norm_g_stat, g_ci_lo)
        if max_jackpot_index < i_ci_lo or i_ci_up < max_jackpot_index:
            return np.array(null_norm_g_stats)
