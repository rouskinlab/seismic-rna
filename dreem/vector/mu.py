from logging import getLogger

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov

from .bits import get_mut_info_bits

logger = getLogger(__name__)


def _calc_obs(p_real: np.ndarray, min_gap: int):
    """ Calculate the probability that a bit vector generated randomly
    from the given mutation rates would not have any mutations closer
    than ```min_gap```. Note that ```mus``` is transposed relative
    to all other uses, so that its shape is (positions x clusters).
    Transposition makes the indexing easier because this function uses
    the positions as the primary axis and clusters as secondary.

    Parameters
    ----------
    p_real: ndarray
        A (positions x clusters) array of the real mutation rates,
        i.e. corrected for the drop-out bias.
    min_gap: int
        Minimum number of non-mutated bases between two mutations;
        must be ≥ 0.

    Returns
    -------
    tuple[ndarray, ndarray]
        - 0: A 1D array of the probability for each cluster that, given
             the real mutation rates for each position in each cluster,
             a randomly generated bit vector coming from the cluster
             would have no two mutations closer than min_gap positions.
        - 1: A 2D (positions x clusters) array of the mutation rates
             that would be observed given the real mutation rates and
             the minimum gap between two mutations.
    """
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    # Determine the number of positions and clusters.
    dims = p_real.shape
    npos, ncls = dims
    if min_gap >= npos:
        raise ValueError(f"min_gap must be < number of positions ({npos}), "
                         f"but got {min_gap}")
    if min_gap == 0:
        # No mutations can be too close, so all observed probabilities
        # are 1.0 and the observed mutation rates equal the real rates.
        return np.ones(ncls), p_real.copy()
    # Compute the real non-mutation rates (q = 1 - p).
    q_real = 1. - p_real
    # Compute the cumulative sums of the log non-mutation rates.
    # Sum logarithms instead of multiply for better numerical stability.
    # The cumulative sums are split into three sections:
    end1begin2 = 1 + min_gap
    end2begin3 = end1begin2 + npos
    # - log_q_cumsum[0: 1 + min_gap]
    #   all 0.0
    log_q_cumsum = np.zeros((end2begin3 + min_gap, ncls), dtype=float)
    # - log_q_cumsum[1 + min_gap: 1 + min_gap + npos]
    #   cumulative sums of log non-mutation rates
    log_q_cumsum[end1begin2: end2begin3] = np.cumsum(np.log(q_real), axis=0)
    # - log_q_cumsum[1 + min_gap + npos: 1 + min_gap + npos + min_gap]
    #   all equal to final cumulative sum, log_q_cumsum[min_gap + npos]
    log_q_cumsum[end2begin3:] = log_q_cumsum[end2begin3 - 1]
    # For each window of (min_gap) positions, find the probability that
    # the window has no mutations, assuming mutations are independent.
    # The log probability is the sum over the window, which equals the
    # cumulative sum at the end of the window minus the cumulative sum
    # one index before the beginning of the window. Then apply np.exp().
    # The dimensions of q_window are (1 + npos + min_gap, ncls).
    q_window = np.exp(log_q_cumsum[min_gap:] - log_q_cumsum[: end2begin3])
    # For each position (j) in the sequence, calculate the probability
    # no_close[j] that no two mutations are too close, given that no two
    # mutations after position (j) are too close.
    # Equivalently, it is the probability that no two mutations from the
    # beginning of the sequence up to position (j) are too close:
    # If position (j) is mutated, then no two mutations up to (j) are
    # too close iff none of the previous (min_gap) positions are mutated
    # (P = pj_qwin[j]) and no two mutations before that window are too
    # close (P = no_close[j - (1 + min_gap)]).
    # If position (j) is not mutated (P = q_real[j]), then no two
    # mutations from the beginning up to (j) are too close iff no two
    # mutations up to (j - 1) are too close (P = no_close[j - 1]).
    # Combining these two situations gives this recurrence relation:
    # no_close[j] = (pj_qwin[j] * no_close[j - (1 + min_gap)]
    #                + q_real[j] * no_close[j - 1])
    # The initial condition is no_close[0] = 1.0 because there are
    # certainly no two mutations that are too close within the first
    # one position in the sequence.
    no_close_prev = np.ones(dims, dtype=float)
    no_close_next = np.ones(dims, dtype=float)
    # Keep track of the probabilities that none of the (min_gap) bases
    # preceding (following) base (j) are mutated and that there are no
    # mutations within (min_gap) positions before (after) those bases.
    no_close_win_prev = np.ones(dims, dtype=float)
    no_close_win_next = np.ones(dims, dtype=float)
    # This recurrence relation has no simple closed-form solution, and
    # likely no closed-form solution at all, so compute it via looping:
    for jp in range(1, npos):
        # Probability that none of the preceding (min_gap) bases are
        # mutated and no mutations before them are too close.
        no_close_win_prev[jp] = (q_window[jp] *
                                 no_close_prev[(jp - end1begin2) % npos])
        # Probability that no two mutations from the beginning to (jp)
        # are too close.
        no_close_prev[jp] = (p_real[jp] * no_close_win_prev[jp] +
                             q_real[jp] * no_close_prev[jp - 1])
    for jn in range(npos - 2, -1, -1):
        # Probability that none of the following (min_gap) bases are
        # mutated and no mutations after them are too close.
        no_close_win_next[jn] = (q_window[jn + end1begin2] *
                                 no_close_next[(jn + end1begin2) % npos])
        # Probability that no two mutations from (jn) to the end are too
        # close.
        no_close_next[jn] = (p_real[jn] * no_close_win_next[jn] +
                             q_real[jn] * no_close_next[jn + 1])
    # The probability that a randomly generated bit vector has no two
    # mutations that are too close is the probability that no two
    # mutations are too close after and including the first position.
    no_close = no_close_next[0]
    # It is also the probability that no two mutations are too close
    # before and including the last position.
    if not np.allclose(no_close, no_close_prev[-1]):
        raise ValueError("Probability of observation failed to converge.")
    # For each position (j), calculate the observed mutation rates given
    # the real mutation rates.
    # Start by calculating the joint probability that a bit vector is
    # observed (i.e. has no mutations that are too close together) and
    # position (j) is mutated: the product of the probabilities that
    # - position (j) is mutated: p_real[j]
    # - no bases within (min_gap) positions before (j) are mutated and
    #   no two mutations before them are too close: no_close_win_prev[j]
    # - no bases within (min_gap) positions after (j) are mutated and no
    #   two mutations after them are too close: no_close_win_next[j]
    # Then compute the conditional probability that position (j) is
    # mutated, given that the bit vector has no two mutations that are
    # too close, by dividing the joint probability by the probability
    # of the condition: no_close.
    p_obs = p_real * no_close_win_prev * no_close_win_next / no_close
    return no_close, p_obs


def denom(mus_real: np.ndarray, min_gap: int):
    """ Return a 1D array of the probability for each cluster that,
    given the real mutation rates for each position in each cluster, a
    randomly generated bit vector coming from the cluster would have no
    two mutations closer than min_gap positions. """
    return _calc_obs(mus_real, min_gap)[0]


def real_to_obs(mus_real: np.ndarray, min_gap: int):
    """ A 2D (positions x clusters) array of the mutation rates that
    would be observed given the real mutation rates and the minimum gap
    between two mutations. """
    return _calc_obs(mus_real, min_gap)[1]


def diff_real_obs(mus_real: np.ndarray, mus_obs: np.ndarray, min_gap: int):
    """ Compute the difference between the mutation rates that would be
    observed if ```mus_real``` were the real mutation rates (including
    unobserved reads), and the actual observed mutation rates.

    Parameters
    ----------
    mus_real: ndarray
        A (positions x clusters) array of the current guesses of each
        cluster's real mutation rates.
    mus_obs: ndarray
        A (positions x clusters) array of the actual observed mutation
        rates, from the weighted average over all bit vectors.
    min_gap: int
        Minimum permitted gap between two mutations.

    Returns
    -------
    ndarray
        A (positions x clusters) array of the difference between each
        expected-to-be-observed and each actual observed mutation rate.
    """
    return real_to_obs(mus_real, min_gap) - mus_obs


def obs_to_real(mus_obs: np.ndarray, min_gap: int,
                mus_guess: np.ndarray | None = None):
    """
    Given observed mutation rates ```mus_obs``` (which do not include
    any reads that dropped out because they had mutations closer than
    ```min_gap``` nt apart), estimate the real mutation rates that
    include these unobserved reads.

    Parameters
    ----------
    mus_obs: ndarray
        A (positions x clusters) array of the observed mutation rates,
        which do not include unobserved reads that dropped out.
    min_gap: int
        Minimum permitted gap between two mutations.
    mus_guess: ndarray
        Initial guess of the real mutation rates. If given, must be the
        same shape as mus_obs. If omitted, defaults to mus_obs, which is
        usually close to the optimal value.

    Returns
    -------
    ndarray
        A (positions x clusters) array of the real mutation rates that
        would be expected to yield the observed mutation rates.
    """
    # Determine the initial guess of the real mutation rates.
    if mus_guess is None:
        mus_guess = mus_obs
    elif mus_guess.shape != mus_obs.shape:
        raise ValueError(f"Dimensions of mus_guess {mus_guess.shape} and "
                         f"mus_obs {mus_obs.shape} differed")
    # Solve for the "real" mutation rates that yield minimal difference
    # between the mutation rates that would be expected to be observed
    # given the real mutation rates (i.e. when reads with mutations too
    # close are removed from the real mutation rates) and the actual
    # observed mutation rates. Use Newton's method, which finds the
    # parameters of a function that make it evaluate to zero, with the
    # Krylov approximation of the Jacobian, which improves performance.
    return newton_krylov(lambda mus_iter: diff_real_obs(mus_iter,
                                                        mus_obs,
                                                        min_gap),
                         mus_guess)


def reads_to_mus_obs(members: pd.DataFrame,
                     muts: pd.DataFrame,
                     refs: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    Calculate the mutation rate at each position in a mutational profile
    for one or more clusters.

    Parameters
    ----------
    members: DataFrame
        Cluster membership: each index (i) is the name of a read, each
        column (k) the name of a cluster, and each value (i, k) the
        likelihood that read (i) came from cluster (k). Any read absent
        from members is ignored when calculating mutations.
    muts: DataDrame
        Mutated bits: each index (i) is the name of a read, each column
        (j) a position in the profile, and each value (i, j) a boolean
        where True means that position (j) of read (i) is mutated and
        False means that it is not.
    refs: DataDrame | None = None
        Reference bits: each index (i) is the name of a read, each
        column (j) a position in the profile, and each value (i, j) a
        boolean where True means that position (j) of read (i) matches
        the reference sequence and False means that it does not.
        If omitted, defaults to the logical not of ```mutbits```.

    Returns
    -------
    DataFrame
        Mutation rates: each index (j) is a position in the profile,
        each column (k) the name of a cluster, and each value (j, k)
        the mutation fraction at position (j) in cluster (k).
    """
    # Ensure dimensions and axes match.
    muts, info = get_mut_info_bits(muts, refs)
    # The observed cluster mutation rate of each position is the number
    # of mutated bits divided by the number of informative bits.
    return (muts.loc[members.index].T.dot(members) /
            info.loc[members.index].T.dot(members))


def mus_obs_to_real(mus_obs: pd.DataFrame, min_gap: int):
    """
    Parameters
    ----------
    mus_obs: DataFrame
        A (positions x clusters) array of the observed mutation rate at
        each position in each cluster.
    min_gap: int
        Minimum distance between mutations. This parameter is used to
        correct the bias in observed mutation rates that is caused by
        dropout of reads with nearby mutations, but NOT to remove any
        rows of mutbits or refbits with nearby mutations. Instead, the
        'members' parameter is responsible for filtering such vectors,
        because it comes from the clustering step in which vectors with
        nearby mutations were already removed.

    Returns
    -------
    DataFrame
        A (positions x clusters) array of the real mutation rate at each
        position in each cluster.
    """
    # Estimate the real mutation rates.
    return pd.DataFrame(obs_to_real(mus_obs.values, min_gap),
                        index=mus_obs.index, columns=mus_obs.columns)


def reads_to_mus_real(members: pd.DataFrame, min_gap: int,
                      muts: pd.DataFrame, refs: pd.DataFrame | None = None):
    """ A convenience function to chain ```reads_to_mus_obs``` and
    ```mus_obs_to_real```. """
    return mus_obs_to_real(reads_to_mus_obs(members, muts, refs), min_gap)