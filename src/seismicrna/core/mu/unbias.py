from logging import getLogger

import numpy as np

logger = getLogger(__name__)

# Maximum allowed mutation rate.
EPSILON = 1.e-6
MAX_MU = 1. - EPSILON


def clip(mus: np.ndarray):
    """ Check if any mutation rates are < 0, ≥ 1, or NaN. If so, then
    fill any NaN values with 0 and clip all values to [0, 1). """
    if mus.size == 0:
        # There is nothing to clip.
        return mus
    mu_min = np.min(mus)
    mu_max = np.max(mus)
    if mu_min < 0. or mu_max > MAX_MU or np.isnan(mu_min) or np.isnan(mu_max):
        clipped = np.clip(np.nan_to_num(mus), 0., MAX_MU)
        logger.warning(f"Mutation rates outside [0, {MAX_MU}]:\n"
                       f"{mus[np.nonzero(mus != clipped)]}")
        return clipped
    return mus


def _calc_obs(mu_adj: np.ndarray, min_gap: int):
    """ Calculate the probability that a bit vector generated randomly
    from the given mutation rates would not have any mutations closer
    than `min_gap`. Note that `mus` is transposed relative
    to all other uses, so that its shape is (positions x clusters).
    Transposition makes the indexing easier because this function uses
    the positions as the primary axis and clusters as secondary.

    Parameters
    ----------
    mu_adj: ndarray
        A (positions x clusters) array of the adjusted mutation rates,
        i.e. corrected for the observer bias.
    min_gap: int
        Minimum number of non-mutated bases between two mutations;
        must be ≥ 0.

    Returns
    -------
    tuple[ndarray, ndarray]
        - A 1D array of the probability for each cluster that, given the
          real mutation rates for each position in each cluster, a
          randomly generated bit vector coming from the cluster would
          have no two mutations closer than `min_gap` positions.
        - A 2D (positions x clusters) array of the mutation rates that
          would be observed given the real mutation rates and the
          minimum gap between two mutations.
    """
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    # Determine the number of positions and clusters.
    dims = mu_adj.shape
    npos, ncls = dims
    # The number of non-mutated bases next to a mutated base cannot
    # exceed the total number of positions minus one.
    if min_gap >= npos:
        min_gap = npos - 1
    if min_gap == 0 or npos == 0 or ncls == 0:
        # No mutations can be too close, so all observed probabilities
        # are 1.0 and the observed mutation rates equal the real rates.
        return np.ones(ncls), mu_adj.copy()
    # Compute the adjusted non-mutation rates (nu = 1 - mu).
    nu_adj = 1. - mu_adj
    # Verify that all positions mutate with probabilty < 1.
    if np.count_nonzero(nu_adj) < nu_adj.size:
        raise ValueError(f"Got mutation rates of 1.0 in {mu_adj}")
    # Compute the cumulative sums of the log non-mutation rates.
    # Sum logarithms instead of multiply for better numerical stability.
    # The cumulative sums are split into three sections:
    end1begin2 = 1 + min_gap
    end2begin3 = end1begin2 + npos
    # - log_nu_cumsum[0: 1 + min_gap]
    #   all 0.0
    log_nu_cumsum = np.zeros((end2begin3 + min_gap, ncls), dtype=float)
    # - log_nu_cumsum[1 + min_gap: 1 + min_gap + npos]
    #   cumulative sums of log non-mutation rates
    log_nu_cumsum[end1begin2: end2begin3] = np.cumsum(np.log(nu_adj), axis=0)
    # - log_nu_cumsum[1 + min_gap + npos: 1 + min_gap + npos + min_gap]
    #   all equal to final cumulative sum, log_nu_cumsum[min_gap + npos]
    log_nu_cumsum[end2begin3:] = log_nu_cumsum[end2begin3 - 1]
    # For each window of (min_gap) positions, find the probability that
    # the window has no mutations, assuming mutations are independent.
    # The log probability is the sum over the window, which equals the
    # cumulative sum at the end of the window minus the cumulative sum
    # one index before the beginning of the window. Then apply np.exp().
    # The dimensions of nu_win are (1 + npos + min_gap, ncls).
    nu_win = np.exp(log_nu_cumsum[min_gap:] - log_nu_cumsum[: end2begin3])
    # For each position (j) in the sequence, calculate the probability
    # f_obs[j] that no two mutations are too close, given that no two
    # mutations after position (j) are too close.
    # Equivalently, it is the probability that no two mutations from the
    # beginning of the sequence up to position (j) are too close:
    # If position (j) is mutated, then no two mutations up to (j) are
    # too close iff none of the previous (min_gap) positions are mutated
    # (P = pj_qwin[j]) and no two mutations before that window are too
    # close (P = f_obs[j - (1 + min_gap)]).
    # If position (j) is not mutated (P = nu_adj[j]), then no two
    # mutations from the beginning up to (j) are too close iff no two
    # mutations up to (j - 1) are too close (P = f_obs[j - 1]).
    # Combining these two situations gives this recurrence relation:
    # f_obs[j] = (pj_qwin[j] * f_obs[j - (1 + min_gap)]
    #                + nu_adj[j] * f_obs[j - 1])
    # The initial condition is f_obs[0] = 1.0 because there are
    # certainly no two mutations that are too close within the first
    # one position in the sequence.
    f_obs_prev = np.ones(dims, dtype=float)
    f_obs_next = np.ones(dims, dtype=float)
    # Keep track of the probabilities that none of the (min_gap) bases
    # preceding (following) base (j) are mutated and that there are no
    # mutations within (min_gap) positions before (after) those bases.
    f_obs_win_prev = np.ones(dims, dtype=float)
    f_obs_win_next = np.ones(dims, dtype=float)
    # This recurrence relation has no simple closed-form solution, and
    # likely no closed-form solution at all, so compute it iteratively:
    for jp in range(1, npos):
        # Probability that none of the preceding (min_gap) bases are
        # mutated and no mutations before them are too close.
        f_obs_win_prev[jp] = (nu_win[jp] *
                              f_obs_prev[(jp - end1begin2) % npos])
        # Probability that no two mutations from the beginning to (jp)
        # are too close.
        f_obs_prev[jp] = (mu_adj[jp] * f_obs_win_prev[jp] +
                          nu_adj[jp] * f_obs_prev[jp - 1])
    for jn in range(npos - 2, -1, -1):
        # Probability that none of the following (min_gap) bases are
        # mutated and no mutations after them are too close.
        f_obs_win_next[jn] = (nu_win[jn + end1begin2] *
                              f_obs_next[(jn + end1begin2) % npos])
        # Probability that no two mutations from (jn) to the end are too
        # close.
        f_obs_next[jn] = (mu_adj[jn] * f_obs_win_next[jn] +
                          nu_adj[jn] * f_obs_next[jn + 1])
    # The probability that a randomly generated bit vector has no two
    # mutations that are too close is the probability that no two
    # mutations are too close after and including the first position.
    f_obs = f_obs_next[0]
    # For each position (j), calculate the observed mutation rates given
    # the real mutation rates.
    # Start by calculating the joint probability that a bit vector is
    # observed (i.e. has no mutations that are too close together) and
    # position (j) is mutated: the product of the probabilities that
    # - position (j) is mutated: p_adj[j]
    # - no bases within (min_gap) positions before (j) are mutated and
    #   no two mutations before them are too close: f_obs_win_prev[j]
    # - no bases within (min_gap) positions after (j) are mutated and no
    #   two mutations after them are too close: f_obs_win_next[j]
    # Then compute the conditional probability that position (j) is
    # mutated, given that the bit vector has no two mutations that are
    # too close, by dividing the joint probability by the probability
    # that no two mutations are too close: f_obs.
    mu_obs = mu_adj * f_obs_win_prev * f_obs_win_next / f_obs
    return f_obs, mu_obs


def _calc_f_obs_clipped(mus_adj: np.ndarray, min_gap: int):
    """ Return a 1D array of the probability for each cluster that,
    given the real mutation rates for each position in each cluster, a
    randomly generated bit vector coming from the cluster would have no
    two mutations closer than min_gap positions. """
    return _calc_obs(mus_adj, min_gap)[0]


def calc_f_obs_numpy(mus_adj: np.ndarray, min_gap: int):
    """ Return a 1D array of the probability for each cluster that,
    given the real mutation rates for each position in each cluster, a
    randomly generated bit vector coming from the cluster would have no
    two mutations closer than min_gap positions. """
    if mus_adj.ndim == 2:
        return _calc_f_obs_clipped(clip(mus_adj), min_gap)
    if mus_adj.ndim == 1:
        return float(np.squeeze(calc_f_obs_numpy(mus_adj.reshape((-1, 1)),
                                                 min_gap),
                                axis=0))
    raise ValueError(f"Expected 1 or 2 dimensions, but got {mus_adj.ndim}")


def _calc_mu_obs(mus_adj: np.ndarray, min_gap: int):
    """ A 2D (positions x clusters) array of the mutation rates that
    would be observed given the real mutation rates and the minimum gap
    between two mutations. """
    return _calc_obs(clip(mus_adj), min_gap)[1]


def _diff_adj_obs(mus_adj: np.ndarray, mus_obs: np.ndarray, min_gap: int):
    """ Compute the difference between the mutation rates that would be
    observed if `mus_adj` were the real mutation rates (including
    unobserved reads), and the actual observed mutation rates.

    Parameters
    ----------
    mus_adj: ndarray
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
    return _calc_mu_obs(mus_adj, min_gap) - mus_obs


def calc_mu_adj_numpy(mus_obs: np.ndarray, min_gap: int,
                      mus_guess: np.ndarray | None = None,
                      f_tol: float = 5e-1,
                      f_rtol: float = 5e-1,
                      x_tol: float = 1e-4,
                      x_rtol: float = 5e-1):
    """
    Given observed mutation rates `mus_obs` (which do not include
    any reads that dropped out because they had mutations closer than
    `min_gap` nt apart), estimate the real mutation rates that
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
    f_tol: float = 1e-4
        Absolute tolerance in residual.
    f_rtol: float = 1e-0
        Relative tolerance in residual.
    x_tol: float = 1e-4
        Absolute tolerance in step.
    x_rtol: float = 1e-0
        Relative tolerance in step.

    Returns
    -------
    ndarray
        A (positions x clusters) array of the real mutation rates that
        would be expected to yield the observed mutation rates.
    """
    # Import scipy here instead of at the top of this module because
    # its import is slow enough to impact global startup time.
    from scipy.optimize import newton_krylov
    if mus_obs.ndim == 1:
        # Expand the mutation rates from a vector into a matrix with one
        # column (cluster), calculate, and squeeze back into a vector.
        return np.squeeze(calc_mu_adj_numpy(mus_obs.reshape(-1, 1), min_gap),
                          axis=1)
    if mus_obs.ndim != 2:
        # Only matrices (positions x clusters) are supported hereafter.
        raise ValueError(f"Expected 1 or 2 dimensions, but got {mus_obs.ndim}")
    # Clip any invalid mutation rates.
    mus_obs = clip(mus_obs)
    if mus_obs.size == 0:
        # If the array is empty, then there are no real mutation rates.
        return mus_obs
    # Determine the initial guess of the real mutation rates.
    if mus_guess is None:
        mus_guess = mus_obs
    else:
        # The dimensions of the guess and the mutation rates must match.
        if mus_guess.shape != mus_obs.shape:
            raise ValueError(f"Dimensions of mus_guess {mus_guess.shape} and "
                             f"mus_obs {mus_obs.shape} differed")
        # Clip any invalid guessed mutation rates.
        mus_guess = clip(mus_guess)
    # Solve for the "real" mutation rates that yield minimal difference
    # between the mutation rates that would be expected to be observed
    # given the real mutation rates (i.e. when reads with mutations too
    # close are removed from the real mutation rates) and the actual
    # observed mutation rates. Use Newton's method, which finds the
    # parameters of a function that make it evaluate to zero, with the
    # Krylov approximation of the Jacobian, which improves performance.
    return clip(newton_krylov(lambda mus_iter: _diff_adj_obs(mus_iter,
                                                             mus_obs,
                                                             min_gap),
                              mus_guess,
                              f_tol=f_tol, f_rtol=f_rtol,
                              x_tol=x_tol, x_rtol=x_rtol))


def calc_prop_adj_numpy(prop_obs: np.ndarray, f_obs: np.ndarray):
    """ Calculate the adjusted proportion of the clusters given their
    observed proportions and the observer bias. """
    weighted_prop_obs = prop_obs / f_obs
    return weighted_prop_obs / np.sum(weighted_prop_obs)

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
