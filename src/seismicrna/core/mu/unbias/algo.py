"""

Mutation Rate Core Module

========================================================================

The functions in this module serve two main purposes:
 1. Adjust mutation rates to correct for observer bias.
 2. Normalize and winsorize mutation rates

------------------------------------------------------------------------

Adjust mutation rates to correct for observer bias

Our lab has found that pairs of mutations in DMS-MaPseq data rarely have
fewer than three non-mutated bases separating them. We suspect that the
reverse transcriptase is prone to falling off or stalling at locations
in the RNA where DMS has methylated bases that are too close.

Regardless of the reason, mutations on nearby bases are not independent,
which violates a core assumption of the Bernoulli mixture model that we
use in the expectation-maximization clustering algorithm. Namely, that
mutations occur independently of each other, such that the likelihood of
observing a bit vector is the product of the mutation rate of each base
that is mutated and one minus the mutation rate ("non-mutation rate") of
each base that is not mutated.

In order to use the Bernoulli mixture model for expectation-maximization
clustering, we modify it such that bases separated by zero, one, or two
other bases are no longer assumed to mutate independently. Specifically,
since pairs of close mutations are rare, we add the assumption that no
mutations separated by fewer than three non-mutated bases may occur, and
exclude the few bit vectors that have such close mutations.

When these bit vectors are assumed to exist in the original population
of RNA molecules but not appear in the observed data, the mutation rates
that are observed will differ from the real, underlying mutation rates,
which would include the unobserved bit vectors. The real mutation rates
must therefore be estimated from the observed mutation rates.

It is relatively easy to estimate the observed mutation rates given the
real mutation rates, but there is no closed-form solution that we know
of to estimate the real mutation rates from the observed mutation rates.
Thus, we use an iterative approach, based on Newton's method for finding
the roots of functions. We initially guess the real mutation rates. Each
iteration, we estimate the mutation rates that would have been observed
given our current guess of the real mutation rates, and then subtract
the mutation rates that were actually observed. This difference would be
zero if we had accurately guessed the real mutation rates. Thus, we use
Newton's method to solve for the real mutation rates that minimize this
difference. The details are described in the comments of this module and
in Tomezsko et al. (2020) (https://doi.org/10.1038/s41586-020-2253-5).

------------------------------------------------------------------------

Normalize and winsorize mutation rates

The mutation rates of an RNA may vary among different samples because of
variations in the chemical probing and mutational profiling procedure.
Thus, it is often helpful to compute "normalized" mutation rates that
can be compared directly across different samples and used to predict
secondary structures.

This module currently provides a simple method for normalizing mutation
rates. First, a specific quantile of the dataset is chosen, such as 0.95
(i.e. the 95th percentile). The mutation rate with this quantile is set
to 1.0, and all other mutation rates are scaled linearly.

If the chosen quantile is less than 1.0, then any mutation rates above
the quantile will be scaled to values greater than 1.0. Since these high
mutation rates may be exceptionally reactive bases, it is often helpful
to cap the normalized mutation rates to a maximum of 1.0. The winsorize
function in this module performs normalization and then sets any value
greater than 1.0 to 1.0.

"""

from logging import getLogger

import numpy as np

logger = getLogger(__name__)


def _clip(p_mut: np.ndarray):
    """ Check if any mutation rates are < 0, ≥ 1, or NaN. If so, then
    fill any NaN values with 0 and clip all values to [0, 1). """
    return np.clip(np.nan_to_num(p_mut), 0., 1.)


def _normalize(p_ends: np.ndarray):
    if p_ends.ndim == 2:
        # Handle the 2-dimensional case.
        return _normalize(np.expand_dims(p_ends, 2))[:, :, 0]
    p_ends = _clip(p_ends)
    if p_ends.size == 0:
        # There is nothing to normalize.
        return p_ends
    # Divide by the sum over both positional axes.
    npos, _, ncls = p_ends.shape
    p_sums = np.zeros(ncls)
    for j in range(npos):
        p_sums += np.sum(p_ends[j, j:], axis=0)
    return p_ends / p_sums


def _adjust_min_gap(num_pos: int, min_gap: int):
    """ Given the number of positions (`npos`) and the desired minimum
    gap between mutations (`min_gap`), find the minimum gap between
    mutations that is no smaller than 0 and, if possible, no larger
    than 1 less than the number of positions.

    Parameters
    ----------
    num_pos: int
        Number of positions in the reference section.
    min_gap: int
        Minimum number of non-mutated bases between two mutations.

    Returns
    -------
    int
        Adjusted minimum gap between mutations.
    """
    return max(min(min_gap, num_pos - 1), 0)


def _calc_p_nomut_window(p_mut_given_span: np.ndarray, min_gap: int):
    """ Given underlying mutation rates (`p_mut_given_span`), find the
    probability of no mutations in each window of size 0 to `min_gap`.

    Parameters
    ----------
    p_mut_given_span: ndarray
        2D (positions x clusters) array of the underlying mutation rates
        (i.e. the probability that a read has a mutation at position (j)
        given that it contains that position).
    min_gap: int
        Minimum number of non-mutated bases between two mutations.

    Returns
    -------
    np.ndarray
        3D (window x positions x clusters) array of the probability that
        (window) consecutive bases, ending at position (position), would
        have 0 mutations at all.
    """
    # Find and validate the dimensions.
    npos, ncls = p_mut_given_span.shape
    min_gap = _adjust_min_gap(npos, min_gap)
    # Determine the shape of the array to return.
    p_nomut_window_shape = min_gap + 1, npos + 1, ncls
    if min_gap == 0:
        # If min_gap == 0, then only size-0 windows must be considered,
        # which always have 0 mutations.
        return np.ones(p_nomut_window_shape)
    # Compute product of the non-mutation rates over each window ranging
    # in size from 0 to min_gap; element [g, j, k] is the product of all
    # non-mutation rates from position (j - g) to (j - 1) for cluster k.
    p_nomut_window = np.empty(p_nomut_window_shape)
    # For window size 0, every window represents 0 positions, for which
    # the probability of having 0 mutations is always 1.
    p_nomut_window[0] = 1.
    # For window size 1, every window represents exactly 1 position, so
    # the product over the window is just the non-mutation rate of the
    # corresponding position.
    p_nomut_window[1, 1:] = 1. - p_mut_given_span
    # For each window size from 2 to min_gap, calculate the probability
    # of there being no mutations in each window of that size.
    for size in range(2, p_nomut_window.shape[0]):
        shift = size - 1
        p_nomut_window[size, size:] = (p_nomut_window[1, size:]
                                       * p_nomut_window[shift, shift: -1])
    return p_nomut_window


def _calc_p_noclose_given_ends(p_mut_given_span: np.ndarray,
                               p_nomut_window: np.ndarray):
    """ Given underlying mutation rates (`p_mut_given_span`), calculate
    the probability that a read starting at position (a) and ending at
    position (b) would have no two mutations too close (i.e. separated
    by fewer than `min_gap` non-mutated positions), for each combination
    of (a) and (b) such that 1 ≤ a ≤ b ≤ L (in biological coordinates)
    or 0 ≤ a ≤ b < L (in Python coordinates).

    Parameters
    ----------
    p_mut_given_span: np.ndarray
        2D (positions x clusters) array of the underlying mutation rates
        (i.e. the probability that a read has a mutation at position (j)
        given that it contains that position).
    p_nomut_window: np.ndarray
        3D (window x positions x clusters) array of the probability that
        (window) consecutive bases, ending at position (position), would
        have zero mutations at all.

    Returns
    -------
    np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a random read starting at position (a) (row) and ending at
        position (b) (column) would have no two mutations too close.
    """
    # Find and validate the dimensions.
    npos, ncls = p_mut_given_span.shape
    inc_gap = p_nomut_window.shape[0]
    min_gap = inc_gap - 1
    if min_gap == 0:
        # If min_gap == 0, then no mutations can be too close, so the
        # probability of not having two mutations too close is 1.
        return np.ones((npos, npos, ncls))
    # For each pair of positions (i, j), find the probability that a
    # random bit vector from (i) to (j), inclusive, would have no two
    # mutations closer than min_gap positions: p_noclose_given_ends[i, j].
    p_noclose_given_ends = np.empty((npos, npos, ncls))
    # First, fill the main diagonal and (min_gap) sub-diagonals with 1
    # because length-1 reads cannot have two mutations too close.
    diag_indexes = np.arange(npos)
    p_noclose_given_ends[diag_indexes, diag_indexes] = 1.
    for size in range(1, inc_gap):
        p_noclose_given_ends[diag_indexes[size:], diag_indexes[:-size]] = 1.
    # The probabilities follow a recurrence relation that is assumed to
    # have no closed-form solution but can be computed via a loop that
    # fills p_noclose_given_ends one column (j) at a time.
    window_indexes = np.minimum(np.arange(npos, 0, -1), min_gap)
    for j in range(1, npos):
        # If position (j) is mutated (probability = p_mut_given_pos[j]),
        # then no two mutations from (i) to (j) are too close iff none
        # of the previous (min_gap) positions, i.e. (j - min_gap) to
        # (j - 1), are mutated (probability = p_nomut_window[j]); and if
        # no two mutations from (i) to (j - (min_gap + 1)) are too close
        # (probability = p_noclose_given_ends[i, j - (min_gap + 1)]).
        p_noclose_given_ends_mutj = (
                p_nomut_window[window_indexes[(npos - j):], j]
                * p_noclose_given_ends[:j, max(j - inc_gap, 0)]
        )
        # Otherwise, no two mutations from (i) to (j) are too close iff
        # the same is true from positions (i) to (j - 1) (probability
        # = p_noclose_given_ends[i, j - 1]).
        p_noclose_given_ends_nomutj = p_noclose_given_ends[:j, j - 1]
        # The probability that no two mutations from (i) to (j) are too
        # close is the weighted sum of these mutually exclusive events.
        p_noclose_given_ends[:j, j] = (
                p_noclose_given_ends_mutj * p_mut_given_span[j]
                + p_noclose_given_ends_nomutj * p_nomut_window[1, j + 1]
        )
    return p_noclose_given_ends


def calc_p_noclose_given_ends_numpy(p_mut_given_span: np.ndarray, min_gap: int):
    """ Given underlying mutation rates (`p_mut_given_span`), calculate
    the probability that a read starting at position (a) and ending at
    position (b) would have no two mutations too close (i.e. separated
    by fewer than `min_gap` non-mutated positions), for each combination
    of (a) and (b) such that 1 ≤ a ≤ b ≤ L (in biological coordinates)
    or 0 ≤ a ≤ b < L (in Python coordinates).

    Parameters
    ----------
    p_mut_given_span: ndarray
        A 2D (positions x clusters) array of the underlying mutation
        rates, i.e. the probability that a read has a mutation at
        position (j) given that it contains position (j).
    min_gap: int
        Minimum number of non-mutated bases between two mutations;
        must be ≥ 0.

    Returns
    -------
    np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a random read starting at position (a) (row) and ending at
        position (b) (column) would have no two mutations too close.
    """
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    if p_mut_given_span.ndim == 2:
        p_mut_given_span = _clip(p_mut_given_span)
        return _calc_p_noclose_given_ends(p_mut_given_span,
                                          _calc_p_nomut_window(p_mut_given_span,
                                                               min_gap))
    if p_mut_given_span.ndim == 1:
        return calc_p_noclose_given_ends_numpy(
            p_mut_given_span.reshape((-1, 1)), min_gap
        ).reshape((p_mut_given_span.size, p_mut_given_span.size))
    raise ValueError("p_mut_given_span must have 1 or 2 dimensions, "
                     f"but got {p_mut_given_span.ndim}")


def _calc_p_mut_given_span_noclose(p_mut_given_span: np.ndarray,
                                   p_ends: np.ndarray,
                                   p_noclose_given_ends: np.ndarray,
                                   p_nomut_window: np.ndarray):
    # Find and validate the dimensions.
    npos, ncls = p_mut_given_span.shape
    if p_ends.shape != (npos, npos):
        raise ValueError(f"Expected p_ends to be shaped {npos, npos}, "
                         f"but got {p_ends.shape}")
    if p_noclose_given_ends.shape != (npos, npos, ncls):
        raise ValueError(
            f"Expected p_noclose_given_ends to be shaped {npos, npos, ncls}, "
            f"but got {p_noclose_given_ends.shape}"
        )
    inc_gap = p_nomut_window.shape[0]
    if p_nomut_window.shape != (inc_gap, npos + 1, ncls):
        raise ValueError(
            f"Expected p_nomut_window to be shaped {inc_gap, npos + 1, ncls}, "
            f"but got {p_nomut_window.shape}"
        )
    min_gap = inc_gap - 1
    if min_gap == 0:
        # If min_gap == 0, then no mutations can be too close, so return
        # the original mutation rates.
        return p_mut_given_span
    # Cache the sum of each row and column of the probabilities that no
    # two mutations are too close, weighted by the end distributions.
    p_noclose_and_ends_rows = np.empty_like(p_mut_given_span)
    p_noclose_and_ends_rows[-1] = 0.
    p_noclose_and_ends_cols = np.empty_like(p_mut_given_span)
    p_noclose_and_ends_cols[-1] = 0.
    for j in range(npos):
        j_prev = j - 1
        p_noclose_and_ends_rows[j] = (
                p_ends[j, j:] @ p_noclose_given_ends[j, j:]
                + p_noclose_and_ends_rows[j_prev]
        )
        p_noclose_and_ends_cols[j] = (
                p_ends[:j, j_prev] @ p_noclose_given_ends[:j, j_prev]
                + p_noclose_and_ends_cols[j_prev]
        )
    # Cache index values for repeated slicing.
    range_indexes = np.arange(npos)
    window_indexes = np.minimum(range_indexes, min_gap)
    # Compute the mutation rates given no two mutations are too close
    # one position (j) at a time.
    p_mut_given_span_noclose = np.empty_like(p_mut_given_span)
    for j in range(npos):
        # Numbers of rows and columns of the submatrix for position (j).
        nrows = j + 1
        ncols = npos - j
        # Probability of no close mutations from (a) to (j - inc_gap).
        row_indexes = range_indexes[:nrows]
        p_noclose_5 = p_noclose_given_ends[
            row_indexes, np.maximum(row_indexes, j - inc_gap)
        ]
        # Probability of no mutations from (j - min_gap) to (j - 1).
        p_nomut_window5 = p_nomut_window[
            window_indexes[j::-1], j
        ]
        # Probability of no mutations from (j + 1) to (j + min_gap).
        window3_indexes = window_indexes[:ncols]
        p_nomut_window3 = p_nomut_window[
            window3_indexes, window3_indexes + nrows
        ]
        # Probability of no close mutations from (j + inc_gap) to (b).
        col_indexes = range_indexes[-ncols:]
        p_noclose_3 = p_noclose_given_ends[
            np.minimum(col_indexes, j + inc_gap), col_indexes
        ]
        # Probability of no close mutations in a read with a mutation
        # at position (j): weighted sum over all such reads.
        p_noclose_given_span_mut = np.einsum("ak,bk,ab->k",
                                             p_noclose_5 * p_nomut_window5,
                                             p_noclose_3 * p_nomut_window3,
                                             p_ends[:nrows, -ncols:])
        # Probability of no close mutations in a read that contains
        # position (j): weighted sum over all such reads.
        p_noclose_given_span = (p_noclose_and_ends_rows[j]
                                - p_noclose_and_ends_cols[j])
        # Probability of a mutation at position (j) given no mutations
        # are too close.
        p_mut_given_span_noclose[j] = (p_mut_given_span[j]
                                       * (p_noclose_given_span_mut
                                          / p_noclose_given_span))
    return p_mut_given_span_noclose


def _calc_p_mut_given_span_noclose_from_params(p_mut_given_span: np.ndarray,
                                               p_ends: np.ndarray,
                                               min_gap: int):
    """
    """
    p_nomut_window = _calc_p_nomut_window(p_mut_given_span, min_gap)
    p_noclose_given_ends = _calc_p_noclose_given_ends(p_mut_given_span,
                                                      p_nomut_window)
    return _calc_p_mut_given_span_noclose(p_mut_given_span,
                                          p_ends,
                                          p_noclose_given_ends,
                                          p_nomut_window)


def _calc_p_ends_given_noclose(p_ends: np.ndarray,
                               p_noclose_given_ends: np.ndarray):
    # Find and validate the dimensions.
    npos = p_ends.shape[0]
    if p_ends.shape != (npos, npos):
        raise ValueError(f"Expected p_ends to be shaped ({npos, npos}), "
                         f"but got {p_ends.shape}")
    ncls = p_noclose_given_ends.shape[-1]
    if p_noclose_given_ends.shape != (npos, npos, ncls):
        raise ValueError(
            f"Expected p_noclose_given_ends to be shaped {npos, npos, ncls}, "
            f"but got {p_noclose_given_ends.shape}"
        )
    # Compute and normalize the probability of each coordinate given no
    # two mutations are too close.
    return _normalize(p_noclose_given_ends * np.expand_dims(p_ends, 2))


def _stack2d_p_mut_p_ends(p_mut_given_span: np.ndarray,
                          p_ends: np.ndarray):
    npos, ncls = p_mut_given_span.shape
    stacked = np.empty((npos, npos + ncls))
    stacked[:, :ncls] = p_mut_given_span
    stacked[:, ncls:] = p_ends
    return stacked


def _unstack2d_p_mut_p_ends(stacked: np.ndarray):
    npos, ncol = stacked.shape
    ncls = ncol - npos
    p_mut_given_span = stacked[:, :ncls]
    p_ends = stacked[:, ncls:]
    return p_mut_given_span, p_ends


def _stack3d_p_mut_p_ends(p_mut_given_span: np.ndarray,
                          p_ends: np.ndarray):
    npos, ncls = p_mut_given_span.shape
    stacked = np.empty((npos + 1, npos, ncls))
    stacked[0] = p_mut_given_span
    stacked[1:] = p_ends
    return stacked


def _unstack3d_p_mut_p_ends(stacked: np.ndarray):
    p_mut_given_span = stacked[0]
    p_ends = stacked[1:]
    return p_mut_given_span, p_ends


def _objective(p_mut_p_ends: np.ndarray,
               min_gap: int,
               p_mut_given_span_noclose: np.ndarray,
               p_ends_given_noclose: np.ndarray):
    """
    """
    # Extract the arguments.
    p_mut_given_span, p_ends = _unstack2d_p_mut_p_ends(p_mut_p_ends)
    #
    p_nomut_window = _calc_p_nomut_window(p_mut_given_span, min_gap)
    p_noclose_given_ends = _calc_p_noclose_given_ends(
        p_mut_given_span, p_nomut_window
    )
    #

    #
    p_mut_given_span_noclose_error = _calc_p_mut_given_span_noclose(
        p_mut_given_span,
        p_ends,
        p_noclose_given_ends,
        p_nomut_window
    ) - p_mut_given_span_noclose
    #
    p_ends_given_noclose_error = _calc_p_ends_given_noclose(
        p_ends, p_noclose_given_ends
    ) - p_ends_given_noclose
    objective = _stack3d_p_mut_p_ends(p_mut_given_span_noclose_error,
                                      p_ends_given_noclose_error)
    print("Objective:", np.sum(np.square(objective)))
    return objective


def _objective_p_ends(p_ends: np.ndarray,
                      min_gap: int,
                      p_mut_given_span: np.ndarray,
                      p_ends_given_noclose: np.ndarray):
    """
    """
    p_nomut_window = _calc_p_nomut_window(p_mut_given_span, min_gap)
    p_noclose_given_ends = _calc_p_noclose_given_ends(p_mut_given_span,
                                                      p_nomut_window)
    #
    return _calc_p_ends_given_noclose(
        p_ends, p_noclose_given_ends
    ) - p_ends_given_noclose


def calc_p_mut_p_ends_numpy(p_mut_given_span_noclose: np.ndarray,
                            p_ends_noclose: np.ndarray,
                            min_gap: int,
                            proportions: np.ndarray | None = None,
                            p_mut_given_span_guess: np.ndarray | None = None,
                            p_ends_guess: np.ndarray | None = None, *,
                            f_tol: float = 5e-1,
                            f_rtol: float = 5e-1,
                            x_tol: float = 1e-4,
                            x_rtol: float = 5e-1):
    """ Given observed mutation rates `mus_obs` (which do not include
    any reads that dropped out because they had mutations closer than
    `min_gap` nt apart), estimate the real mutation rates that include
    these unobserved reads.

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
    f_tol: float = 5e-1
        Absolute tolerance in residual.
    f_rtol: float = 5e-1
        Relative tolerance in residual.
    x_tol: float = 1e-4
        Absolute tolerance in step.
    x_rtol: float = 5e-1
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
    # Validate the dimensionality of the arguments.
    if p_mut_given_span_noclose.ndim == 1:
        # If p_mut_given_span_noclose has 1 dimension, then convert it
        # to 2 dimensions, calculate, and convert back to 1 dimension.
        p_mut_given_span, p_ends = calc_p_mut_p_ends_numpy(
            np.expand_dims(p_mut_given_span_noclose, 1),
            p_ends_noclose,
            min_gap,
            proportions,
            p_mut_given_span_guess,
            p_ends_guess,
            f_tol=f_tol,
            f_rtol=f_rtol,
            x_tol=x_tol,
            x_rtol=x_rtol,
        )
        return p_mut_given_span[:, 0], p_ends
    elif p_mut_given_span_noclose.ndim != 2:
        raise ValueError(f"p_mut_given_span_noclose must have 1-2 dimensions, "
                         f"but got {p_mut_given_span_noclose.ndim}")
    if p_ends_noclose.ndim == 2:
        # If p_mut_given_span_noclose has 1 dimension, then convert it
        # to 2 dimensions, calculate, and convert back to 1 dimension.
        return calc_p_mut_p_ends_numpy(
            p_mut_given_span_noclose,
            np.expand_dims(p_ends_noclose, 2),
            min_gap,
            proportions,
            p_mut_given_span_guess,
            p_ends_guess,
            f_tol=f_tol,
            f_rtol=f_rtol,
            x_tol=x_tol,
            x_rtol=x_rtol,
        )
    elif p_ends_noclose.ndim != 3:
        raise ValueError(f"p_ends_noclose must have 2-3 dimensions, "
                         f"but got {p_ends_noclose.ndim}")
    # Validate the dimensions of the arguments.
    npos, ncls = p_mut_given_span_noclose.shape
    if p_ends_noclose.shape != (npos, npos, ncls):
        raise ValueError("The dimensions of p_ends_noclose must be (positions "
                         f"x positions x clusters) {npos, npos, ncls}, "
                         f"but got {p_ends_noclose.shape}")
    if proportions is not None and proportions.shape != (ncls,):
        raise ValueError("If given, proportions must have dimensions "
                         f"{ncls,}, but got {proportions.shape}")
    if (p_mut_given_span_guess is not None
            and (p_mut_given_span_guess.shape != (npos, ncls))):
        raise ValueError("If given, p_mut_given_span_guess must have "
                         "the same dimensions as p_mut_given_span_noclose, "
                         f"but got {p_mut_given_span_guess.shape} "
                         f"(≠ {p_mut_given_span_noclose.shape})")
    if p_ends_guess is not None and p_ends_guess.shape != (npos, npos):
        raise ValueError("If given, p_ends_guess must have "
                         "the same dimensions as p_ends_noclose, "
                         f"but got {p_ends_guess.shape} "
                         f"(≠ {p_ends_noclose.shape})")
    # Validate and normalize the values.
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    min_gap = _adjust_min_gap(npos, min_gap)
    p_mut_given_span_noclose = _clip(p_mut_given_span_noclose)
    p_ends_noclose = _normalize(p_ends_noclose)
    if min_gap == 0:
        # No two mutations can be too close.
        return p_mut_given_span_noclose, p_ends_noclose
    # Determine the initial guess of the real mutation rates.
    if p_mut_given_span_guess is None or p_ends_guess is None:
        # If either guess is missing, then use the argument values.
        initial = _stack2d_p_mut_p_ends(p_mut_given_span_noclose,
                                        np.average(p_ends_noclose,
                                                   axis=2,
                                                   weights=proportions))
        if p_mut_given_span_guess is not None:
            logger.warning("Got guess for p_mut_given_span but not p_ends")
        if p_ends_guess is not None:
            logger.warning("Got guess for p_ends but not p_mut_given_span")
    else:
        # Use the explicit initial guess.
        initial = _stack2d_p_mut_p_ends(_clip(p_mut_given_span_guess),
                                        _normalize(p_ends_guess))
    # Solve for the "real" mutation rates that yield minimal difference
    # between the mutation rates that would be expected to be observed
    # given the real mutation rates (i.e. when reads with mutations too
    # close are removed from the real mutation rates) and the actual
    # observed mutation rates. Use Newton's method, which finds the
    # parameters of a function that make it evaluate to zero, with the
    # Krylov approximation of the Jacobian, which improves performance.
    p_mut_given_span, p_ends = _unstack2d_p_mut_p_ends(newton_krylov(
        lambda p_mut_p_ends: _objective(p_mut_p_ends,
                                        min_gap,
                                        p_mut_given_span_noclose,
                                        p_ends_noclose),
        initial,
        f_tol=f_tol,
        f_rtol=f_rtol,
        x_tol=x_tol,
        x_rtol=x_rtol,
    ))
    return _clip(p_mut_given_span), _normalize(p_ends)


def calc_p_mut_given_span(p_mut_given_span_noclose: np.ndarray,
                          min_gap: int,
                          p_ends: np.ndarray,
                          initial: np.ndarray | None = None, *,
                          f_tol: float = 5e-1,
                          f_rtol: float = 5e-1,
                          x_tol: float = 1e-4,
                          x_rtol: float = 5e-1):
    """ Given observed mutation rates `mus_obs` (which do not include
    any reads that dropped out because they had mutations closer than
    `min_gap` nt apart), estimate the real mutation rates that include
    these unobserved reads.

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
    f_tol: float = 5e-1
        Absolute tolerance in residual.
    f_rtol: float = 5e-1
        Relative tolerance in residual.
    x_tol: float = 1e-4
        Absolute tolerance in step.
    x_rtol: float = 5e-1
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
    # Validate the dimensionality of the arguments.
    if p_mut_given_span_noclose.ndim == 1:
        # If p_mut_given_span_noclose has 1 dimension, then convert it
        # to 2 dimensions, calculate, and convert back to 1 dimension.
        return calc_p_mut_given_span(
            np.expand_dims(p_mut_given_span_noclose, 1),
            min_gap,
            p_ends,
            initial,
            f_tol=f_tol,
            f_rtol=f_rtol,
            x_tol=x_tol,
            x_rtol=x_rtol,
        )[:, 0]
    elif p_mut_given_span_noclose.ndim != 2:
        raise ValueError(f"p_mut_given_span_noclose must have 1-2 dimensions, "
                         f"but got {p_mut_given_span_noclose.ndim}")
    # Validate the dimensions of the arguments.
    if initial is None:
        # If no initial value is given, then use the mutation rate among
        # reads with no two mutations too close.
        initial = p_mut_given_span_noclose
    elif initial.shape != p_mut_given_span_noclose.shape:
        raise ValueError("If given, initial must have the same dimensions as "
                         f"p_mut_given_span_noclose, but got {initial.shape} "
                         f"≠ {p_mut_given_span_noclose.shape}")
    # Validate and normalize the values.
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    min_gap = _adjust_min_gap(p_mut_given_span_noclose.shape[0], min_gap)
    p_mut_given_span_noclose = _clip(p_mut_given_span_noclose)
    p_ends = _normalize(p_ends)
    if min_gap == 0:
        # No two mutations can be too close.
        return p_mut_given_span_noclose

    # Find the mutation rates including reads with two mutations too
    # close by using the Newton-Krylov method to make the difference
    # between the observed and the theoretical mutation rates equal 0.

    def objective(p_mut_given_span: np.ndarray):
        return _calc_p_mut_given_span_noclose_from_params(
            p_mut_given_span,
            p_ends,
            min_gap
        ) - p_mut_given_span_noclose

    return _clip(newton_krylov(objective,
                               initial,
                               f_tol=f_tol,
                               f_rtol=f_rtol,
                               x_tol=x_tol,
                               x_rtol=x_rtol))


def calc_prop_adj_numpy(prop_obs: np.ndarray, f_obs: np.ndarray):
    """ Calculate the adjusted proportion of the clusters given their
    observed proportions and the observer bias. """
    weighted_prop_obs = prop_obs / f_obs
    return weighted_prop_obs / np.sum(weighted_prop_obs)

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
