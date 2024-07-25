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

import warnings
from logging import getLogger
from typing import Iterable

import numpy as np
from numba import jit, NumbaPerformanceWarning

from ..array import find_dims, triangular

logger = getLogger(__name__)

# Disable performance warnings from Numba.
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

# Define dimensions
READS = "reads"
POSITIONS = "positions"
CLUSTERS = "clusters"


@jit()
def _clip(x: np.ndarray):
    """ Fill NaN with 0, infinity with 1, and restrict all values to the
    interval [0, 1].

    Parameters
    ----------
    x: np.ndarray
        Values to fill and clip.

    Returns
    -------
    np.ndarray
        Array of the same shape as `x` with all values ≥ 0 and ≤ 1.
    """
    return np.clip(np.nan_to_num(x), 0., 1.)


@jit()
def _normalize(x: np.ndarray):
    """ Normalize the values to sum to 1, or if they sum to 0, then
    return an array with the same value for each element. """
    x_sum = np.sum(x)
    if x_sum == 0.:
        if x.size == 0:
            # Handle the edge case of size zero, which would raise a
            # ZeroDivisionError if handled in the next branch.
            return np.ones_like(x)
        # If the sum of the input array is 0 and the array has at least
        # 1 element, then return an array of the same size as the input
        # where all elements are equal and sum to 1.
        return np.full_like(x, 1. / x.size)
    # Divide each element by the sum of all elements so the resulting
    # array has the same size as the input and all elements sum to 1.
    return x / x_sum


@jit()
def _triu_log(a: np.ndarray):
    """ Calculate the logarithm of the upper triangle(s) of array `a`.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` has at least 2 dimensions.
    -   The first two dimensions of `a` have equal length.

    Parameters
    ----------
    a: np.ndarray
        Array of whose upper triangle to compute the logarithm.

    Returns
    -------
    np.ndarray
        Logarithm of the upper triangle(s) of `a`.
    """
    log = np.empty_like(a)
    for j in range(a.shape[0]):
        log[j, j:] = np.log(a[j, j:])
    return log


def triu_log(a: np.ndarray):
    """ Calculate the logarithm of the upper triangle(s) of array `a`.
    In the result, elements below the main diagonal are undefined.

    Parameters
    ----------
    a: np.ndarray
        Array (≥ 2 dimensions) of whose upper triangle to compute the
        logarithm; the first 2 dimensions must have equal lengths.

    Returns
    -------
    np.ndarray
        Logarithm of the upper triangle(s) of `a`.
    """
    if a.ndim < 2:
        raise ValueError(f"a must have ≥ 2 dimensions, but got {a.ndim}")
    npos = a.shape[0]
    if a.shape[:2] != (npos, npos):
        raise ValueError("The first 2 dimensions of a must have equal lengths, "
                         f"but got dimensions {a.shape}")
    return _triu_log(a)


@jit()
def _triu_sum(a: np.ndarray):
    """ Calculate the sum over the upper triangle(s) of array `a`.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` has at least 2 dimensions.
    -   The first two dimensions of `a` have equal length.

    Parameters
    ----------
    a: np.ndarray
        Array whose upper triangle to sum.

    Returns
    -------
    np.ndarray
        Sum of the upper triangle(s), with the same shape as the third
        and subsequent dimensions of `a`.
    """
    a_sum = np.zeros(a.shape[2:])
    # Sum over axes 0 and 1.
    for j in range(a.shape[0]):
        a_sum += a[j, j:].sum(axis=0)
    return a_sum


@jit()
def _triu_cumsum(a: np.ndarray):
    """ Calculate the cumulative sum from the upper right corner of `a`
    to every other element in the upper right triangle.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` has at least 2 dimensions.
    -   The first two dimensions of `a` have equal length.

    Parameters
    ----------
    a: np.ndarray
        Array whose upper triangle to sum.

    Returns
    -------
    np.ndarray
        Cumulative sum, with the same shape as `a`.
    """
    a_cumsum = np.empty_like(a)
    if a_cumsum.size == 0:
        return a_cumsum
    npos = a.shape[0]
    extras = a.shape[2:]
    # Numba does not (as of v0.59.0) support np.cumsum with the optional
    # argument "axis", so need to replicate this functionality manually.
    # First, fill the first row, starting from the top right element.
    a_cumsum[0, -1] = a[0, -1]
    for col in range(npos - 2, -1, -1):
        a_cumsum[0, col] = a[0, col] + a_cumsum[0, col + 1]
    # Then, fill each subsequent row.
    for row in range(1, npos):
        row_cumsum = np.zeros(extras, dtype=a.dtype)
        for col in range(npos - 1, row - 1, -1):
            row_cumsum += a[row, col]
            a_cumsum[row, col] = row_cumsum + a_cumsum[row - 1, col]
    return a_cumsum


@jit()
def _triu_div(numer: np.ndarray, denom: np.ndarray):
    """ Divide the upper triangles of `numer` and `denom`.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `numer` has at least 2 dimensions.
    -   The first two dimensions of `numer` have equal length.
    -   `denom` has the same first 2 dimensions as `numer`.
    -   `numer` and `denom` can be broadcast to each other.
    -   No value in the upper triangle of `denom` is 0.

    Parameters
    ----------
    numer: np.ndarray
        Numerator of division.
    denom: np.ndarray
        Denominator of division.

    Returns
    -------
    np.ndarray
        Quotient of the upper triangles; values below the main diagonal
        are undefined.
    """
    quotient = np.empty(np.broadcast_shapes(numer.shape, denom.shape))
    for j in range(numer.shape[0]):
        quotient[j, j:] = numer[j, j:] / denom[j, j:]
    return quotient


# @jit()
def _triu_norm(a: np.ndarray):
    """ Normalize the upper triangle of array `a` to sum to 1.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` has at least 2 dimensions.
    -   The first two dimensions of `a` have equal length.
    -   The elements of the upper triangle of `a` do not sum to 0.

    Parameters
    ----------
    a: np.ndarray
        Array to normalize.

    Returns
    -------
    np.ndarray
        Array of the same shape as `a` but scaled so that the
    """
    # Calculate the sum over axes 0 and 1.
    a_sum = _triu_sum(a)
    # Normalize by dividing by that sum, ignoring division by zero,
    # which is handled subsequently.
    a_norm = _triu_div(a, np.broadcast_to(a_sum, a.shape))
    # Determine which sums to normalize by are zero.
    a_sum_zero = a_sum == 0.
    if np.any(np.atleast_1d(a_sum_zero)):
        # Count the elements over which normalization occurred.
        n_norm = triangular(a.shape[0])
        if n_norm == 0:
            # Handle the edge case of size zero, which would raise a
            # ZeroDivisionError if handled in the next branch.
            return np.ones_like(a_norm)
        # For each sum that was zero, set all elements to the same value
        # such that they sum to 1.
        fill = 1. / n_norm
        if a_sum_zero.ndim:
            for indexes in zip(*np.nonzero(a_sum_zero), strict=True):
                a_norm[(slice(None),) * 2 + indexes] = fill
        else:
            a_norm[:, :] = fill
    return a_norm


@jit()
def _triu_mul(factor1: np.ndarray, factor2: np.ndarray):
    """ Multiply the upper triangles of `numer` and `denom`.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `factor1` has at least 2 dimensions.
    -   The first two dimensions of `factor1` have equal length.
    -   `factor2` has the same first 2 dimensions as `factor1`.
    -   `factor1` and `factor2` can be broadcast to each other.

    Parameters
    ----------
    factor1: np.ndarray
        Factor 1.
    factor2: np.ndarray
        Factor 2.

    Returns
    -------
    np.ndarray
        Product of the upper triangles; values below the main diagonal
        are undefined.
    """
    product = np.empty(np.broadcast_shapes(factor1.shape, factor2.shape))
    for j in range(factor1.shape[0]):
        product[j, j:] = factor1[j, j:] * factor2[j, j:]
    return product


@jit()
def _triu_dot(a: np.ndarray, b: np.ndarray):
    """ Dot product of `a` and `b` over their first 2 dimensions.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` and `b` both have at least 2 dimensions.
    -   The first and second dimensions of `a` are equal.
    -   The first and second dimensions of `b` are equal.

    Parameters
    ----------
    a: np.ndarray
        Array 1.
    b: np.ndarray
        Array 2.

    Returns
    -------
    np.ndarray
        Dot product of `a` and `b` over their first 2 dimensions.
    """
    dot = np.zeros(np.broadcast_shapes(a.shape, b.shape)[2:])
    # Dot product over axes 0 and 1.
    for j in range(a.shape[0]):
        dot += np.sum(a[j, j:] * b[j, j:], axis=0)
    return dot


@jit()
def _triu_allclose(a: np.ndarray,
                   b: np.ndarray,
                   rtol: float = 1.e-3,
                   atol: float = 1.e-6):
    """ Whether the upper triangles of `a` and `b` are all close.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` and `b` both have at least 2 dimensions.
    -   The first and second dimensions of `a` are equal.
    -   The first and second dimensions of `b` are equal.

    Parameters
    ----------
    a: np.ndarray
        Array 1.
    b: np.ndarray
        Array 2.
    rtol: float = 1.0e-3
        Relative tolerance.
    atol: float = 1.0e-6
        Absolute tolerance.

    Returns
    -------
    bool
        Whether all elements of the upper triangles of `a` and `b` are
        close using the function `np.allclose`.
    """
    for j in range(a.shape[0]):
        if not np.allclose(a[j, j:], b[j, j:], rtol=rtol, atol=atol):
            return False
    return True


@jit()
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


@jit()
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


@jit()
def _calc_p_noclose_given_ends(p_mut_given_span: np.ndarray,
                               p_nomut_window: np.ndarray):
    """ Given underlying mutation rates (`p_mut_given_span`), calculate
    the probability that a read starting at position (a) and ending at
    position (b) would have no two mutations too close, for each (a) and
    (b) where 1 ≤ a ≤ b ≤ L (biological coordinates) or 0 ≤ a ≤ b < L
    (Python coordinates).

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
    # Fill the main diagonal and the (min_gap) diagonals below the main
    # diagonal with ones.
    for j in range(npos):
        p_noclose_given_ends[j, max(j - min_gap, 0): (j + 1)] = 1.
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


@jit()
def _calc_rectangular_sum(array: np.ndarray):
    """ For each element of the main diagonal, calculate the sum over
    the rectangular array from that element to the upper right corner.
    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `array` has at least 2 dimensions.
    -   The first and second dimensions of `array` have equal lengths.

    Parameters
    ----------
    array: np.ndarray
        Array of at least two dimensions for which to calculate the sum
        of each rectangular array from each element on the main diagonal
        to the upper right corner.

    Returns
    -------
    np.ndarray
        Array with all but the first dimension of `array` indicating the
        sum of the array from each element on the main diagonal to the
        upper right corner of `array`.
    """
    npos = array.shape[0]
    dims = array.shape[1:]
    # For each position, calculate the total area of all rows up to and
    # including that position and of all columns up to but not including
    # that position.
    rows_area = np.empty(dims)
    cols_area = np.empty(dims)
    if npos > 0:
        # Initialize the first element of each array.
        rows_area[0] = array[0].sum(axis=0)
        cols_area[0] = 0.
        # Calculate every successive element.
        for j in range(1, npos):
            i = j - 1
            rows_area[j] = rows_area[i] + array[j, j:].sum(axis=0)
            cols_area[j] = cols_area[i] + array[:j, i].sum(axis=0)
    # For each position, the area of the rectangular array from that
    # position to the upper right corner equals the area of the rows up
    # to and including that position minus the area of the columns up to
    # but not including that position.
    return rows_area - cols_area


@jit()
def _calc_p_mut_given_span_noclose(p_mut_given_span: np.ndarray,
                                   p_ends: np.ndarray,
                                   p_noclose_given_ends: np.ndarray,
                                   p_nomut_window: np.ndarray):
    """ Calculate the mutation rates of only reads with no two mutations
    too close. """
    # Find the dimensions.
    npos, ncls = p_mut_given_span.shape
    inc_gap = p_nomut_window.shape[0]
    min_gap = inc_gap - 1
    if min_gap == 0:
        # If min_gap == 0, then no mutations can be too close, so return
        # the original mutation rates.
        return p_mut_given_span
    # Calculate the probability that a read spanning each position would
    # have no two mutations too close.
    p_noclose_given_span = _calc_rectangular_sum(p_noclose_given_ends
                                                 * p_ends[:, :, np.newaxis])
    # Compute the mutation rates given no two mutations are too close
    # one position (j) at a time.
    p_mut_given_span_noclose = np.nan_to_num(
        p_mut_given_span / p_noclose_given_span
    )
    for j in range(npos):
        nrows = j + 1
        ncols = npos - j
        # For each starting position (a), calculate the probability of
        # no two mutations too close from (a) to (j).
        p_noclose_a = np.empty((ncls, nrows))
        for a in range(nrows):
            # Number of bases 5' of (j) that must have no mutations:
            # the smallest of (min_gap) and (j - a).
            nomut = min(min_gap, j - a)
            # Probability of no close mutations, (a) to (j - 1 - nomut):
            #   p_noclose_given_ends[a, max(j - 1 - nomut, a)]
            # Probability of no mutations over the (nomut)-sized window
            # from (j - nomut) to (j - 1).
            #   p_nomut_window[nomut, j]
            # Probability of no close mutations from (a) to (j).
            p_noclose_a[:, a] = (
                    p_noclose_given_ends[a, max(j - 1 - nomut, a)]
                    * p_nomut_window[nomut, j]
            )
        # For each ending position (b), calculate the probability of
        # no two mutations too close from (j) to (b).
        p_noclose_b = np.empty((ncls, ncols))
        for b in range(j, npos):
            # Number of bases 3' of (j) that must have no mutations:
            # the smallest of (min_gap) and (b - j).
            nomut = min(min_gap, b - j)
            # Probability of no mutations over the (nomut)-sized window
            # from (j + 1) to (j + nomut):
            #   p_nomut_window[nomut, j + 1 + nomut]
            # Probability of no close mutations, (j + 1 + nomut) to (b):
            #   p_noclose_given_ends[min(j + 1 + nomut, b), b]
            # Probability of no close mutations from (j) to (b).
            p_noclose_b[:, b - j] = (
                    p_nomut_window[nomut, nrows + nomut]
                    * p_noclose_given_ends[min(nrows + nomut, b), b]
            )
        # The probability that a read has a mutation at position (j)
        # given that it has no two mutations too close is the weighted
        # sum of such probabilities for all (a) and (b).
        for k in range(ncls):
            p_mut_given_span_noclose[j, k] *= (p_noclose_a[k]
                                               @ p_ends[:nrows, -ncols:]
                                               @ p_noclose_b[k])
    return p_mut_given_span_noclose


def _slice_p_ends(p_ends: np.ndarray,
                  p_ends_cumsum: np.ndarray,
                  end5: int,
                  end3: int):
    """ Slice a matrix of end coordinate probabilities. """
    p_ends_slice = p_ends[end5: end3, end5: end3].copy()
    if p_ends_slice.size > 0:
        p_ends_slice[0, :-1] = p_ends[:end5 + 1, end5: end3 - 1].sum(axis=0)
        p_ends_slice[1:, -1] = p_ends[end5 + 1: end3, end3 - 1:].sum(axis=1)
        p_ends_slice[0, -1] = p_ends_cumsum[end5, end3 - 1]
    return p_ends_slice


def _find_split_positions(p_mut: np.ndarray, min_gap: int, threshold: float):
    """ Find the positions (at or below the threshold) at which to split
    the mutation rates. """
    npos, ncls = p_mut.shape
    min_gap = _adjust_min_gap(npos, min_gap)
    if min_gap == 0 or ncls == 0:
        return np.array([], dtype=int)
    # Count the cumulative number of mutation rates below the threshold.
    cum_below = np.cumsum(np.concatenate([np.zeros(min_gap, dtype=bool),
                                          p_mut.max(axis=1) <= threshold]))
    # Label every position for which it and the previous (min_gap - 1)
    # positions are all below the threshold.
    win_below = cum_below[min_gap:] - cum_below[:-min_gap] == min_gap
    # Locate each stretch of consecutive positions where the previous
    # (min_gap - 1) positions are all below the threshold; the first
    # (min_gap - 1) positions can be ignored because they cannot have
    # (min_gap - 1) preceding positions.
    diff_win_below = np.diff(win_below.astype(int)[(min_gap - 1):])
    first_pos = np.flatnonzero(diff_win_below == 1) + 1
    last_pos = np.flatnonzero(diff_win_below == -1) + min_gap
    # Collect all boundary positions and find the unique ones.
    return np.unique(np.concatenate([first_pos, last_pos]))


def _split_positions(p_mut: np.ndarray,
                     p_mut_init: np.ndarray,
                     p_ends: np.ndarray,
                     split_pos: np.ndarray):
    """ Split `p_mut` and `p_ends` where at least `min_gap` consecutive
    positions are no greater than `threshold` in every cluster. """
    dims = find_dims([(POSITIONS, CLUSTERS),
                      (POSITIONS, CLUSTERS),
                      (POSITIONS, POSITIONS)],
                     [p_mut, p_mut_init, p_ends],
                     ["p_mut", "p_mut_init", "p_ends"])
    n_splits, = split_pos.shape
    if n_splits == 0:
        return [p_mut], [p_mut_init], [p_ends]
    npos = dims[POSITIONS]
    p_mut_split = np.split(p_mut, split_pos)
    p_mut_init_split = np.split(p_mut_init, split_pos)
    p_ends_cumsum = _triu_cumsum(p_ends)
    p_ends_split = [_slice_p_ends(p_ends, p_ends_cumsum, 0, split_pos[0])]
    for i in range(n_splits - 1):
        p_ends_split.append(_slice_p_ends(p_ends,
                                          p_ends_cumsum,
                                          split_pos[i],
                                          split_pos[i + 1]))
    p_ends_split.append(_slice_p_ends(p_ends,
                                      p_ends_cumsum,
                                      split_pos[-1],
                                      npos))
    return p_mut_split, p_mut_init_split, p_ends_split


def _calc_p_mut_given_span(p_mut_given_span_observed: np.ndarray,
                           min_gap: int,
                           p_ends: np.ndarray,
                           init_p_mut_given_span: np.ndarray, *,
                           quick_unbias: bool = True,
                           quick_unbias_thresh: float = 0.,
                           f_tol: float = 1.e-4,
                           x_rtol: float = 1.e-3):
    """ Calculate the underlying mutation rates including for reads with
    two mutations too close based on the observed mutation rates. """
    # Validate the dimensionality of the arguments.
    if p_mut_given_span_observed.ndim == 1:
        # If p_mut_given_span_noclose has 1 dimension, then convert it
        # to 2 dimensions, calculate, and convert back to 1 dimension.
        return _calc_p_mut_given_span(
            p_mut_given_span_observed[:, np.newaxis],
            min_gap,
            p_ends,
            init_p_mut_given_span,
            f_tol=f_tol,
            x_rtol=x_rtol,
        )[:, 0]
    dims = find_dims([(POSITIONS, CLUSTERS)],
                     [p_mut_given_span_observed],
                     ["p_mut_given_span_observed"],
                     nonzero=True)
    min_gap = _adjust_min_gap(dims[POSITIONS], min_gap)
    if min_gap == 0:
        # No two mutations can be too close.
        return p_mut_given_span_observed

    if quick_unbias:
        # Split the mutation rates and end coordinate distributions,
        # calculate them separately, and reassemble them.
        p_mut_given_span_list = list()
        for p_mut_split, p_mut_init_split, p_ends_split in zip(
                *_split_positions(p_mut_given_span_observed,
                                  init_p_mut_given_span,
                                  p_ends,
                                  _find_split_positions(
                                      p_mut_given_span_observed,
                                      min_gap,
                                      quick_unbias_thresh
                                  )),
                strict=True
        ):
            if p_mut_split.max(initial=0.) > quick_unbias_thresh:
                p_mut_given_span_list.append(_calc_p_mut_given_span(
                    p_mut_split,
                    min_gap,
                    p_ends_split,
                    p_mut_init_split,
                    quick_unbias=False,
                    quick_unbias_thresh=quick_unbias_thresh,
                    f_tol=f_tol,
                    x_rtol=x_rtol,
                ))
            else:
                # No mutation rate exceeds the threshold.
                p_mut_given_span_list.append(p_mut_split)
        return (np.concatenate(p_mut_given_span_list, axis=0)
                if p_mut_given_span_list
                else np.empty((dims[POSITIONS], dims[CLUSTERS])))

    # Use the Newton-Krylov method to solve for the total mutation rates
    # (including reads with two mutations too close) that result in zero
    # difference between theoretically observed and actually observed
    # mutation rates.

    def objective(p_mut_given_span: np.ndarray):
        """ Calculate the difference between the observed mutation rates
        and the mutation rates of reads with no two mutations too close
        assuming the true mutaiton rates are `p_mut_given_span`. """
        p_nomut_window = _calc_p_nomut_window(p_mut_given_span, min_gap)
        p_noclose_given_ends = _calc_p_noclose_given_ends(p_mut_given_span,
                                                          p_nomut_window)
        return p_mut_given_span_observed - _calc_p_mut_given_span_noclose(
            p_mut_given_span,
            p_ends,
            p_noclose_given_ends,
            p_nomut_window
        )

    # Import scipy here instead of at the top of this module because
    # its import is slow enough to impact global startup time.
    from scipy.optimize import newton_krylov, NoConvergence

    try:
        return _clip(newton_krylov(objective,
                                   init_p_mut_given_span,
                                   f_tol=f_tol,
                                   x_rtol=x_rtol))
    except NoConvergence as error:
        logger.warning(f"Failed to unbias mutation rates: {error}")
        return init_p_mut_given_span


def _calc_p_ends(p_ends_observed: np.ndarray,
                 p_noclose_given_ends: np.ndarray,
                 p_mut_given_span: np.ndarray,
                 p_clust: np.ndarray):
    """ Calculate the proportion of total reads with each pair of 5' and
    3' coordinates.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   Every value in the upper triangle of `p_ends_observed` is
        ≥ 0 and ≤ 1; no values below the main diagonal are used.
    -   The upper triangle of `p_ends_observed` sums to 1.
    -   `min_gap` is a non-negative integer.
    -   Every value in `p_mut_given_span` is ≥ 0 and ≤ 1.

    Parameters
    ----------
    p_ends_observed: np.ndarray
        3D (positions x positions x clusters) array of the proportion of
        observed reads in each cluster beginning at the row position and
        ending at the column position.
    p_noclose_given_ends: np.ndarray
        3D (positions x positions x clusters) array of the pobabilities
        that a read with 5' and 3' coordinates corresponding to the row
        and column would have no two mutations too close.
    p_mut_given_span: np.ndarray
        2D (positions x clusters) array of the total mutation rate at
        each position in each cluster.
    p_clust: np.ndarray
        1D (clusters) array of the proportion of each cluster.

    Returns
    -------
    np.ndarray
        2D (positions x positions) array of the proportion of reads
        beginning at the row position and ending at the column position.
        This array is assumed to be identical for all clusters.
    """
    # Validate the dimensionality of the arguments.
    if p_ends_observed.ndim == 2:
        # If p_ends_observed has 2 dimensions, then promote it to 3
        # dimensions first.
        return _calc_p_ends(p_ends_observed[:, :, np.newaxis],
                            p_noclose_given_ends,
                            p_mut_given_span,
                            p_clust)
    dims = find_dims([(POSITIONS, POSITIONS, CLUSTERS),
                      (POSITIONS, POSITIONS, CLUSTERS),
                      (POSITIONS, CLUSTERS),
                      (CLUSTERS,)],
                     [p_ends_observed,
                      p_noclose_given_ends,
                      p_mut_given_span,
                      p_clust],
                     ["p_ends_observed",
                      "p_noclose_given_ends",
                      "p_mut_given_span",
                      "p_clust"],
                     nonzero=True)
    # Calculate the proportion of total reads that would have each
    # pair of end coordinates.
    p_ends = _triu_norm(_triu_div(p_ends_observed, p_noclose_given_ends))
    # Return a consensus distribution among all clusters.
    if dims[CLUSTERS] == 1:
        return p_ends[:, :, 0]
    return np.average(p_ends, axis=2, weights=p_clust)


def calc_p_noclose(p_ends: np.ndarray,
                   p_noclose_given_ends: np.ndarray):
    """ Calculate the proportion of each cluster considering only reads
    with no two mutations too close. """
    # Validate the dimensionality of the arguments.
    find_dims([(POSITIONS, POSITIONS), (POSITIONS, POSITIONS, CLUSTERS)],
              [p_ends, p_noclose_given_ends],
              ["p_ends", "p_noclose_given_ends"],
              nonzero=True)
    # Compute the weighted sum of the probabilities that reads from each
    # cluster would have no two mutations too close.
    return _triu_dot(p_ends[:, :, np.newaxis], p_noclose_given_ends)


def calc_p_clust(p_clust_observed: np.ndarray, p_noclose: np.ndarray):
    """ Cluster proportion among all reads.

    Parameters
    ----------
    p_clust_observed: np.ndarray
        Proportion of each cluster among reads with no two mutations too
        close.
        1D (clusters)
    p_noclose: np.ndarray
        Probability that a read from each cluster would have no two
        mutations too close.
        1D (clusters)

    Returns
    -------
    np.ndarray
        Proportion of each cluster among all reads.
        1D (clusters)
    """
    # Validate the dimensions.
    find_dims([(CLUSTERS,), (CLUSTERS,)],
              [p_clust_observed, p_noclose],
              ["p_clust_observed", "p_noclose"],
              nonzero=True)
    # The cluster proportions among all reads are obtained by weighting
    # each cluster proportion among reads with no two mutations too
    # close by the reciprocal of the probability that no two mutations
    # are too close in that cluster, then normalizing so the sum is 1.
    return _normalize(p_clust_observed / p_noclose)


def calc_p_clust_given_noclose(p_clust: np.ndarray, p_noclose: np.ndarray):
    """ Cluster proportions among reads with no two mutations too close.

    Parameters
    ----------
    p_clust: np.ndarray
        Proportion of each cluster among all reads.
        1D (clusters)
    p_noclose: np.ndarray
        Probability that a read from each cluster would have no two
        mutations too close.
        1D (clusters)

    Returns
    -------
    np.ndarray
        Proportion of each cluster among reads with no two mutations too
        close.
        1D (clusters)
    """
    # Validate the dimensions.
    find_dims([(CLUSTERS,), (CLUSTERS,)],
              [p_clust, p_noclose],
              ["p_clust", "p_noclose"],
              nonzero=True)
    # The cluster proportions among reads with no two mutations too
    # close are obtained by weighting each cluster proportion by the
    # probability that no two mutations are too close in that cluster,
    # then normalizing so the sum is 1.
    return _normalize(p_clust * p_noclose)


def calc_p_noclose_given_ends(p_mut_given_span: np.ndarray, min_gap: int):
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
    if p_mut_given_span.ndim == 2:
        p_mut_given_span = _clip(p_mut_given_span)
        return _calc_p_noclose_given_ends(p_mut_given_span,
                                          _calc_p_nomut_window(p_mut_given_span,
                                                               min_gap))
    if p_mut_given_span.ndim == 1:
        return calc_p_noclose_given_ends(
            p_mut_given_span[:, np.newaxis], min_gap
        ).reshape((p_mut_given_span.size, p_mut_given_span.size))
    raise ValueError("p_mut_given_span must have 1 or 2 dimensions, "
                     f"but got {p_mut_given_span.ndim}")


def calc_rectangluar_sum(array: np.ndarray):
    """ For each element of the main diagonal, calculate the sum over
    the rectangular array from that element to the upper right corner.

    Parameters
    ----------
    array: np.ndarray
        Array of at least two dimensions for which to calculate the sum
        of each rectangular array from each element on the main diagonal
        to the upper right corner. The first two dimensions must have
        equal lengths.

    Returns
    -------
    np.ndarray
        Array with all but the first dimension of `array` indicating the
        sum of the array from each element on the main diagonal to the
        upper right corner of `array`.
    """
    # Validate the dimensions.
    find_dims([(POSITIONS, POSITIONS, None)], [array])
    return _calc_rectangular_sum(array)


def calc_params(p_mut_given_span_observed: np.ndarray,
                p_ends_observed: np.ndarray,
                p_clust_observed: np.ndarray,
                min_gap: int,
                guess_p_mut_given_span: np.ndarray | None = None,
                guess_p_ends: np.ndarray | None = None,
                guess_p_clust: np.ndarray | None = None, *,
                prenormalize: bool = True,
                max_iter: int = 128,
                convergence_thresh: float = 1.e-4,
                **kwargs):
    """ Calculate the three sets of parameters based on observed data.

    Parameters
    ----------
    p_mut_given_span_observed: np.ndarray
        Observed probability that each position is mutated given that no
        two mutations are too close:
        2D array (positions x clusters)
    p_ends_observed: np.ndarray
        Observed proportion of reads aligned with each pair of 5' and 3'
        end coordinates given that no two mutations are too close:
        3D array (positions x positions x clusters)
    p_clust_observed: np.ndarray
        Observed proportion of reads in each cluster given that no two
        mutations are too close:
        1D array (clusters)
    min_gap: int
        Minimum number of non-mutated bases between two mutations. Must
        be a non-negative integer.
    guess_p_mut_given_span: np.ndarray | None = None
        Initial guess for the probability that each position is mutated.
        If given, must be a 2D array (positions x clusters); defaults to
        `p_mut_given_span_observed`.
    guess_p_ends: np.ndarray | None = None
        Initial guess for the proportion of total reads aligned to each
        pair of 5' and 3' end coordinates. If given, must be a 2D array
        (positions x positions); defaults to `p_ends_observed`.
    guess_p_clust: np.ndarray | None = None
        Initial guess for the proportion of total reads in each cluster.
        If given, must be a 1D array (clusters); defaults to
        `p_clust_observed`.
    prenormalize: bool = True
        Fill missing values in `guess_p_mut_given_span`, `guess_p_ends`,
        and `guess_p_clust`, and clip every value to be ≥ 0 and ≤ 1.
        Ensure the proportions in `guess_p_clust` and the upper triangle
        of `guess_p_ends` sum to 1.
    max_iter: int = 128
        Maximum number of iterations in which to refine the parameters.
    convergence_thresh: float = 1.e-4
        Convergence threshold based on the root-mean-square difference
        in mutation rates between consecutive iterations.
    **kwargs
        Additional keyword arguments for `_calc_p_mut_given_span`.
    """
    # Validate the dimensions.
    dims = find_dims([(POSITIONS, CLUSTERS),
                      (POSITIONS, POSITIONS, CLUSTERS),
                      (CLUSTERS,)],
                     [p_mut_given_span_observed,
                      p_ends_observed,
                      p_clust_observed],
                     ["p_mut_given_span_observed",
                      "p_ends_observed",
                      "p_clust_observed"],
                     nonzero=True)
    # Normalize the values.
    min_gap = _adjust_min_gap(dims[POSITIONS], min_gap)
    if prenormalize:
        p_mut_given_span_observed = _clip(p_mut_given_span_observed)
        p_ends_observed = _triu_norm(_clip(p_ends_observed))
        p_clust_observed = _normalize(_clip(p_clust_observed))
    # Determine the initial guess for the mutation rates.
    if guess_p_mut_given_span is None:
        # If no initial guess was given, then use the mutation rates of
        # the observed reads.
        guess_p_mut_given_span = p_mut_given_span_observed
    elif guess_p_mut_given_span.shape != p_mut_given_span_observed.shape:
        raise ValueError("If given, guess_p_mut_given_span must have the same "
                         "dimensions as p_mut_given_span_observed, but got "
                         f"{guess_p_mut_given_span.shape} "
                         f"≠ {p_mut_given_span_observed.shape}")
    elif prenormalize:
        # Ensure the initial mutation rates are clipped.
        guess_p_mut_given_span = _clip(guess_p_mut_given_span)
    # Determine the initial guess for the cluster proportions.
    if guess_p_clust is None:
        # If no initial guess was given, then use the proportions of the
        # observed reads.
        guess_p_clust = p_clust_observed
    elif guess_p_clust.shape != p_clust_observed.shape:
        raise ValueError(
            "If given, guess_p_clust must have the same dimensions as "
            f"p_clust_observed, but got {guess_p_clust.shape} "
            f"≠ {p_clust_observed.shape}"
        )
    elif prenormalize:
        # Ensure the initial cluster proportions sum to 1.
        guess_p_clust = _normalize(_clip(guess_p_clust))
    # Determine the initial guess for the read coordinate distribution.
    if guess_p_ends is None:
        # If no initial guess was given, then use the coordinates of
        # the observed reads.
        if dims[CLUSTERS] == 1:
            guess_p_ends = p_ends_observed[:, :, 0]
        else:
            guess_p_ends = np.average(p_ends_observed,
                                      axis=2,
                                      weights=guess_p_clust)
    elif guess_p_ends.shape != p_ends_observed.shape[:2]:
        raise ValueError(
            "If given, guess_p_ends must have dimensions (positions x "
            f"positions), but got {guess_p_ends.shape} "
            f"≠ {p_ends_observed.shape[:2]}"
        )
    elif prenormalize:
        # Ensure the initial coordinate distributions are normalized.
        guess_p_ends = _triu_norm(_clip(guess_p_ends))
    # Iteratively update the mutation rates and read coordinates.
    for i in range(max_iter):
        # Update the mutation rates using the read coordinates.
        next_p_mut_given_span = _calc_p_mut_given_span(p_mut_given_span_observed,
                                                       min_gap,
                                                       guess_p_ends,
                                                       guess_p_mut_given_span,
                                                       **kwargs)
        # Compute the RMSD change in mutation rates.
        rmsd_p_mut_given_span = np.sqrt(np.mean(np.square(
            next_p_mut_given_span - guess_p_mut_given_span
        )))
        # Update guess_p_mut_given_span for the next iteration.
        guess_p_mut_given_span = next_p_mut_given_span
        # Check for convergence using the RMSD change.
        if rmsd_p_mut_given_span <= convergence_thresh:
            break
        # Compute the probability that reads with each end coordinates
        # would have no two mutations too close.
        p_noclose_given_ends = _calc_p_noclose_given_ends(
            guess_p_mut_given_span,
            _calc_p_nomut_window(guess_p_mut_given_span, min_gap)
        )
        # Update the distribution of read end coordinates and clusters
        # using the mutation rates.
        guess_p_ends = _calc_p_ends(p_ends_observed,
                                    p_noclose_given_ends,
                                    guess_p_mut_given_span,
                                    guess_p_clust)
        guess_p_clust = calc_p_clust(p_clust_observed,
                                     calc_p_noclose(guess_p_ends,
                                                    p_noclose_given_ends))
    else:
        logger.warning("Mutation rates and distribution of read coordinates "
                       f"failed to converge in {max_iter} iterations.")
    return guess_p_mut_given_span, guess_p_ends, guess_p_clust


def calc_p_ends_given_noclose(p_ends: np.ndarray,
                              p_noclose_given_ends: np.ndarray):
    """ Calculate the proportion of reads with no two mutations too
    close with each pair of 5' and 3' coordinates.

    Assumptions
    -----------
    -   `p_ends` has 2 dimensions: (positions x clusters)
    -   Every value in the upper triangle of `p_ends` is ≥ 0 and ≤ 1;
        no values below the main diagonal are used.
    -   The upper triangle of `p_ends` sums to 1.
    -   `min_gap` is a non-negative integer.
    -   `p_mut_given_span` has 2 dimensions:
        (positions x clusters)
    -   Every value in `p_mut_given_span` is ≥ 0 and ≤ 1.
    -   There is at least 1 cluster.

    Parameters
    ----------
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads in
        each cluster beginning at the row position and ending at the
        column position.
    p_noclose_given_ends: np.ndarray
        3D (positions x positions x clusters) array of the pobabilities
        that a read with 5' and 3' coordinates corresponding to the row
        and column would have no two mutations too close.

    Returns
    -------
    np.ndarray
        3D (positions x positions x clusters) array of the proportion of
        reads without mutations too close, beginning at the row position
        and ending at the column position, in each cluster.
    """
    # Validate the dimensions of the arguments.
    find_dims([(POSITIONS, POSITIONS), (POSITIONS, POSITIONS, CLUSTERS)],
              [p_ends, p_noclose_given_ends],
              ["p_ends", "p_noclose_given_ends"],
              nonzero=True)
    # Calculate the proportion of total reads that would have each
    # pair of end coordinates.
    return _triu_norm(_triu_mul(p_ends[:, :, np.newaxis], p_noclose_given_ends))


@jit()
def _calc_p_ends_observed(npos: int,
                          end5s: np.ndarray,
                          end3s: np.ndarray,
                          weights: np.ndarray):
    """ Calculate the proportion of each pair of 5'/3' end coordinates
    in `end5s` and `end3s`, weighted by `weights`.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `npos` is a non-negative integer.
    -   `end5s` and `end3s` are each a 1D array of integers whose length
        equals the number of reads.
    -   `weights` is a 2D array of floats whose length and width are the
        numbers of reads and clusters, respectively.
    -   All integers in `end5s` and `end3s` are ≥ 0 and < `npos`.

    Parameters
    ----------
    npos: int
        Number of positions.
    end5s: np.ndarray
        5' end coordinates of the reads: 1D array (reads)
    end3s: np.ndarray
        3' end coordinates of the reads: 1D array (reads)
    weights: np.ndarray
        Number of times each read occurs in each cluster:
        2D array (reads x clusters)

    Returns
    -------
    np.ndarray
        Fraction of reads with each 5' (row) and 3' (column) coordinate:
        3D array (positions x positions x clusters)
    """
    nreads, ncls = weights.shape
    p_ends = np.zeros((npos, npos, ncls))
    for i in range(nreads):
        p_ends[end5s[i], end3s[i]] += weights[i]
    return p_ends


def calc_p_ends_observed(npos: int,
                         end5s: np.ndarray,
                         end3s: np.ndarray,
                         weights: np.ndarray | None = None,
                         check_values: bool = True):
    """ Calculate the proportion of each pair of 5'/3' end coordinates
    observed in `end5s` and `end3s`, optionally weighted by `weights`.

    Parameters
    ----------
    npos: int
        Number of positions.
    end5s: np.ndarray
        5' end coordinates of the reads: 1D array (reads)
    end3s: np.ndarray
        3' end coordinates of the reads: 1D array (reads)
    weights: np.ndarray | None = None
        Number of times each read occurs in each cluster:
        2D array (reads x clusters)
    check_values: bool = True
        Check that `end5s`, `end3s`, and `weights` are all valid.

    Returns
    -------
    np.ndarray
        Fraction of reads with each 5' (row) and 3' (column) coordinate:
        3D array (positions x positions x clusters)
    """
    if npos < 0:
        raise ValueError(f"Number of positions must be ≥ 0, but got {npos}")
    # Validate the dimensions.
    if weights is None:
        # Assume all reads are equally likely and there is one cluster.
        return calc_p_ends_observed(npos,
                                    end5s,
                                    end3s,
                                    weights=np.ones_like(end5s),
                                    check_values=check_values)
    if weights.ndim == 1:
        # There is one cluster: return a 2D array for that cluster.
        return calc_p_ends_observed(npos,
                                    end5s,
                                    end3s,
                                    weights=weights[:, np.newaxis],
                                    check_values=check_values)[:, :, 0]
    find_dims([(READS,), (READS,), (READS, CLUSTERS)],
              [end5s, end3s, weights],
              ["end5s", "end3s", "weights"])
    if check_values:
        # Validate the values.
        if end5s.min(initial=npos - 1) < 0:
            raise ValueError(
                f"All end5s must be ≥ 0, but got {end5s[end5s < 0]}"
            )
        if end5s.max(initial=0) >= npos:
            raise ValueError(
                f"All end5s must be < {npos}, but got {end5s[end5s >= npos]}"
            )
        if end3s.min(initial=npos - 1) < 0:
            raise ValueError(
                f"All end3s must be ≥ 0, but got {end3s[end3s < 0]}"
            )
        if end3s.max(initial=0) >= npos:
            raise ValueError(
                f"All end3s must be < {npos}, but got {end3s[end3s >= npos]}"
            )
        if weights.min(initial=1.) < 0.:
            raise ValueError(
                f"All weights must be ≥ 0, but got {weights[weights < 0.]}"
            )
    # Call the compiled function.
    return _calc_p_ends_observed(npos, end5s, end3s, weights)


def calc_params_observed(n_pos_total: int,
                         k: int,
                         unmasked_pos: Iterable[int],
                         muts_per_pos: Iterable[np.ndarray],
                         end5s: np.ndarray,
                         end3s: np.ndarray,
                         counts_per_uniq: np.ndarray,
                         resps: np.ndarray):
    """ Calculate the observed estimates of the parameters.

    Parameters
    ----------
    n_pos_total: int
        Total number of positions in the section.
    k: int
        Number of clusters.
    unmasked_pos: Iterable[int]
        Unmasked positions; must be zero-indexed with respect to the
        5' end of the section.
    muts_per_pos: Iterable[np.ndarray]
        For each unmasked position, numbers of all reads with a mutation
        at that position.
    end5s: np.ndarray
        Coordinates of the 5' ends of all reads; must be 0-indexed with
        respect to the 5' end of the section.
    end3s: np.ndarray
        Coordinates of the 3' ends of all reads; must be 0-indexed with
        respect to the 5' end of the section.
    counts_per_uniq: np.ndarray
        Number of times each unique read occurs.
    resps: np.ndarray
        Cluster memberships of each read: 2D array (reads x clusters)

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    # Count each unique read in each cluster.
    # 2D (unique reads x clusters)
    n_each_read_each_clust = counts_per_uniq[:, np.newaxis] * resps
    # Count the total number of reads in each cluster.
    # 1D (clusters)
    n_reads_per_clust = np.sum(n_each_read_each_clust, axis=0)
    # Calculate the observed proportion of reads in each cluster.
    # 1D (clusters)
    p_clust_observed = n_reads_per_clust / n_reads_per_clust.sum()
    # Calculate the proportion of each unique read in each cluster.
    # 2D (unique reads x clusters)
    p_each_read_each_clust = n_each_read_each_clust / n_reads_per_clust
    # Calculate the proportion of observed reads with each pair of
    # 5' and 3' end coordinates.
    # 3D (all positions x all positions x clusters)
    p_ends_observed = calc_p_ends_observed(n_pos_total,
                                           end5s,
                                           end3s,
                                           p_each_read_each_clust,
                                           check_values=False)
    # Count the observed reads that cover each position.
    # 2D (all positions x clusters)
    n_reads_per_pos = (calc_rectangluar_sum(p_ends_observed)
                       * n_reads_per_clust)
    # Count the observed mutations at each position.
    # 2D (all positions x clusters)
    n_muts_per_pos = np.zeros((n_pos_total, k))
    for j, mut_reads in zip(unmasked_pos, muts_per_pos, strict=True):
        # Calculate the number of mutations at each position in each
        # cluster by summing the count-weighted likelihood that each
        # read with a mutation at (j) came from the cluster.
        n_muts_per_pos[j] = counts_per_uniq[mut_reads] @ resps[mut_reads]
    # Calculate the observed mutation rate at each position.
    # 2D (all positions x clusters)
    with np.errstate(invalid="ignore"):
        # For any positions with 0 divided by 0, replace the quotient
        # with 0.
        p_mut_observed = np.nan_to_num(n_muts_per_pos / n_reads_per_pos)
    return p_mut_observed, p_ends_observed, p_clust_observed

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
