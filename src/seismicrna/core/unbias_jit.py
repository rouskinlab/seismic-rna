"""Numba-jitted helpers for core.unbias.

Separated from ``unbias`` so that importing ``unbias`` does not import
numba (slow). ``unbias``'s public wrappers import these lazily via a
proxy; see ``unbias._jit``. Importing this module imports numba.
"""

import warnings

import numpy as np
from numba import njit, NumbaPerformanceWarning

# Disable performance warnings from Numba.
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

# Maximum allowed mutation rate (a mutation rate of 1.0 will break this
# algorithm because it involves the log of 1 minus the mutation rate).
MAX_MU = 0.999999


@njit()
def clip_jit(x: np.ndarray):
    """Fill NaN with 0, infinity with 1, and restrict all values to the
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
    return np.clip(np.nan_to_num(x), 0.0, MAX_MU)


@njit()
def normalize_jit(x: np.ndarray):
    """Normalize the values to sum to 1, or if they sum to 0, then
    return an array with the same value for each element."""
    x_sum = np.sum(x)
    if x_sum == 0.0:
        if x.size == 0:
            # Handle the edge case of size zero, which would raise a
            # ZeroDivisionError if handled in the next branch.
            return np.ones_like(x)
        # If the sum of the input array is 0 and the array has at least
        # 1 element, then return an array of the same size as the input
        # where all elements are equal and sum to 1.
        return np.full_like(x, 1.0 / x.size)
    # Divide each element by the sum of all elements so the resulting
    # array has the same size as the input and all elements sum to 1.
    return x / x_sum


@njit()
def triu_sum_jit(a: np.ndarray):
    """Calculate the sum over the upper triangle(s) of array `a`.

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


@njit()
def triu_cumsum_jit(a: np.ndarray):
    """Calculate the cumulative sum from the upper right corner of `a`
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


@njit()
def triu_div_jit(numer: np.ndarray, denom: np.ndarray):
    """Divide the upper triangles of `numer` and `denom`.

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


@njit()
def triu_mul_jit(factor1: np.ndarray, factor2: np.ndarray):
    """Multiply the upper triangles of `numer` and `denom`.

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


@njit()
def triu_dot_jit(a: np.ndarray, b: np.ndarray):
    """Dot product of `a` and `b` over their first 2 dimensions.

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


@njit()
def adjust_min_gap_jit(num_pos: int, min_gap: int):
    """Given the number of positions (`npos`) and the desired minimum
    gap between mutations (`min_gap`), find the minimum gap between
    mutations that is no smaller than 0 and, if possible, no larger
    than 1 less than the number of positions.

    Parameters
    ----------
    num_pos: int
        Number of positions in the reference region.
    min_gap: int
        Minimum number of non-mutated bases between two mutations.

    Returns
    -------
    int
        Adjusted minimum gap between mutations.
    """
    return max(min(min_gap, num_pos - 1), 0)


@njit()
def calc_p_nomut_window_jit(p_mut_given_span: np.ndarray, min_gap: int):
    """Given underlying mutation rates (`p_mut_given_span`), find the
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
        3D (window x positions + 1 x clusters) array of the probability that
        (window) consecutive bases, ending at position (position), would
        have 0 mutations at all.
    """
    # Find and validate the dimensions.
    npos, ncls = p_mut_given_span.shape
    min_gap = adjust_min_gap_jit(npos, min_gap)
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
    p_nomut_window[0] = 1.0
    # For window size 1, every window represents exactly 1 position, so
    # the product over the window is just the non-mutation rate of the
    # corresponding position.
    p_nomut_window[1, 1:] = 1.0 - p_mut_given_span
    # For each window size from 2 to min_gap, calculate the probability
    # of there being no mutations in each window of that size.
    for size in range(2, p_nomut_window.shape[0]):
        shift = size - 1
        p_nomut_window[size, size:] = (
            p_nomut_window[1, size:] * p_nomut_window[shift, shift:-1]
        )
    return p_nomut_window


@njit()
def calc_p_noclose_given_ends_jit(
    p_mut_given_span: np.ndarray, p_nomut_window: np.ndarray
):
    """Given underlying mutation rates (`p_mut_given_span`), calculate
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
        p_noclose_given_ends[j, max(j - min_gap, 0) : (j + 1)] = 1.0
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
            p_nomut_window[window_indexes[(npos - j) :], j]
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


@njit()
def calc_rectangular_sum_jit(array: np.ndarray):
    """For each element of the main diagonal, calculate the sum over
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
        cols_area[0] = 0.0
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


@njit()
def calc_rectangular_sum_weighted_jit(array: np.ndarray, weights: np.ndarray):
    """Like `_calc_rectangular_sum`, but multiplies each element of
    `array` by the corresponding scalar in `weights` on the fly, avoiding
    the allocation of a full product array.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `array` has exactly 3 dimensions (positions x positions x clusters).
    -   `weights` has exactly 2 dimensions (positions x positions).
    -   The first and second dimensions of both arrays have equal lengths.
    """
    npos = array.shape[0]
    ncls = array.shape[2]
    rows_area = np.empty((npos, ncls))
    cols_area = np.empty((npos, ncls))
    if npos > 0:
        rows_area[0] = weights[0] @ array[0]
        cols_area[0] = 0.0
        for j in range(1, npos):
            i = j - 1
            rows_area[j] = rows_area[i] + weights[j, j:] @ array[j, j:]
            # Use an element loop for the column sum to avoid non-contiguous
            # '@': array[a, i] is a contiguous 1D row of the last dimension,
            # whereas array[:j, i] as a whole is a non-contiguous 2D slice.
            col_contrib = np.zeros(ncls)
            for a in range(j):
                col_contrib += weights[a, i] * array[a, i]
            cols_area[j] = cols_area[i] + col_contrib
    return rows_area - cols_area


@njit()
def calc_p_mut_given_span_dropped_jit(
    p_mut_given_span: np.ndarray,
    p_ends: np.ndarray,
    p_noclose_given_ends: np.ndarray,
    p_nomut_window: np.ndarray,
):
    """Calculate the mutation rates after dropping reads with two
    mutations that are too close.

    Parameters
    ----------
    p_mut_given_span: np.ndarray
        2D (positions x clusters) array of the underlying mutation rates
        (i.e. the probability that a read has a mutation at position (j)
        given that it contains that position).
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads in
        each cluster beginning at the row position and ending at the
        column position.
    p_noclose_given_ends: np.ndarray
        3D (positions x positions x clusters) array of the probabilities
        that a read with 5' and 3' coordinates corresponding to the row
        and column would have no two mutations too close.
    p_nomut_window: np.ndarray
        3D (window x positions x clusters) array of the probability that
        (window) consecutive bases, ending at position (position), would
        have zero mutations at all.

    Returns
    -------
    np.ndarray
        2D (positions x clusters) array of the mutation rate among reads
        with no two mutations too close per position per cluster.
    """
    # Find the dimensions.
    npos, ncls = p_mut_given_span.shape
    min_gap = p_nomut_window.shape[0] - 1
    if min_gap <= 0:
        # If min_gap <= 0, then no mutations can be too close, so return
        # the original mutation rates.
        return p_mut_given_span
    # Calculate the probability that a read spanning each position would
    # have no two mutations too close.
    p_noclose_given_span = calc_rectangular_sum_weighted_jit(
        p_noclose_given_ends, p_ends
    )
    p_mut_given_span_noclose = np.nan_to_num(p_mut_given_span / p_noclose_given_span)
    # Preallocate work buffers (clusters x positions) at maximum size to
    # avoid per-iteration heap allocation. Row-leading layout [(cls, pos)]
    # means buf[k, :n] is always a 1D contiguous row prefix, which keeps
    # all '@' operands typed as 1D 'C' by Numba (no non-contiguous slices).
    p_noclose_a = np.empty((ncls, npos))
    p_noclose_b = np.empty((ncls, npos))
    right_mat = np.empty((ncls, npos))
    # Compute the mutation rates given no two mutations are too close
    # one position (j) at a time.
    for j in range(npos):
        nrows = j + 1
        ncols = npos - j
        # For each starting position (a), calculate the probability of
        # no two mutations too close from (a) to (j).
        for a in range(nrows):
            # Number of bases 5' of (j) that must have no mutations:
            # the smaller of (min_gap) and (j - a).
            nomut = min(min_gap, j - a)
            # Probability of no close mutations, (a) to (j - 1 - nomut):
            #   p_noclose_given_ends[a, max(j - 1 - nomut, a)]
            # Probability of no mutations over the (nomut)-sized window
            # from (j - nomut) to (j - 1):
            #   p_nomut_window[nomut, j]
            # Probability of no close mutations from (a) to (j).
            p_noclose_a[:, a] = (
                p_noclose_given_ends[a, max(j - 1 - nomut, a)]
                * p_nomut_window[nomut, j]
            )
        # For each ending position (b), calculate the probability of
        # no two mutations too close from (j) to (b).
        for b in range(j, npos):
            # Number of bases 3' of (j) that must have no mutations:
            # the smaller of (min_gap) and (b - j).
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
        # Build right_mat[k, a] = p_ends[a, j:] · p_noclose_b[k, :ncols]
        # using only 1D row-prefix slices so every '@' operand is typed
        # as 1D 'C' by Numba (a 2D slice would be typed 'A' even when
        # contiguous at runtime, triggering a performance warning).
        for k in range(ncls):
            for a in range(nrows):
                right_mat[k, a] = p_ends[a, j:] @ p_noclose_b[k, :ncols]
            p_mut_given_span_noclose[j, k] *= (
                p_noclose_a[k, :nrows] @ right_mat[k, :nrows]
            )
    return p_mut_given_span_noclose


@njit()
def calc_p_mut_given_span_merged_jit(
    p_mut_given_span: np.ndarray, p_ends: np.ndarray, min_gap: int
) -> np.ndarray:
    """Calculate the mutation rates after merging mutations closer than
    min_gap positions.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    - `p_mut_given_span` is a 2D (positions x clusters) array of the
      underlying mutation rates.
    - All values of `p_mut_given_span` are between 0 and 1, inclusive.
    - `p_ends` is a 2D (positions x positions) array of the proportion
      of reads in each cluster beginning at the row position and ending
      at the column position. (Only the upper triangle of `p_ends` is
      used, so values below the main diagonal can be ignored.)
    - All values in the upper triangle of `p_ends` are ≥ 0.
    - The sum of the upper triangle of `p_ends` equals 1.
    - `min_gap` is a non-negative integer.

    Parameters
    ----------
    p_mut_given_span: np.ndarray
        2D (positions x clusters) array of the underlying mutation rates
        (i.e. the probability that a read has a mutation at position (j)
        given that it contains that position).
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads in
        each cluster beginning at the row position and ending at the
        column position.
    min_gap: int
        Minimum number of non-mutated bases between two mutations;
        must be ≥ 0.

    Returns
    -------
    np.ndarray
        2D (positions x clusters) array of the mutation rate among reads
        with no two mutations too close per position per cluster.
    """
    if min_gap <= 0:
        # No positions will be merged.
        return p_mut_given_span.copy()
    npos = p_mut_given_span.shape[0]
    ncls = p_mut_given_span.shape[1]
    # Initialize the probability that a read spans each position.
    p_span = np.zeros(npos, dtype=p_mut_given_span.dtype)
    # Initialize the joint probability that reads both span and have a
    # mutation at each position after merging.
    p_joint_mut_span_merged = np.zeros_like(p_mut_given_span)
    # Preallocate a single buffer for the merged rates; reused every b.
    p_merged_buf = np.empty_like(p_mut_given_span)
    # Preallocate a running window sum over clusters.
    window_sum = np.empty(ncls)
    # For a given end position b, merged probabilities do not depend on
    # the start position a, so compute them once for each end position b
    # from 0 to (npos - 1).
    for b in range(npos):
        b_inc = b + 1
        # Calculate the joint probability that a read ends at position b
        # and spans each position from 0 to b (inclusive).
        p_joint_span_end_b = np.cumsum(p_ends[:b_inc, b])
        # If no reads end at position b (p_joint_span_end_b[b] == 0),
        # then there is nothing to do, so skip to the next b.
        if p_joint_span_end_b[b] > 0.0:
            # Accumulate the probability that a read spans each position
            # by adding the joint probabilities for each end position b.
            p_span[:b_inc] += p_joint_span_end_b
            # Calculate the merged rates for reads ending at position b
            # using an inlined right-to-left recurrence.
            # The 3'-most position b has no neighbours on its right
            # within [0, b], so its merged rate equals the raw rate.
            p_merged_buf[b, :] = p_mut_given_span[b, :]
            # Initialize the running window sum for j = b-1.
            # The window at position j covers p_merged_buf[j+1:j+min_gap+1]
            # clipped to [0, b]; for j = b-1 that is just {p_merged_buf[b]}.
            window_sum[:] = p_merged_buf[b, :]
            # Fill positions b-1 down to 0 using the running window sum,
            # keeping update cost O(ncls) per step instead of O(min_gap * ncls).
            for j in range(b - 1, -1, -1):
                p_merged_buf[j, :] = p_mut_given_span[j, :] * (1.0 - window_sum[:])
                # Slide the window: add p_merged_buf[j], drop position
                # j + min_gap if it is still within [0, b].
                window_sum[:] += p_merged_buf[j, :]
                out_pos = j + min_gap
                if out_pos <= b:
                    window_sum[:] -= p_merged_buf[out_pos, :]
            # Accumulate the joint probability that a read has a mutation
            # at each position after merging and ends at position b.
            for j in range(b_inc):
                wt = p_joint_span_end_b[j]
                p_joint_mut_span_merged[j, :] += p_merged_buf[j, :] * wt
    # Calculate the mutation rates after merging by dividing the joint
    # probabilities of mutation and spanning by the probabilities of
    # spanning. If a position is not spanned by any reads, the mutation
    # rate is set to zero.
    return np.nan_to_num(p_joint_mut_span_merged / p_span[:, np.newaxis])


@njit()
def calc_p_ends_observed_jit(
    npos: int, end5s: np.ndarray, end3s: np.ndarray, weights: np.ndarray
):
    """Calculate the proportion of each pair of 5'/3' end coordinates
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
        5' ends (0-indexed) of the reads: 1D array (reads)
    end3s: np.ndarray
        3' ends (0-indexed) of the reads: 1D array (reads)
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
