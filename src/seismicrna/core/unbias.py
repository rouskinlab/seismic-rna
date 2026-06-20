from __future__ import annotations
from typing import Iterable


from .arg.cli import MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE
from .array import find_dims, triangular
from .logs import logger
from .validate import (
    require_isinstance,
    require_equal,
    require_atleast,
    require_atmost,
    require_between,
)


# Define dimensions
READS = "reads"
UNIQUE_READS = "unique reads"
POSITIONS = "positions"
CLUSTERS = "clusters"
POSITIONS_PLUS_1 = "positions + 1"
WINDOW = "window"
SIZE = "size"


def require_square_atleast2d(name: str, array: np.ndarray):
    """Require the input to be a NumPy NDArray with ≥ 2 dimensions,
    and the first and second dimensions to be of equal length."""
    import numpy as np

    require_isinstance(name, array, np.ndarray)
    require_atleast(f"{name}.ndim", array.ndim, 2)
    require_equal(
        f"{name}.shape[0]", array.shape[0], array.shape[1], f"{name}.shape[1]"
    )


def require_same_square_atleast2d(a: np.ndarray, b: np.ndarray):
    """Require a and b to each be a NumPy NDArray with ≥ 2 dimensions,
    with the first and second dimensions of a and b all equal."""
    require_square_atleast2d("a", a)
    require_square_atleast2d("b", b)
    require_equal("a.shape[:2]", a.shape[:2], b.shape[:2], "b.shape[:2]")


def triu_norm(a: np.ndarray):
    """Normalize the upper triangle of array `a` to sum to 1.

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
        Array of the same shape as `a` but scaled so that the upper
        triangle sums to 1.
    """
    import numpy as np
    from .unbias_jit import triu_sum_jit, triu_div_jit
    # Calculate the sum over axes 0 and 1.
    a_sum = triu_sum_jit(a)
    # Normalize by dividing by that sum, ignoring division by zero,
    # which is handled subsequently.
    a_norm = triu_div_jit(a, np.broadcast_to(a_sum, a.shape))
    # Determine which sums to normalize by are zero.
    a_sum_zero = a_sum == 0.0
    if np.any(np.atleast_1d(a_sum_zero)):
        # Count the elements over which normalization occurred.
        n_norm = triangular(a.shape[0])
        if n_norm == 0:
            # Handle the edge case of size zero, which would raise a
            # ZeroDivisionError if handled in the next branch.
            return np.ones_like(a_norm)
        # For each sum that was zero, set all elements to the same value
        # such that they sum to 1.
        fill = 1.0 / n_norm
        if a_sum_zero.ndim > 0:
            for indexes in zip(*np.nonzero(a_sum_zero), strict=True):
                a_norm[(slice(None),) * 2 + indexes] = fill
        else:
            a_norm[:, :] = fill
    return a_norm


def triu_dot(a: np.ndarray, b: np.ndarray):
    """Dot product of `a` and `b` over their first 2 dimensions.

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
    from .unbias_jit import triu_dot_jit
    require_same_square_atleast2d(a, b)
    return triu_dot_jit(a, b)


def calc_p_nomut_window(p_mut_given_span: np.ndarray, min_gap: int):
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
        3D (window x positions + 1 x clusters) array of the probability
        that (window) consecutive bases, ending at position (position),
        would have 0 mutations at all.
    """
    import numpy as np
    from .unbias_jit import calc_p_nomut_window_jit
    require_isinstance("p_mut_given_span", p_mut_given_span, np.ndarray)
    require_equal("p_mut_given_span.ndim", p_mut_given_span.ndim, 2)
    require_atleast("min_gap", min_gap, 0, classes=int)
    return calc_p_nomut_window_jit(p_mut_given_span, min_gap)


def calc_p_noclose_given_ends(p_mut_given_span: np.ndarray, p_nomut_window: np.ndarray):
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
    from .unbias_jit import calc_p_noclose_given_ends_jit
    dims = find_dims(
        [(POSITIONS, CLUSTERS), (WINDOW, POSITIONS_PLUS_1, CLUSTERS)],
        [p_mut_given_span, p_nomut_window],
        ["p_mut_given_span", "p_nomut_window"],
    )
    require_equal(
        "p_mut_given_span.shape[0]",
        dims[POSITIONS],
        dims[POSITIONS_PLUS_1] - 1,
        "p_nomut_window.shape[1] - 1",
    )
    return calc_p_noclose_given_ends_jit(p_mut_given_span, p_nomut_window)


def calc_p_noclose_given_ends_auto(p_mut_given_span: np.ndarray, min_gap: int):
    """Given underlying mutation rates (`p_mut_given_span`), calculate
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
    import numpy as np
    from .unbias_jit import (
        clip_jit,
        calc_p_noclose_given_ends_jit,
        calc_p_nomut_window_jit,
    )
    require_isinstance("p_mut_given_span", p_mut_given_span, np.ndarray)
    require_atleast("min_gap", min_gap, 0, classes=int)
    if p_mut_given_span.ndim == 2:
        p_mut_given_span = clip_jit(p_mut_given_span)
        return calc_p_noclose_given_ends_jit(
            p_mut_given_span, calc_p_nomut_window_jit(p_mut_given_span, min_gap)
        )
    if p_mut_given_span.ndim == 1:
        return calc_p_noclose_given_ends_auto(
            p_mut_given_span[:, np.newaxis], min_gap
        ).reshape((p_mut_given_span.size, p_mut_given_span.size))
    raise ValueError(
        f"p_mut_given_span.ndim must equal 1 or 2, but got {p_mut_given_span.ndim}"
    )


def calc_rectangular_sum(array: np.ndarray):
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
    from .unbias_jit import calc_rectangular_sum_jit
    require_square_atleast2d("array", array)
    return calc_rectangular_sum_jit(array)


def calc_p_mut_given_span_dropped(
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
    from .unbias_jit import calc_p_mut_given_span_dropped_jit
    dims = find_dims(
        [
            (POSITIONS, CLUSTERS),
            (POSITIONS, POSITIONS),
            (POSITIONS, POSITIONS, CLUSTERS),
            (WINDOW, POSITIONS_PLUS_1, CLUSTERS),
        ],
        [p_mut_given_span, p_ends, p_noclose_given_ends, p_nomut_window],
        ["p_mut_given_span", "p_ends", "p_noclose_given_ends", "p_nomut_window"],
    )
    require_equal(
        "p_mut_given_span.shape[0]",
        dims[POSITIONS],
        dims[POSITIONS_PLUS_1] - 1,
        "p_nomut_window.shape[1] - 1",
    )
    return calc_p_mut_given_span_dropped_jit(
        p_mut_given_span, p_ends, p_noclose_given_ends, p_nomut_window
    )


def _calc_p_mut_given_span_biased(
    p_mut_given_span: np.ndarray, p_ends: np.ndarray, min_gap: int, mut_collisions: str
):
    """Calculate the biased mutation rates after dropping reads with
    two mutations too close or merging mutations too close.

    Parameters
    ----------
    p_mut_given_span: np.ndarray
        2D (positions x clusters) array of the underlying mutation rates.
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads
        beginning at the row position and ending at the column position.
    min_gap: int
        Minimum number of non-mutated bases between two mutations.
    mut_collisions: str
        Method for handling reads with two mutations that are too close;
        must be either "drop" or "merge".

    Returns
    -------
    np.ndarray
        2D (positions x clusters) array of the biased mutation rates.
    """
    from .unbias_jit import (
        calc_p_nomut_window_jit,
        calc_p_noclose_given_ends_jit,
        calc_p_mut_given_span_dropped_jit,
        calc_p_mut_given_span_merged_jit,
    )
    if mut_collisions == MUT_COLLISIONS_DROP:
        # Use the read-dropping method.
        p_nomut_window = calc_p_nomut_window_jit(p_mut_given_span, min_gap)
        p_noclose_given_ends = calc_p_noclose_given_ends_jit(
            p_mut_given_span, p_nomut_window
        )
        return calc_p_mut_given_span_dropped_jit(
            p_mut_given_span, p_ends, p_noclose_given_ends, p_nomut_window
        )
    if mut_collisions == MUT_COLLISIONS_MERGE:
        # Use the mutation-merging method.
        return calc_p_mut_given_span_merged_jit(p_mut_given_span, p_ends, min_gap)
    raise ValueError(
        f"mut_collisions must be either {repr(MUT_COLLISIONS_DROP)} "
        f"or {repr(MUT_COLLISIONS_MERGE)}, but got {repr(mut_collisions)}"
    )


def slice_p_ends(p_ends: np.ndarray, p_ends_cumsum: np.ndarray, end5: int, end3: int):
    """Slice a matrix of end coordinate probabilities to the subregion
    from `end5` to `end3`, redistributing out-of-region probability mass
    onto the boundary cells.

    Parameters
    ----------
    p_ends: np.ndarray
        2D (positions x positions) array of end coordinate proportions.
    p_ends_cumsum: np.ndarray
        Cumulative sum of the upper triangle of `p_ends` (from the
        upper-right corner), as computed by `_jit._triu_cumsum`.
    end5: int
        Index of the 5' boundary of the slice (inclusive).
    end3: int
        Index one past the 3' boundary of the slice (exclusive).

    Returns
    -------
    np.ndarray
        Sliced (end3 - end5) x (end3 - end5) array of end coordinate
        proportions.
    """
    p_ends_slice = p_ends[end5:end3, end5:end3].copy()
    if p_ends_slice.size > 0:
        p_ends_slice[0, :-1] = p_ends[: end5 + 1, end5 : end3 - 1].sum(axis=0)
        p_ends_slice[1:, -1] = p_ends[end5 + 1 : end3, end3 - 1 :].sum(axis=1)
        p_ends_slice[0, -1] = p_ends_cumsum[end5, end3 - 1]
    return p_ends_slice


def find_split_positions(p_mut: np.ndarray, min_gap: int, threshold: float):
    """Find positions at which to split the mutation rates into
    independent segments, defined as stretches of at least `min_gap`
    consecutive positions all at or below `threshold`.

    Parameters
    ----------
    p_mut: np.ndarray
        2D (positions x clusters) array of mutation rates.
    min_gap: int
        Minimum length of a low-mutation stretch to use as a split.
    threshold: float
        Mutation rate at or below which a position is considered
        low-mutation.

    Returns
    -------
    np.ndarray
        1D array of split positions (integers).
    """
    import numpy as np
    from .unbias_jit import adjust_min_gap_jit
    npos, ncls = p_mut.shape
    min_gap = adjust_min_gap_jit(npos, min_gap)
    if min_gap == 0 or ncls == 0:
        return np.array([], dtype=int)
    # Count the cumulative number of mutation rates below the threshold.
    cum_below = np.cumsum(
        np.concatenate([np.zeros(min_gap, dtype=bool), p_mut.max(axis=1) <= threshold])
    )
    # Label every position for which it and the previous (min_gap - 1)
    # positions are all below the threshold.
    win_below = cum_below[min_gap:] - cum_below[:-min_gap] == min_gap
    # Locate each stretch of consecutive positions where the previous
    # (min_gap - 1) positions are all below the threshold; the first
    # (min_gap - 1) positions can be ignored because they cannot have
    # (min_gap - 1) preceding positions.
    diff_win_below = np.diff(win_below.astype(int)[(min_gap - 1) :])
    first_pos = np.flatnonzero(diff_win_below == 1) + 1
    last_pos = np.flatnonzero(diff_win_below == -1) + min_gap
    # Collect all boundary positions and find the unique ones.
    return np.unique(np.concatenate([first_pos, last_pos]))


def _split_positions(
    p_mut: np.ndarray, p_mut_init: np.ndarray, p_ends: np.ndarray, split_pos: np.ndarray
):
    """Split `p_mut`, `p_mut_init`, and `p_ends` at the given split
    positions.

    Parameters
    ----------
    p_mut: np.ndarray
        2D (positions x clusters) array of observed mutation rates.
    p_mut_init: np.ndarray
        2D (positions x clusters) array of initial mutation rate
        estimates.
    p_ends: np.ndarray
        2D (positions x positions) array of end coordinate proportions.
    split_pos: np.ndarray
        1D array of positions at which to split.

    Returns
    -------
    tuple[list, list, list]
        Lists of split subarrays for p_mut, p_mut_init, and p_ends.
    """
    import numpy as np
    from .unbias_jit import triu_cumsum_jit
    dims = find_dims(
        [(POSITIONS, CLUSTERS), (POSITIONS, CLUSTERS), (POSITIONS, POSITIONS)],
        [p_mut, p_mut_init, p_ends],
        ["p_mut", "p_mut_init", "p_ends"],
    )
    (n_splits,) = split_pos.shape
    if n_splits == 0:
        return [p_mut], [p_mut_init], [p_ends]
    npos = dims[POSITIONS]
    p_mut_split = np.split(p_mut, split_pos)
    p_mut_init_split = np.split(p_mut_init, split_pos)
    p_ends_cumsum = triu_cumsum_jit(p_ends)
    p_ends_split = [slice_p_ends(p_ends, p_ends_cumsum, 0, split_pos[0])]
    for i in range(n_splits - 1):
        p_ends_split.append(
            slice_p_ends(p_ends, p_ends_cumsum, split_pos[i], split_pos[i + 1])
        )
    p_ends_split.append(slice_p_ends(p_ends, p_ends_cumsum, split_pos[-1], npos))
    return p_mut_split, p_mut_init_split, p_ends_split


def _calc_p_mut_given_span(
    p_mut_given_span_observed: np.ndarray,
    p_ends: np.ndarray,
    min_gap: int,
    mut_collisions: str,
    init_p_mut_given_span: np.ndarray,
    f_tol: float,
    x_rtol: float,
):
    """Calculate the underlying mutation rates including for reads with
    two mutations too close based on the observed mutation rates.
    Do not validate the argument types or values, since it is assumed
    that calc_p_mut_given_span has already done so."""
    import numpy as np
    from .unbias_jit import clip_jit
    # Use the Newton-Krylov method to solve for the total mutation rates
    # (including reads with two mutations too close) that result in zero
    # difference between theoretically observed and actually observed
    # mutation rates.

    def objective(p_mut_given_span: np.ndarray):
        # Clip the iterate to [0, 1] before evaluating the biased model.
        # The merge formula is numerically unstable outside [0, 1]:
        # errors grow multiplicatively as the recursion proceeds from
        # the 3' end to the 5' end, which can cause Newton-Krylov to
        # spend a huge number of iterations chasing a runaway residual.
        return p_mut_given_span_observed - _calc_p_mut_given_span_biased(
            clip_jit(p_mut_given_span), p_ends, min_gap, mut_collisions
        )

    # Import scipy here instead of at the top of this module because
    # its import is slow enough to impact global startup time.
    from scipy.optimize import newton_krylov, NoConvergence

    try:
        return clip_jit(
            newton_krylov(
                objective,
                init_p_mut_given_span,
                f_tol=f_tol,
                x_rtol=x_rtol,
                maxiter=1000,
            )
        )
    except (NoConvergence, ValueError, ZeroDivisionError, OverflowError) as error:
        # The Newton-Krylov solver could not unbias the mutation rates.
        # Return the original mutation rates as a fallback.
        logger.warning(error)
        return init_p_mut_given_span


def calc_p_mut_given_span(
    p_mut_given_span_observed: np.ndarray,
    p_ends: np.ndarray,
    min_gap: int,
    mut_collisions: str,
    init_p_mut_given_span: np.ndarray,
    *,
    quick_unbias: bool = True,
    quick_unbias_thresh: float = 0.0,
    f_tol: float = 1.0e-4,
    x_rtol: float = 1.0e-3,
):
    """Calculate the underlying mutation rates including for reads with
    two mutations too close based on the observed mutation rates.

    Parameters
    ----------
    p_mut_given_span_observed: np.ndarray
        Observed mutation rates (among reads with no two mutations too
        close): 2D array (positions x clusters).
    p_ends: np.ndarray
        Proportion of reads with each pair of 5'/3' end coordinates:
        2D array (positions x positions).
    min_gap: int
        Minimum number of non-mutated bases between two mutations;
        must be ≥ 0.
    mut_collisions: str
        Method for handling reads with two mutations that are too close;
        must be either "drop" or "merge".
    init_p_mut_given_span: np.ndarray
        Initial guess for the underlying mutation rates:
        2D array (positions x clusters).
    quick_unbias: bool = True
        Use the approximation method to reduce time complexity from
        O(n^3) to O(n) by splitting the region at low-mutation sites.
    quick_unbias_thresh: float = 0.
        Threshold mutation rate below which positions are considered
        low-mutation for the purpose of splitting.
    f_tol: float = 1.e-4
        Tolerance on the residual for the Newton-Krylov solver.
    x_rtol: float = 1.e-3
        Relative tolerance on the step size for the Newton-Krylov solver.

    Returns
    -------
    np.ndarray
        Underlying (unbiased) mutation rates:
        2D array (positions x clusters).
    """
    import numpy as np
    from .unbias_jit import adjust_min_gap_jit
    # Validate the argument types, values, and dimensions (of arrays).
    require_isinstance(
        "p_mut_given_span_observed", p_mut_given_span_observed, np.ndarray
    )
    if p_mut_given_span_observed.ndim == 1:
        # If p_mut_given_span_noclose has 1 dimension, then convert it
        # to 2 dimensions, calculate, and convert back to 1 dimension.
        return calc_p_mut_given_span(
            p_mut_given_span_observed[:, np.newaxis],
            p_ends,
            min_gap,
            mut_collisions,
            init_p_mut_given_span,
            f_tol=f_tol,
            x_rtol=x_rtol,
        )[:, 0]
    dims = find_dims(
        [(POSITIONS, CLUSTERS), (POSITIONS, CLUSTERS), (POSITIONS, POSITIONS)],
        [p_mut_given_span_observed, init_p_mut_given_span, p_ends],
        ["p_mut_given_span_observed", "init_p_mut_given_span", "p_ends"],
        nonzero=True,
    )
    require_atleast("min_gap", min_gap, 0, classes=int)
    min_gap = adjust_min_gap_jit(dims[POSITIONS], min_gap)
    if min_gap == 0:
        # No two mutations can be too close.
        return p_mut_given_span_observed
    require_atleast("f_tol", f_tol, 0.0, classes=float)
    require_atleast("x_rtol", x_rtol, 0.0, classes=float)
    # Decide whether to use the approximation method.
    if quick_unbias:
        # Use the approximation method to reduce the time complexity
        # from O(n^3) to O(n) where n is the number of positions.
        require_between(
            "quick_unbias_thresh", quick_unbias_thresh, 0.0, 1.0, classes=float
        )
        # Split the mutation rates and end coordinate distributions,
        # calculate them separately, and reassemble them.
        p_mut_given_span_list = list()
        for p_mut_split, p_mut_init_split, p_ends_split in zip(
            *_split_positions(
                p_mut_given_span_observed,
                init_p_mut_given_span,
                p_ends,
                find_split_positions(
                    p_mut_given_span_observed, min_gap, quick_unbias_thresh
                ),
            ),
            strict=True,
        ):
            if p_mut_split.max(initial=0.0) > quick_unbias_thresh:
                p_mut_given_span_list.append(
                    _calc_p_mut_given_span(
                        p_mut_split,
                        p_ends_split,
                        min_gap,
                        mut_collisions,
                        p_mut_init_split,
                        f_tol,
                        x_rtol,
                    )
                )
            else:
                # No mutation rate exceeds the threshold.
                p_mut_given_span_list.append(p_mut_split)
        return (
            np.concatenate(p_mut_given_span_list, axis=0)
            if p_mut_given_span_list
            else np.empty((dims[POSITIONS], dims[CLUSTERS]))
        )
    else:
        # Do not use the approximation method.
        return _calc_p_mut_given_span(
            p_mut_given_span_observed,
            p_ends,
            min_gap,
            mut_collisions,
            init_p_mut_given_span,
            f_tol,
            x_rtol,
        )


def calc_p_ends(
    p_ends_observed: np.ndarray,
    p_noclose_given_ends: np.ndarray,
    p_mut_given_span: np.ndarray,
    p_clust: np.ndarray,
):
    """Calculate the proportion of total reads with each pair of 5' and
    3' coordinates.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   Every value in the upper triangle of `p_ends_observed` is
        ≥ 0 and ≤ 1; no values below the main diagonal are used.
    -   The upper triangle of `p_ends_observed` sums to 1.
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
    import numpy as np
    from .unbias_jit import triu_div_jit
    # Validate the dimensionality of the arguments.
    require_isinstance("p_ends_observed", p_ends_observed, np.ndarray)
    if p_ends_observed.ndim == 2:
        # If p_ends_observed has 2 dimensions, then promote it to 3
        # dimensions first.
        return calc_p_ends(
            p_ends_observed[:, :, np.newaxis],
            p_noclose_given_ends,
            p_mut_given_span,
            p_clust,
        )
    dims = find_dims(
        [
            (POSITIONS, POSITIONS, CLUSTERS),
            (POSITIONS, POSITIONS, CLUSTERS),
            (POSITIONS, CLUSTERS),
            (CLUSTERS,),
        ],
        [p_ends_observed, p_noclose_given_ends, p_mut_given_span, p_clust],
        ["p_ends_observed", "p_noclose_given_ends", "p_mut_given_span", "p_clust"],
        nonzero=True,
    )
    # Calculate the proportion of total reads that would have each
    # pair of end coordinates.
    p_ends = triu_norm(triu_div_jit(p_ends_observed, p_noclose_given_ends))
    # Return a consensus distribution among all clusters.
    if dims[CLUSTERS] == 1:
        return p_ends[:, :, 0]
    return np.average(p_ends, axis=2, weights=p_clust)


def calc_p_noclose_given_clust(p_ends: np.ndarray, p_noclose_given_ends: np.ndarray):
    """Calculate the probability that a read from each cluster would
    have no two mutations too close.

    Parameters
    ----------
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads
        beginning at the row position and ending at the column position.
    p_noclose_given_ends: np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a read with 5' and 3' coordinates corresponding to the row
        and column would have no two mutations too close.

    Returns
    -------
    np.ndarray
        1D (clusters) array of the probability that a read from each
        cluster would have no two mutations too close.
    """
    import numpy as np
    from .unbias_jit import triu_dot_jit
    # Validate the dimensionality of the arguments.
    find_dims(
        [(POSITIONS, POSITIONS), (POSITIONS, POSITIONS, CLUSTERS)],
        [p_ends, p_noclose_given_ends],
        ["p_ends", "p_noclose_given_ends"],
        nonzero=True,
    )
    # Compute the weighted sum of the probabilities that reads from each
    # cluster would have no two mutations too close.
    return triu_dot_jit(p_ends[:, :, np.newaxis], p_noclose_given_ends)


def calc_p_clust(p_clust_observed: np.ndarray, p_noclose_given_clust: np.ndarray):
    """Cluster proportion among all reads.

    Parameters
    ----------
    p_clust_observed: np.ndarray
        Proportion of each cluster among reads with no two mutations too
        close.
        1D (clusters)
    p_noclose_given_clust: np.ndarray
        Probability that a read from each cluster would have no two
        mutations too close.
        1D (clusters)

    Returns
    -------
    np.ndarray
        Proportion of each cluster among all reads.
        1D (clusters)
    """
    from .unbias_jit import normalize_jit
    # Validate the dimensions.
    find_dims(
        [(CLUSTERS,), (CLUSTERS,)],
        [p_clust_observed, p_noclose_given_clust],
        ["p_clust_observed", "p_noclose_given_clust"],
        nonzero=True,
    )
    # The cluster proportions among all reads are obtained by weighting
    # each cluster proportion among reads with no two mutations too
    # close by the reciprocal of the probability that no two mutations
    # are too close in that cluster, then normalizing so the sum is 1.
    return normalize_jit(p_clust_observed / p_noclose_given_clust)


def calc_p_clust_given_noclose(p_clust: np.ndarray, p_noclose_given_clust: np.ndarray):
    """Cluster proportions among reads with no two mutations too close.

    Parameters
    ----------
    p_clust: np.ndarray
        Proportion of each cluster among all reads.
        1D (clusters)
    p_noclose_given_clust: np.ndarray
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
    from .unbias_jit import normalize_jit
    # Validate the dimensions.
    find_dims(
        [(CLUSTERS,), (CLUSTERS,)],
        [p_clust, p_noclose_given_clust],
        ["p_clust", "p_noclose_given_clust"],
        nonzero=True,
    )
    # The cluster proportions among reads with no two mutations too
    # close are obtained by weighting each cluster proportion by the
    # probability that no two mutations are too close in that cluster,
    # then normalizing so the sum is 1.
    return normalize_jit(p_clust * p_noclose_given_clust)


def calc_p_noclose(p_clust: np.ndarray, p_noclose_given_clust: np.ndarray):
    """Probability that any read would have two mutations too close.

    Parameters
    ----------
    p_clust: np.ndarray
        Proportion of each cluster among all reads.
        1D (clusters)
    p_noclose_given_clust: np.ndarray
        Probability that a read from each cluster would have no two
        mutations too close.
        1D (clusters)

    Returns
    -------
    float
        Probability that any read would have no two mutations too close.
    """
    import numpy as np
    find_dims(
        [(CLUSTERS,), (CLUSTERS,)],
        [p_clust, p_noclose_given_clust],
        ["p_clust", "p_noclose_given_clust"],
    )
    return float(np.vdot(p_clust, p_noclose_given_clust))


def calc_params(
    p_mut_given_span_observed: np.ndarray,
    p_ends_observed: np.ndarray,
    p_clust_observed: np.ndarray,
    min_gap: int,
    mut_collisions: str,
    guess_p_mut_given_span: np.ndarray | None = None,
    guess_p_ends: np.ndarray | None = None,
    guess_p_clust: np.ndarray | None = None,
    *,
    prenormalize: bool = True,
    max_iter: int = 128,
    convergence_thresh: float = 1.0e-4,
    **kwargs,
):
    """Calculate the three sets of parameters based on observed data.

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
    mut_collisions: str
        Method for handling reads with two mutations that are too close;
        must be either "drop" or "merge".
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
    import numpy as np
    from .unbias_jit import (
        adjust_min_gap_jit,
        clip_jit,
        normalize_jit,
        calc_p_noclose_given_ends_jit,
        calc_p_nomut_window_jit,
    )
    # Validate the dimensions.
    dims = find_dims(
        [(POSITIONS, CLUSTERS), (POSITIONS, POSITIONS, CLUSTERS), (CLUSTERS,)],
        [p_mut_given_span_observed, p_ends_observed, p_clust_observed],
        ["p_mut_given_span_observed", "p_ends_observed", "p_clust_observed"],
        nonzero=True,
    )
    # Normalize the values.
    require_atleast("min_gap", min_gap, 0, classes=int)
    min_gap = adjust_min_gap_jit(dims[POSITIONS], min_gap)
    if prenormalize:
        p_mut_given_span_observed = clip_jit(p_mut_given_span_observed)
        p_ends_observed = triu_norm(clip_jit(p_ends_observed))
        p_clust_observed = normalize_jit(clip_jit(p_clust_observed))
    # Determine the initial guess for the mutation rates.
    if guess_p_mut_given_span is not None:
        require_isinstance("guess_p_mut_given_span", guess_p_mut_given_span, np.ndarray)
        require_equal(
            "guess_p_mut_given_span.shape",
            guess_p_mut_given_span.shape,
            p_mut_given_span_observed.shape,
            "p_mut_given_span_observed.shape",
        )
        if prenormalize:
            # Ensure the initial mutation rates are clipped.
            guess_p_mut_given_span = clip_jit(guess_p_mut_given_span)
    else:
        # If no initial guess was given, then use the mutation rates of
        # the observed reads.
        guess_p_mut_given_span = p_mut_given_span_observed
    # Determine the initial guess for the cluster proportions.
    if guess_p_clust is not None:
        require_isinstance("guess_p_clust", guess_p_clust, np.ndarray)
        require_equal(
            "guess_p_clust.shape",
            guess_p_clust.shape,
            p_clust_observed.shape,
            "p_clust_observed.shape",
        )
        if prenormalize:
            # Ensure the initial cluster proportions sum to 1.
            guess_p_clust = normalize_jit(clip_jit(guess_p_clust))
    else:
        # If no initial guess was given, then use the proportions of the
        # observed reads.
        guess_p_clust = p_clust_observed
    # Determine the initial guess for the read coordinate distribution.
    if guess_p_ends is not None:
        require_isinstance("guess_p_ends", guess_p_ends, np.ndarray)
        require_equal(
            "guess_p_ends.shape",
            guess_p_ends.shape,
            p_ends_observed.shape[:2],
            "p_ends_observed.shape[:2]",
        )
        if prenormalize:
            # Ensure the initial coordinate distributions are normalized.
            guess_p_ends = triu_norm(clip_jit(guess_p_ends))
    else:
        # If no initial guess was given, then use the coordinates of
        # the observed reads.
        if dims[CLUSTERS] == 1:
            guess_p_ends = p_ends_observed[:, :, 0]
        else:
            guess_p_ends = np.average(p_ends_observed, axis=2, weights=guess_p_clust)
    if mut_collisions == MUT_COLLISIONS_DROP:
        # Iteratively update the mutation rates and read coordinates.
        require_atleast("max_iter", max_iter, 1, classes=int)
        for _ in range(max_iter):
            # Update the mutation rates using the read coordinates.
            next_p_mut_given_span = calc_p_mut_given_span(
                p_mut_given_span_observed,
                guess_p_ends,
                min_gap,
                mut_collisions,
                guess_p_mut_given_span,
                **kwargs,
            )
            # Compute the RMSD change in mutation rates.
            rmsd_p_mut_given_span = np.sqrt(
                np.mean(np.square(next_p_mut_given_span - guess_p_mut_given_span))
            )
            # Update guess_p_mut_given_span for the next iteration.
            guess_p_mut_given_span = next_p_mut_given_span
            # Check for convergence using the RMSD change.
            if rmsd_p_mut_given_span <= convergence_thresh:
                break
            # Compute the probability that reads with each pair of end
            # coordinates would have no two mutations too close.
            p_noclose_given_ends = calc_p_noclose_given_ends_jit(
                guess_p_mut_given_span,
                calc_p_nomut_window_jit(guess_p_mut_given_span, min_gap),
            )
            # Update the distribution of read end coordinates and
            # clusters using the mutation rates.
            guess_p_ends = calc_p_ends(
                p_ends_observed,
                p_noclose_given_ends,
                guess_p_mut_given_span,
                guess_p_clust,
            )
            guess_p_clust = calc_p_clust(
                p_clust_observed,
                calc_p_noclose_given_clust(guess_p_ends, p_noclose_given_ends),
            )
        else:
            logger.warning(
                "Mutation rates and distribution of end coordinates "
                "failed to converge in {} iterations",
                max_iter,
            )
    else:
        # If reads are not being dropped, then only the mutation rates
        # need to be updated.
        guess_p_mut_given_span = calc_p_mut_given_span(
            p_mut_given_span_observed,
            guess_p_ends,
            min_gap,
            mut_collisions,
            guess_p_mut_given_span,
            **kwargs,
        )
    return guess_p_mut_given_span, guess_p_ends, guess_p_clust


def calc_p_ends_given_clust_noclose(
    p_ends: np.ndarray, p_noclose_given_ends: np.ndarray
):
    """Calculate the proportion of reads with no two mutations too
    close with each pair of 5' and 3' coordinates.

    Assumptions
    -----------
    -   `p_ends` has 2 dimensions: (positions x positions)
    -   Every value in the upper triangle of `p_ends` is ≥ 0 and ≤ 1;
        no values below the main diagonal are used.
    -   The upper triangle of `p_ends` sums to 1.
    -   `p_noclose_given_ends` has 3 dimensions:
        (positions x positions x clusters)
    -   Every value in `p_noclose_given_ends` is ≥ 0 and ≤ 1.
    -   There is at least 1 cluster.

    Parameters
    ----------
    p_ends: np.ndarray
        2D (positions x positions) array of the proportion of reads in
        each cluster beginning at the row position and ending at the
        column position.
    p_noclose_given_ends: np.ndarray
        3D (positions x positions x clusters) array of the probabilities
        that a read with 5' and 3' coordinates corresponding to the row
        and column would have no two mutations too close.

    Returns
    -------
    np.ndarray
        3D (positions x positions x clusters) array of the proportion of
        reads without mutations too close, beginning at the row position
        and ending at the column position, in each cluster.
    """
    import numpy as np
    from .unbias_jit import triu_mul_jit
    # Validate the dimensions of the arguments.
    find_dims(
        [(POSITIONS, POSITIONS), (POSITIONS, POSITIONS, CLUSTERS)],
        [p_ends, p_noclose_given_ends],
        ["p_ends", "p_noclose_given_ends"],
        nonzero=True,
    )
    # Calculate the proportion of total reads that would have each
    # pair of end coordinates.
    return triu_norm(triu_mul_jit(p_ends[:, :, np.newaxis], p_noclose_given_ends))


def calc_p_ends_given_noclose(
    p_ends_given_clust_noclose: np.ndarray, p_clust_given_noclose: np.ndarray
):
    """Calculate the probability that a read would have each pair of
    5'/3' ends and no two mutations too close.

    Parameters
    ----------
    p_ends_given_clust_noclose: np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a read from each cluster has each pair of 5'/3' ends given
        that it has no two mutations too close.
    p_clust_given_noclose: np.ndarray
        1D (clusters) array of the probability that a read comes from
        each cluster given that it has no two mutations too close.

    Returns
    -------
    np.ndarray
        2D (positions x positions) array of the probability that a read
        with no two mutations too close has each pair of 5'/3' ends,
        regardless of the cluster to which it belongs.
    """
    # Validate the dimensions of the arguments.
    find_dims(
        [(POSITIONS, POSITIONS, CLUSTERS), (CLUSTERS,)],
        [p_ends_given_clust_noclose, p_clust_given_noclose],
        ["p_ends_given_clust_noclose", "p_clust_given_noclose"],
        nonzero=True,
    )
    return triu_norm(p_ends_given_clust_noclose @ p_clust_given_noclose)


def calc_p_clust_given_ends_noclose(
    p_ends_given_clust_noclose: np.ndarray, p_clust_given_noclose: np.ndarray
):
    """Calculate the probability that a read with each pair of 5'/3'
    ends and no two mutations too close came from each cluster.

    Parameters
    ----------
    p_ends_given_clust_noclose: np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a read from each cluster has each pair of 5'/3' ends given
        that it has no two mutations too close.
    p_clust_given_noclose: np.ndarray
        1D (clusters) array of the probability that a read comes from
        each cluster given that it has no two mutations too close.

    Returns
    -------
    np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a read with each pair of 5'/3' ends and no two mutations
        too close comes from each cluster.
    """
    import numpy as np
    from .unbias_jit import triu_div_jit
    # Validate the dimensions of the arguments.
    find_dims(
        [(POSITIONS, POSITIONS, CLUSTERS), (CLUSTERS,)],
        [p_ends_given_clust_noclose, p_clust_given_noclose],
        ["p_ends_given_clust_noclose", "p_clust_given_noclose"],
        nonzero=True,
    )
    p_ends_given_noclose = calc_p_ends_given_noclose(
        p_ends_given_clust_noclose, p_clust_given_noclose
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        # Apply Bayes' theorem to calculate p_clust_given_ends_noclose.
        return p_clust_given_noclose * triu_div_jit(
            p_ends_given_clust_noclose, p_ends_given_noclose[:, :, np.newaxis]
        )


def calc_p_ends_observed(
    npos: int,
    end5s: np.ndarray,
    end3s: np.ndarray,
    weights: np.ndarray | None = None,
    check_values: bool = True,
):
    """Calculate the proportion of each pair of 5'/3' end coordinates
    observed in `end5s` and `end3s`, optionally weighted by `weights`.

    Parameters
    ----------
    npos: int
        Number of positions.
    end5s: np.ndarray
        5' ends (0-indexed) of the reads: 1D array (reads)
    end3s: np.ndarray
        3' ends (0-indexed) of the reads: 1D array (reads)
    weights: np.ndarray | None
        Number of times each read occurs in each cluster:
        2D array (reads x clusters)
    check_values: bool
        Check that `end5s`, `end3s`, and `weights` are all valid.

    Returns
    -------
    np.ndarray
        Fraction of reads with each 5' (row) and 3' (column) coordinate:
        3D array (positions x positions x clusters)
    """
    import numpy as np
    from .unbias_jit import calc_p_ends_observed_jit
    require_atleast("npos", npos, 1, classes=int)
    # Validate the dimensions.
    if weights is None:
        # Assume all reads are equally likely and there is one cluster.
        return calc_p_ends_observed(
            npos, end5s, end3s, weights=np.ones_like(end5s), check_values=check_values
        )
    require_isinstance("weights", weights, np.ndarray)
    if weights.ndim == 1:
        # There is one cluster: return a 2D array for that cluster.
        return calc_p_ends_observed(
            npos,
            end5s,
            end3s,
            weights=weights[:, np.newaxis],
            check_values=check_values,
        )[:, :, 0]
    find_dims(
        [(READS,), (READS,), (READS, CLUSTERS)],
        [end5s, end3s, weights],
        ["end5s", "end3s", "weights"],
    )
    if check_values:
        # Validate the values.
        require_atleast("min(end5s)", end5s.min(initial=0), 0)
        require_atmost("max(end5s)", end5s.max(initial=0), npos - 1, "npos - 1")
        require_atleast("min(end3s)", end3s.min(initial=0), 0)
        require_atmost("max(end3s)", end3s.max(initial=0), npos - 1, "npos - 1")
        require_atleast("min(weights)", weights.min(initial=0.0), 0.0)
    # Call the compiled function.
    return triu_norm(calc_p_ends_observed_jit(npos, end5s, end3s, weights))


def calc_n_reads_per_pos(p_ends_observed: np.ndarray, n_reads_per_clust: np.ndarray):
    """Calculate the number of reads covering each position per cluster.

    Parameters
    ----------
    p_ends_observed: np.ndarray
        Observed proportion of reads with each pair of 5'/3' end
        coordinates: 3D array (positions x positions x clusters).
    n_reads_per_clust: np.ndarray
        Total number of reads in each cluster: 1D array (clusters).

    Returns
    -------
    np.ndarray
        Number of reads covering each position per cluster:
        2D array (positions x clusters).
    """
    find_dims(
        [(POSITIONS, POSITIONS, CLUSTERS), (CLUSTERS,)],
        [p_ends_observed, n_reads_per_clust],
        ["p_ends_observed", "n_reads_per_clust"],
    )
    return calc_rectangular_sum(p_ends_observed) * n_reads_per_clust


def calc_params_observed(
    n_pos_total: int,
    unmasked_pos: Iterable[int],
    muts_per_pos: Iterable[np.ndarray],
    end5s: np.ndarray,
    end3s: np.ndarray,
    counts_per_uniq: np.ndarray,
    resps: np.ndarray,
):
    """Calculate the observed estimates of the parameters.

    Parameters
    ----------
    n_pos_total: int
        Total number of positions in the region.
    unmasked_pos: Iterable[int]
        Unmasked positions; must be zero-indexed with respect to the
        5' end of the region.
    muts_per_pos: Iterable[np.ndarray]
        For each unmasked position, numbers of all reads with a mutation
        at that position.
    end5s: np.ndarray
        5' end of every unique read; must be 0-indexed with respect to
        the 5' end of the region.
    end3s: np.ndarray
        3' end of every unique read; must be 0-indexed with
        respect to the 5' end of the region.
    counts_per_uniq: np.ndarray
        Number of times each unique read occurs.
    resps: np.ndarray
        Cluster memberships of each read: 2D array (reads x clusters)

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    import numpy as np
    dims = find_dims(
        [(READS,), (READS,), (READS,), (READS, CLUSTERS)],
        [end5s, end3s, counts_per_uniq, resps],
        ["end5s", "end3s", "counts_per_uniq", "resps"],
    )
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
    p_ends_observed = calc_p_ends_observed(
        n_pos_total, end5s, end3s, p_each_read_each_clust, check_values=False
    )
    # Count the observed reads that cover each position.
    # 2D (all positions x clusters)
    n_reads_per_pos = calc_n_reads_per_pos(p_ends_observed, n_reads_per_clust)
    # Count the observed mutations at each position.
    # 2D (all positions x clusters)
    n_muts_per_pos = np.zeros((n_pos_total, dims[CLUSTERS]))
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
