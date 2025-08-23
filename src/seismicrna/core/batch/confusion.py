import numpy as np
import pandas as pd
from numba import jit

from ..array import intersect1d_unique_sorted
from ..logs import logger
from ..seq import POS_NAME
from ..validate import require_isinstance, require_between

POSITION_A = "Position A"
POSITION_B = "Position B"


def init_confusion_matrix(pos_index: pd.Index,
                          clusters: pd.Index | None = None,
                          min_gap: int = 0):
    """ For every pair of positions, initialize the confusion matrix:

    +----+----+
    | AB | AO | A.
    +----+----+
    | OB | OO | O.
    +----+----+
      .B   .O   ..

    And return .., A., .B, AB in that order.
    """
    positions = pos_index.get_level_values(POS_NAME)
    idx_a, idx_b = np.triu_indices(positions.size, k=1)
    pos_a = positions[idx_a]
    pos_b = positions[idx_b]
    select = pos_b - pos_a > min_gap
    pos_pairs = pd.MultiIndex.from_arrays((pos_a[select], pos_b[select]),
                                          names=[POSITION_A, POSITION_B])
    if clusters is None:
        # There are no clusters.
        array_type = pd.Series
        kwargs = dict(data=0, index=pos_pairs)
    elif isinstance(clusters, pd.Index):
        # There are clusters.
        array_type = pd.DataFrame
        kwargs = dict(data=0., index=pos_pairs, columns=clusters)
    else:
        raise TypeError(clusters)
    # Create the four elements of the confusion matrix.
    n = array_type(**kwargs)
    a = array_type(**kwargs)
    b = array_type(**kwargs)
    ab = array_type(**kwargs)
    return n, a, b, ab


@jit
def _count_intersection(x: np.ndarray, y: np.ndarray):
    """ Count how many elements are in both x and y, assuming x and y
    are both NumPy arrays where all elements are unique and sorted from
    smallest to largest. """
    nx = x.size
    ny = y.size
    i = 0
    j = 0
    count = 0
    while i < nx and j < ny:
        xi = x[i]
        yj = y[j]
        if xi < yj:
            i += 1
        elif xi > yj:
            j += 1
        else:
            count += 1
            i += 1
            j += 1
    return count


def _count_intersection_weighted(x: np.ndarray,
                                 y: np.ndarray,
                                 read_weights: np.ndarray):
    return read_weights[intersect1d_unique_sorted(x, y)].sum(axis=0)


def calc_confusion_matrix(pos_index: pd.Index,
                          covering_reads: dict[int, np.ndarray],
                          mutated_reads: dict[int, np.ndarray],
                          read_weights: pd.DataFrame | None = None,
                          min_gap: int = 0):
    """ For every pair of positions, calculate the confusion matrix:

    +----+----+
    | AB | AO | A.
    +----+----+
    | OB | OO | O.
    +----+----+
      .B   .O   ..

    And return .., A., .B, AB in that order.
    """
    logger.routine("Began calculating confusion matrix")
    for pos in pos_index.get_level_values(POS_NAME):
        # Confirm all mutated reads are also covered.
        assert np.all(np.searchsorted(covering_reads[pos],
                                      mutated_reads[pos])
                      < covering_reads[pos].size)
    if read_weights is not None:
        require_isinstance("read_weights", read_weights, pd.DataFrame)
        if read_weights.index.size > 0:
            max_read_num = int(read_weights.index.max())
        else:
            max_read_num = -1
        rwv = np.empty((max_read_num + 1, read_weights.columns.size),
                       dtype=float)
        rwv[read_weights.index] = read_weights
    else:
        rwv = None
    # Initialize the confusion matrix.
    n, a, b, ab = init_confusion_matrix(pos_index,
                                        (read_weights.columns
                                         if read_weights is not None
                                         else None),
                                        min_gap=min_gap)
    # Cache the values for faster access.
    nv = n.values
    av = a.values
    bv = b.values
    abv = ab.values
    # For each pair of positions, count the reads in each category.
    for i, (pos5, pos3) in enumerate(n.index):
        # This method of counting the intersection is faster than matrix
        # multiplication, even though it needs to loop over every pair.
        # For speed, use x.values[i] instead of x.at[(pos5, pos3)].
        if read_weights is None:
            nv[i] = _count_intersection(covering_reads[pos5],
                                        covering_reads[pos3])
            av[i] = _count_intersection(covering_reads[pos3],
                                        mutated_reads[pos5])
            bv[i] = _count_intersection(covering_reads[pos5],
                                        mutated_reads[pos3])
            abv[i] = _count_intersection(mutated_reads[pos5],
                                         mutated_reads[pos3])
        else:
            nv[i] = _count_intersection_weighted(covering_reads[pos5],
                                                 covering_reads[pos3],
                                                 rwv)
            av[i] = _count_intersection_weighted(covering_reads[pos3],
                                                 mutated_reads[pos5],
                                                 rwv)
            bv[i] = _count_intersection_weighted(covering_reads[pos5],
                                                 mutated_reads[pos3],
                                                 rwv)
            abv[i] = _count_intersection_weighted(mutated_reads[pos5],
                                                  mutated_reads[pos3],
                                                  rwv)
    logger.routine("Ended calculating confusion matrix")
    return n, a, b, ab


def validate_confusion_matrix(n: pd.Series | pd.DataFrame,
                              a: pd.Series | pd.DataFrame,
                              b: pd.Series | pd.DataFrame,
                              ab: pd.Series | pd.DataFrame):
    if isinstance(n, pd.Series):
        require_isinstance("a", a, pd.Series)
        require_isinstance("b", b, pd.Series)
        require_isinstance("ab", ab, pd.Series)
    elif isinstance(n, pd.DataFrame):
        require_isinstance("a", a, pd.DataFrame)
        require_isinstance("b", b, pd.DataFrame)
        require_isinstance("ab", ab, pd.DataFrame)
    if np.any(n < 0):
        raise ValueError("n < 0")
    if np.any(a < 0):
        raise ValueError("a < 0")
    if np.any(b < 0):
        raise ValueError("b < 0")
    if np.any(ab < 0):
        raise ValueError("a_and_b < 0")
    if np.any(ab < a + b - n):
        raise ValueError("a_and_b < a + b - n")
    if np.any(a > n):
        raise ValueError("a > n")
    if np.any(b > n):
        raise ValueError("b > n")
    if np.any(ab > a):
        raise ValueError("a_and_b > a")
    if np.any(ab > b):
        raise ValueError("a_and_b > b")


def calc_confusion_phi(n: pd.Series | pd.DataFrame,
                       a: pd.Series | pd.DataFrame,
                       b: pd.Series | pd.DataFrame,
                       ab: pd.Series | pd.DataFrame, *,
                       min_cover: int | float = 1,
                       tol: float = 1.e6,
                       validate: bool = True):
    """ Calculate the phi correlation coefficient for a 2x2 matrix.

    +----+----+
    | AB | AO | A.
    +----+----+
    | OB | OO | O.
    +----+----+
      .B   .O   ..

    where
      A. = AB + AO
      .B = AB + OB
      .. = A. + O. = .B + .O

    Parameters
    ----------
    n: pd.Series | pd.DataFrame
        Observations in total (..)
    a: pd.Series | pd.DataFrame
        Observations for which A is true, regardless of B (A.)
    b: pd.Series | pd.DataFrame
        Observations for which B is true, regardless of A (.B)
    ab: pd.Series | pd.DataFrame
        Observations for which A and B are both true (AB)
    min_cover: float
        Set phi values with coverage < min_cov to NaN
    tol: float
        Tolerance for phi coefficients outside [-1, 1]
    validate: bool
        Validate the confusion matrix first

    Returns
    -------
    float | np.ndarray | pd.Series | pd.DataFrame
        Phi correlation coefficient
    """
    if validate:
        validate_confusion_matrix(n, a, b, ab)
    if min_cover < 0:
        raise ValueError("min_cover < 0")
    if tol < 0:
        raise ValueError("tol < 0")
    with np.errstate(divide="ignore", invalid="ignore"):
        # Use logarithms to prevent numerical overflow.
        log_axb = np.log(a) + np.log(b)
        log_denom = (log_axb + np.log(n - a) + np.log(n - b)) / 2
        log_phi_pos = np.log(n) + np.log(ab) - log_denom
        log_phi_neg = log_axb - log_denom
        phi = np.exp(log_phi_pos) - np.exp(log_phi_neg)
    phi.mask(np.asarray(n < min_cover), inplace=True)
    # Ensure all phi values are in [-1, 1]; values may be slightly out
    # of bounds due to floating-point rounding, especially after the
    # logarithmic transformation.
    phi_lo = np.asarray(phi < -1.)
    if np.any(phi_lo):
        if np.any(phi < -1. - tol):
            raise ValueError("phi < -1")
        phi.mask(phi_lo, other=-1., inplace=True)
    phi_hi = np.asarray(phi > 1.)
    if np.any(phi_hi):
        if np.any(phi > 1. + tol):
            raise ValueError("phi > 1")
        phi.mask(phi_hi, other=1., inplace=True)
    return phi


def calc_confusion_pvals(n: pd.Series,
                         a: pd.Series,
                         b: pd.Series,
                         ab: pd.Series, *,
                         validate: bool = True):
    """ Calculate the p-value of each element of the confusion matrix
    using a two-sided Fisher exact test. """
    if validate:
        validate_confusion_matrix(n, a, b, ab)
    # Calculate the p-values of the overlap being ≤ (le) or ≥ (ge) the
    # observed overlap (ab).
    from scipy.stats import hypergeom
    dist = hypergeom(n, b, a)
    pvals_le = dist.cdf(ab)
    pvals_ge = dist.sf(ab - 1)
    # Approximate the overall two-sided p-value as 2x the smaller tail,
    # capped at 1.
    pvals = np.minimum(2.0 * np.minimum(pvals_le, pvals_ge), 1.0)
    return pvals


def label_significant_pvals(pvals: np.ndarray | pd.Series | pd.DataFrame,
                            alpha: float):
    """ Label the significant p-values using the Benjamini-Hochberg
    method. """
    require_between("alpha", alpha, 0., 1., inclusive=False, classes=float)
    require_isinstance("pvals", pvals, (np.ndarray, pd.Series, pd.DataFrame))
    if pvals.ndim == 2:
        # If pvals is 2D, then calculate each column separately.
        values = np.stack([label_significant_pvals(np.asarray(pvals)[:, i],
                                                   alpha)
                           for i in range(pvals.shape[1])],
                          axis=1)
        if isinstance(pvals, pd.DataFrame):
            return pd.DataFrame(values, pvals.index, pvals.columns)
        return values
    elif pvals.ndim != 1:
        raise ValueError(f"pvals must be 1D or 2D, but is {pvals.ndim}D")
    if np.any(np.logical_or(pvals < 0, pvals > 1)):
        raise ValueError("pvals are not all in [0, 1]")
    # Count and rank the valid p-values.
    is_valid = ~np.isnan(pvals)
    pvals_valid = pvals[is_valid]
    m = pvals_valid.size
    if m == 0:
        # No p-values are valid, so none are significant.
        return np.logical_and(False, pvals.astype(bool))
    ranks = np.empty_like(pvals_valid, dtype=int)
    ranks[np.argsort(pvals_valid, kind="stable")] = np.arange(1, m + 1)
    # Determine the largest p-value that is ≤ (alpha / m) * rank.
    pvals_less_or_equal = pvals_valid[pvals_valid <= (alpha / m) * ranks]
    if pvals_less_or_equal.size > 0:
        max_significant_pval = np.max(pvals_less_or_equal)
    else:
        max_significant_pval = np.nan
    # All p-values less than or equal to that value are significant.
    is_significant = pvals <= max_significant_pval
    is_significant[~is_valid] = False
    return is_significant
