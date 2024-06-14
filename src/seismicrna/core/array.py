from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from numba import jit

from .types import UINT_NBYTES, fit_uint_type, get_uint_type

rng = np.random.default_rng()

MISSING = -1


def _unpack_tuple(items: Any):
    """ If items is a length-1 tuple, then return its single item;
    otherwise, return the items unchanged. """
    if isinstance(items, tuple) and len(items) == 1:
        return items[0]
    return items


def list_naturals(n: int):
    """ List natural numbers up to and including `n`. """
    return np.arange(1, n + 1)


def check_naturals(values: np.ndarray, what: str = "values"):
    """ Raise ValueError if the values are not monotonically increasing
    natural numbers. """
    length = get_length(values, what)
    if not np.array_equal(values, np.arange(1, length + 1)):
        raise ValueError(f"{what} must be numbered 1 to {length}, "
                         f"but got {values}")
    return np.asarray(values, dtype=int)


def get_length(array: np.ndarray, what: str = "array") -> int:
    if not isinstance(array, np.ndarray):
        raise TypeError(
            f"{what} must be ndarray, but got {type(array).__name__}"
        )
    if array.ndim != 1:
        raise ValueError(f"{what} must have 1 dimension, but got {array.ndim}")
    length, = array.shape
    return length


@jit()
def _fill_inverse_fwd(inverse: np.ndarray, default: int):
    """ Fill missing indexes in `inverse` in forward order. """
    fill = default
    for i in range(inverse.size):
        inv = inverse[i]
        if inv == MISSING:
            inverse[i] = fill
        else:
            fill = inv


@jit()
def _fill_inverse_rev(inverse: np.ndarray, default: int):
    """ Fill missing indexes in `inverse` in reverse order. """
    fill = default
    for i in range(inverse.size - 1, -1, -1):
        inv = inverse[i]
        if inv == MISSING:
            inverse[i] = fill
        else:
            fill = inv


def calc_inverse(target: np.ndarray,
                 require: int = -1,
                 fill: bool = False,
                 fill_rev: bool = False,
                 fill_default: int | None = None,
                 verify: bool = True,
                 what: str = "array"):
    """ Calculate the inverse of `target`, such that if element i of
    `target` has value x, then element x of the inverse has value i.

    >>> list(calc_inverse(np.array([3, 2, 7, 5, 1])))
    [-1, 4, 1, 0, -1, 3, -1, 2]
    >>> list(calc_inverse(np.arange(5)))
    [0, 1, 2, 3, 4]

    Parameters
    ----------
    target: np.ndarray
        Target values; must be a 1-dimensional array of non-negative
        integers with no duplicate values.
    require: int = -1
        Require the inverse to contain all indexes up to and including
        `require` (i.e. that its length is at least `require` + 1);
        ignored if `require` is -1; must be ≥ -1.
    fill: bool = False
        Fill missing indexes (that do not appear in `target`).
    fill_rev: bool = False
        Fill missing indexes in reverse order instead of forward order;
        only used if `fill` is True.
    fill_default: int | None = None
        Value with which to fill before the first non-missing value has
        been encountered; if `fill_rev` is True, defaults to the length
        of `target`, otherwise to -1.
    verify: bool = True
        Verify that all target values are unique, non-negative integers.
        If this is incorrect, then if `verify` is True, then ValueError
        will be raised; and if False, then the results of this function
        will be incorrect. Always set to True unless you have already
        verified that `target` is unique, non-negative integers.
    what: str = "array"
        What to name the array (only used for error messages).

    Returns
    -------
    np.ndarray
        Inverse of `target`.
    """
    length = get_length(target, what)
    if length == 0:
        # If target is empty, then return an empty inverse.
        return np.full(require + 1, -1)
    if verify:
        uniq, counts = np.unique(target, return_counts=True)
        if uniq.size > 0:
            # Verify that all values in target are non-negative.
            if get_length(uniq, what) > 0 > uniq[0]:
                raise ValueError(
                    f"{what} has negative values: {uniq[uniq < 0]}"
                )
            # Verify that all values in target are unique.
            if counts.max() > 1:
                raise ValueError(
                    f"{what} has repeated values: {uniq[counts > 1]}"
                )
    # Create a 1-dimensional array whose length is one greater than the
    # maximum value of target, so that the array has every index in the
    # range [0, max(target)]; initialize all elements to be missing.
    max_value = max(target.max(), require)
    inverse = np.full(max_value + 1, MISSING)
    # For each value n with index i in target, set the value at index n
    # of inverse to i.
    inverse[target] = np.arange(length)
    if fill:
        # Fill missing values in inverse.
        if fill_rev:
            _fill_inverse_rev(inverse,
                              (fill_default
                               if fill_default is not None
                               else length))
        else:
            _fill_inverse_fwd(inverse,
                              (fill_default
                               if fill_default is not None
                               else -1))
    return inverse


def locate_elements(collection: np.ndarray,
                    *elements: np.ndarray,
                    what: str = "collection",
                    verify: bool = True):
    """ Find the index at which each element of `elements` occurs in
    `collection`.

    >>> list(locate_elements(np.array([4, 1, 2, 7, 5, 3]), np.array([5, 2, 5])))
    [4, 2, 4]

    Parameters
    ----------
    collection: np.ndarray
        Collection in which to find each element in `elements`; must be
        a 1-dimensional array of non-negative integers with no duplicate
        values.
    *elements: np.ndarray
        Elements to find; must be a 1-dimensional array that is a subset
        of `collection`, although duplicate values are permitted.
    what: str = "collection"
        What to name the collection (only used for error messages).
    verify: bool = True
        Verify that all values in `collection` are unique, non-negative
        integers and that all items in `elements` are in `collections`.

    Returns
    -------
    np.ndarray
        Index of each element of `elements` in `collections`.
    """
    if verify:
        for e in elements:
            if get_length(extras := e[np.isin(e, collection, invert=True)]):
                raise ValueError(f"Elements {extras} are not in {what}")
    inverse = calc_inverse(collection, verify=verify, what=what)
    return _unpack_tuple(tuple(inverse[e] for e in elements))


def ensure_same_length(arr1: np.ndarray,
                       arr2: np.ndarray,
                       what1: str = "array1",
                       what2: str = "array2"):
    if (len1 := get_length(arr1, what1)) != (len2 := get_length(arr2, what2)):
        raise ValueError(
            f"Lengths differ between {what1} ({len1}) and {what2} ({len2})"
        )
    return len1


def ensure_order(array1: np.ndarray,
                 array2: np.ndarray,
                 what1: str = "array1",
                 what2: str = "array2",
                 gt_eq: bool = False):
    """ Ensure that `array1` is ≤ or ≥ `array2`, element-wise.

    Parameters
    ----------
    array1: np.ndarray
        Array 1 (same length as `array2`).
    array2: np.ndarray
        Array 2 (same length as `array1`).
    what1: str = "array1"
        What `array1` contains (only used for error messages).
    what2: str = "array2"
        What `array2` contains (only used for error messages).
    gt_eq: bool = False
        Ensure `array1 ≥ array2` if True, otherwise `array1 ≤ array2`.

    Returns
    -------
    int
        Shared length of `array1` and `array2`.
    """
    length = ensure_same_length(array1, array2, what1, what2)
    ineq_func, ineq_sign = (np.less, "<") if gt_eq else (np.greater, ">")
    if np.any(is_err := ineq_func(array1, array2)):
        index = pd.Index(np.arange(length)[is_err])
        errors = pd.DataFrame.from_dict(
            {what1: pd.Series(array1[is_err], index=index),
             what2: pd.Series(array2[is_err], index=index)}
        )
        raise ValueError(f"Got {what1} {ineq_sign} {what2}:\n{errors}")
    return length


def sanitize_values(values: Iterable[int],
                    lower_limit: int,
                    upper_limit: int,
                    whats: str = "values"):
    """ Validate and sort values, and return them as an array. """
    # Convert the values to an array and ensure it is one-dimensional.
    if not isinstance(values, (np.ndarray, list)):
        values = list(values)
    array = np.asarray(values)
    n_values = get_length(array, whats)
    if n_values == 0:
        # The array is empty.
        return np.array([], dtype=get_uint_type(min(UINT_NBYTES)))
    # Find and sort the unique values.
    array = np.unique(array)
    if array.size != n_values:
        raise ValueError(f"Got non-unique {whats}")
    # Validate the minimum and maximum values.
    min_value = array[0]
    max_value = array[-1]
    if min_value < lower_limit:
        raise ValueError(f"All {whats} must be ≥ {lower_limit}, but got "
                         f"{array[array < lower_limit]}")
    if max_value > upper_limit:
        raise ValueError(f"All {whats} must be ≤ {upper_limit}, but got "
                         f"{array[array > upper_limit]}")
    # Return the array as the smallest data type that will fit the data.
    return np.asarray(array, dtype=fit_uint_type(max_value))


def stochastic_round(values: np.ndarray):
    """ Round values to integers stochastically, so that the probability
    of rounding up equals the mantissa. """
    ints = np.floor(values)
    mantissas = values - ints
    return np.asarray(ints, dtype=int) + (rng.random(values.shape) < mantissas)


def find_dims(dims: Sequence[Sequence[str | None]],
              arrays: Sequence[np.ndarray],
              names: Sequence[str] | None = None,
              nonzero: Iterable[str] | bool = False):
    """ Check the dimensions of the arrays.

    Parameters
    ----------

    """
    # Ensure that nonzero is either True or a set of str.
    if nonzero is False:
        nonzero = set()
    elif nonzero is not True:
        nonzero = set(map(str, nonzero))
    # Verify there are the same number of arrays, dimensions, and names.
    if (n := len(arrays)) != len(dims):
        raise ValueError("The numbers of arrays and dimensions must equal, "
                         f"but got {n} array(s) and {len(dims)} dimension(s)")
    if names is not None:
        if len(names) != n:
            raise ValueError("The numbers of arrays and names must equal, "
                             f"but got {n} array(s) and {len(names)} name(s)")
    else:
        names = [f"array{i}" for i in range(n)]
    # Check the dimensions of the arrays.
    sizes = dict()
    for array, dim, name in zip(arrays, dims, names, strict=True):
        if not isinstance(array, np.ndarray):
            raise TypeError(f"Each array must be a NumPy NDArray, "
                            f"but got {type(array).__name__} for {repr(name)}")
        # Count the named and extra dimensions for this array.
        if len(dim) > 0 and dim[-1] is None:
            # The last dimension is None, so subtract it from the number
            # of named dimensions.
            n_named = len(dim) - 1
            # Extra dimensions in the array are allowed.
            extras = True
        else:
            # All dimensions are named.
            n_named = len(dim)
            # Extra dimensions in the array are forbidden.
            extras = False
        # Verify the array has a valid number of dimensions.
        if array.ndim != n_named:
            if not extras:
                raise ValueError(f"Array {repr(name)} must have {n_named} "
                                 f"dimension(s), but got {array.ndim}")
            if array.ndim < n_named:
                raise ValueError(f"Array {repr(name)} must have ≥ {n_named} "
                                 f"dimension(s), but got {array.ndim}")
        # Check each named dimension of the array.
        for i in range(n_named):
            if not isinstance(dim[i], str):
                raise TypeError("The name of each dimension must be str, "
                                f"but got {type(dim[i]).__name__}")
            # Get the size of this dimension in the array.
            size = array.shape[i]
            if (other_size := sizes.get(dim[i])) is not None:
                # A dimension of this name was already encountered.
                if size != other_size:
                    raise ValueError("Got multiple sizes for dimension "
                                     f"{repr(dim[i])}: {other_size} ≠ {size}")
            else:
                # This is the first time this dimension was encountered.
                # Validate the size.
                if not isinstance(size, int):
                    raise TypeError(f"Size of dimension {repr(dim[i])} must "
                                    f"be int, but got {type(size).__name__}")
                min_size = int(nonzero is True or dim[i] in nonzero)
                if size < min_size:
                    raise ValueError(f"Size of dimension {repr(dim[i])} must "
                                     f"be ≥ {min_size}, but got {size}")
                sizes[dim[i]] = size
    # Check if any dimensions in nonzero were not defined.
    if nonzero is not True and (unknown := nonzero - set(sizes)):
        raise ValueError(f"Unknown dimensions for nonzero: {unknown}")
    # Return the size of each dimension.
    return sizes


# Use @jit() because triangular is called by other jitted functions.
@jit()
def triangular(n: int):
    """ The `n` th triangular number (`n` ≥ 0): number of items in an
    equilateral triangle with `n` items on each side.

    Parameters
    ----------
    n: int
        Index of the triangular number to return; equivalently, the side
        length of the equilateral triangle.

    Returns
    -------
    int
        The triangular number with index `n`; equivalently, the number
        of items in the equilateral triangle of side length `n`.
    """
    return (n * n + n) // 2


def find_true_dists(booleans: np.ndarray):
    """ Find the distance to each True element in a boolean array. """
    # Initialize the distances to the maximum (the length of booleans).
    length = get_length(booleans, "booleans")
    dists = np.full(length, length)
    # Locate each True (nonzero) element of booleans.
    trues = np.flatnonzero(booleans)
    if trues.size > 0:
        n2 = trues[0]
        dists[:n2] = np.arange(n2, 0, -1)
        for i in range(trues.size - 1):
            a = trues[i]
            b = trues[i + 1]
            d = b - a
            n2 = d // 2
            n1 = b - n2
            dists[a: n1] = np.arange(d - n2)
            dists[n1: b] = np.arange(n2, 0, -1)
        n1 = trues[-1]
        dists[n1:] = np.arange(length - n1)
    return dists
