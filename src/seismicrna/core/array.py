from typing import Iterable

import numpy as np
import pandas as pd

from .types import UINT_NBYTES, fit_uint_type, get_uint_type


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
        raise TypeError(f"array must be ndarray, but got {type(array).__name__}")
    if array.ndim != 1:
        raise ValueError(f"{what} must have 1 dimension, but got {array.ndim}")
    length, = array.shape
    return length


def get_inverse(target: np.ndarray, what: str = "array"):
    """ Map integers in [0, max(target)] to their 0-based indexes in
    `target`, or to -1 if not in `target`. """
    uniq, counts = np.unique(target, return_counts=True)
    if uniq.size > 0:
        # Verify that all values in target are non-negative.
        if get_length(uniq, what) > 0 > uniq[0]:
            raise ValueError(f"{what} has negative values: {uniq[uniq < 0]}")
        # Verify that all values in target are unique.
        if np.max(counts) > 1:
            raise ValueError(f"{what} has repeated values: {uniq[counts > 1]}")
        # Initialize a 1-dimensional array whose length is one more than
        # the maximum value of target, so that the array has every index
        # in the range [0, max(target)].
        inverse_size = np.max(target) + 1
    else:
        inverse_size = 0
    # Initialize all elements in the array to -1, a placeholder value.
    inverse = np.full(inverse_size, -1)
    # For each value n with index i in target, set the value at index n
    # of inverse to i.
    inverse[target] = np.arange(get_length(target, what))
    return inverse


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
