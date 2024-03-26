from typing import Iterable

import numpy as np
import pandas as pd

from ..seq import BASE_NAME
from ..types import fit_uint_type, get_uint_type, UINT_NBYTES

# Missing identifiers.
NO_READ = -1

# Indexes of read and batch numbers.
READ_NUM = "Read Number"
BATCH_NUM = "Batch Number"
RB_INDEX_NAMES = [BATCH_NUM, READ_NUM]

# Indexes of read end coordinates.
END5_COORD = "5' End"
END3_COORD = "3' End"
END_COORDS = END5_COORD, END3_COORD


def list_batch_nums(num_batches: int):
    """ List the batch numbers. """
    return list(range(num_batches))


def get_length(array: np.ndarray, what: str = "array") -> int:
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
        index = pd.Index(np.arange(length)[is_err], name=READ_NUM)
        errors = pd.DataFrame.from_dict({what1: pd.Series(array1[is_err],
                                                          index=index),
                                         what2: pd.Series(array2[is_err],
                                                          index=index)})
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


def sanitize_pos(positions: Iterable[int], seq_length: int):
    """ Validate and sort positions, and return them as an array. """
    return sanitize_values(positions, 1, seq_length, "positions")


def has_mids(mid5s: np.ndarray | None, mid3s: np.ndarray | None):
    """ Whether the 5' and 3' middle coordinates exist. """
    if isinstance(mid5s, np.ndarray) and isinstance(mid3s, np.ndarray):
        # They exist if they are both arrays and have the same length.
        ensure_same_length(mid5s,
                           mid3s,
                           "5' middle coordinates",
                           "3' middle coordinates")
        return True
    if mid5s is None and mid3s is None:
        # They do not exist if they are both None.
        return False
    # All other cases are invalid.
    raise ValueError(
        f"Inconsistent 5' and 3' middle coordinates:\n{mid5s}\nand\n{mid3s}"
    )


def contiguous_mates(mid5s: np.ndarray | None, mid3s: np.ndarray | None):
    """ Return whether the two mates form a contiguous read. """
    ensure_same_length(mid5s,
                       mid3s,
                       "5' middle coordinates",
                       "3' middle coordinates")
    return np.less_equal(mid5s, mid3s + 1)


def sanitize_ends(max_pos: int,
                  end5s: list[int] | np.ndarray,
                  mid5s: list[int] | np.ndarray | None,
                  mid3s: list[int] | np.ndarray | None,
                  end3s: list[int] | np.ndarray):
    pos_dtype = fit_uint_type(max_pos)
    end5s = np.asarray(end5s, pos_dtype)
    end3s = np.asarray(end3s, pos_dtype)
    # Verify 5' end ≥ min position
    ensure_order(end5s,
                 np.broadcast_to(1, end5s.shape),
                 "5' end positions",
                 f"minimum position (1)",
                 gt_eq=True)
    # Verify 5' end ≤ 3' end
    ensure_order(end5s,
                 end3s,
                 "5' end positions",
                 "3' end positions")
    # Verify 3' end ≤ max position
    ensure_order(end3s,
                 np.broadcast_to(max_pos, end3s.shape),
                 "3' end positions",
                 f"maximum position ({max_pos})")
    if has_mids(mid5s, mid3s):
        # Verify 5' end ≤ 5' mid
        ensure_order(end5s, mid5s, "5' end positions", "5' middle positions")
        # Verify 5' mid ≤ 3' end
        ensure_order(mid5s, end3s, "5' middle positions", "3' end positions")
        # Verify 5' end ≤ 3' mid
        ensure_order(end5s, mid3s, "5' end positions", "3' middle positions")
        # Verify 3' mid ≤ 3' end
        ensure_order(mid3s, end3s, "3' middle positions", "3' end positions")
        # Check which reads are made of contiguous mates, i.e. the mates
        # overlap or at least abut with no gap between.
        is_contiguous = contiguous_mates(mid5s, mid3s)
        if np.all(is_contiguous):
            # If all mates are contiguous (which is common in short read
            # sequencing data and hence worth optimizing), nullify the
            # 5' and 3' middle coordinates to reduce space.
            mid5s = None
            mid3s = None
        else:
            # For contiguous mates, set the 5' and 3' middle coordinates
            # to the 5' and 3' ends of the contiguous region they cover.
            mid5s[is_contiguous] = end5s[is_contiguous]
            mid3s[is_contiguous] = end3s[is_contiguous]
    return end5s, mid5s, mid3s, end3s


def stack_end_coords(end5s: np.ndarray, end3s: np.ndarray):
    """ Return the 5' and 3' ends as one 2D array. """
    ensure_same_length(end5s, end3s, "5' end coordinates", "3' end coordinates")
    return np.stack([end5s, end3s], axis=1)


def count_base_types(base_pos_index: pd.Index):
    """ Return the number of each type of base in the index of positions
    and bases. """
    base_types, counts = np.unique(base_pos_index.get_level_values(BASE_NAME),
                                   return_counts=True)
    return pd.Series(counts, base_types)


def iter_base_types(base_pos_index: pd.Index):
    """ For each type of base in the index of positions and bases, yield
    the positions in the index with that type of base. """
    bases, inverse = np.unique(base_pos_index.get_level_values(BASE_NAME),
                               return_inverse=True)
    for i, base in enumerate(bases):
        yield base, base_pos_index[inverse == i]


def iter_windows(pos_nums: np.ndarray, size: int):
    """ Yield the positions in each window of size positions of the
    section. """
    if size < 1:
        raise ValueError(f"size must be ≥ 1, but got {size}")
    if get_length(pos_nums, "pos_nums") > 0:
        # Create a Series with the position numbers as its index.
        pos_series = pd.Series(True, pos_nums)
        min_pos = pos_series.index[0]
        max_pos = pos_series.index[-1]
        # Define the 5' and 3' ends of the window.
        win5 = min_pos
        win3 = min(win5 + (size - 1), max_pos)
        # Yield the positions in each window.
        while win3 <= max_pos:
            yield (win5, win3), pos_series.loc[win5: win3].index.values
            win5 += 1
            win3 += 1

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
