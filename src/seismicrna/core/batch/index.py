from typing import Iterable

import numpy as np
import pandas as pd

from ..array import (ensure_order,
                     ensure_same_length,
                     get_length,
                     sanitize_values)
from ..seq import BASE_NAME
from ..types import fit_uint_type

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


def get_num_segments(ends: np.ndarray) -> int:
    """ Number of segments for the given end coordinates. """
    if not isinstance(ends, np.ndarray):
        raise TypeError(f"ends must be ndarray, but got {type(ends).__name__}")
    if ends.ndim != 2:
        raise ValueError(f"ends must have 2 dimensions, but got {ends.ndim}")
    _, num_ends = ends.shape
    num_segs, odd = divmod(num_ends, 2)
    if odd:
        raise ValueError(
            f"Number of end coordinates must be even, but got {num_ends}"
        )
    return num_segs


def split_ends(ends: np.ndarray):
    """ Split an array of end positions into 5' and 3' ends. """
    # Although num_ends will just be the 2nd dimension (axis 1) of ends,
    # use get_num_segments to validate the shape of ends.
    num_ends = get_num_segments(ends) * 2
    end5s = ends[:, np.arange(0, num_ends, 2)]
    end3s = ends[:, np.arange(1, num_ends, 2)]
    return end5s, end3s


def sanitize_pos(positions: Iterable[int], min_pos: int, max_pos: int):
    """ Validate and sort positions, and return them as an array. """
    return sanitize_values(positions, min_pos, max_pos, "positions")


def sanitize_ends(ends: list[tuple[list | np.ndarray, list | np.ndarray]],
                  min_pos: int,
                  max_pos: int,
                  check_values: bool = True):
    """ Sanitize end coordinates.

    Parameters
    ----------
    ends: list[tuple[list | np.ndarray, list | np.ndarray]]
        End coordinates (1-indexed). Each item corresponds to a segment
        of the read and must be a tuple of its 5' and 3' coordinates.
    min_pos: int
        Minimum allowed value of a position.
    max_pos: int
        Maximum allowed value of a position.
    check_values: bool = True
        Whether to check the bounds of the values, which is the most
        expensive operation in this function. Can be set to False if the
        only desired effect is to ensure the output is a positive, even
        number of arrays in the proper data type.

    Returns
    -------
    list[tuple[np.ndarray, np.ndarray]]
        Sanitized end coordinates: positive, even number of arrays, all
        the same length and with the most efficient data type, and if
        `check_values` is True then all with valid coordinates.
    """
    if not ends:
        raise ValueError("Got no pairs of end coordinates")
    # Convert all end coordinates into arrays, if not already.
    ends = [(np.asarray(e5), np.asarray(e3)) for e5, e3 in ends]
    # Ensure all arrays have the same length.
    lengths = [get_length(ej, f"segment {i} {j}' coordinates")
               for i, seg in enumerate(ends)
               for j, ej in zip("53", seg, strict=True)]
    if len(set(lengths)) != 1:
        raise ValueError("Arrays of end coordinates must have equal lengths, "
                         f"but got lengths {lengths}")
    if check_values:
        # Validate every pair of end coordinates.
        length = lengths[0]
        for i, (e5, e3) in enumerate(ends):
            # All coordinates must be ≥ 1.
            ensure_order(e5,
                         np.broadcast_to(min_pos, length),
                         f"segment {i} 5' coordinates",
                         f"minimum position ({min_pos})",
                         gt_eq=True)
            # End 1 coordinates must be ≤ end 2 coordinates.
            ensure_order(e5,
                         e3,
                         f"segment {i} 5' coordinates",
                         f"segment {i} 3' coordinates")
            # All coordinates must be ≤ max position
            ensure_order(e3,
                         np.broadcast_to(max_pos, length),
                         f"segment {i} 3' coordinates",
                         f"maximum position ({max_pos})")
    # Convert all ends into arrays of the most efficient data type.
    # Conversion must be done after checking the values against max_pos.
    # Otherwise, values greater than the data type can accomodate would
    # be converted to their modulo relative to the maximum value of the
    # data type, which would mean checking the wrong values.
    pos_dtype = fit_uint_type(max_pos)
    return [(np.asarray(e5, pos_dtype), np.asarray(e3, pos_dtype))
            for e5, e3 in ends]


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
