import numpy as np
import pandas as pd

from ..array import get_length
from ..seq import BASE_NAME

# Indexes of read and batch numbers.
READ_NUM = "Read Number"
BATCH_NUM = "Batch Number"
RB_INDEX_NAMES = [BATCH_NUM, READ_NUM]


def list_batch_nums(num_batches: int):
    """ List the batch numbers. """
    return list(range(num_batches))


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
