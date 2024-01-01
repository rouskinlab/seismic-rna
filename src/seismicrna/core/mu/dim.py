import numpy as np
import pandas as pd


def count_pos(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Count the positions in an array of mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    int
        Number of positions in the array of mutation rates.
    """
    if mus.ndim == 0:
        raise ValueError("A 0-D array has no positional axis")
    return mus.shape[0]


def counts_pos(*mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Count the positions in each array of mutation rates.

    Parameters
    ----------
    *mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Groups of mutation rates; each can contain multiple sets as the
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    tuple[int, ...]
        Number of positions in each array of mutation rates.
    """
    return tuple(map(count_pos, mus))


def counts_pos_consensus(*mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Find the number of positions in every array of mutation rates;
    every array must have the same number of positions.

    Parameters
    ----------
    *mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Groups of mutation rates; each can contain multiple sets as the
        columns of a multidimensional array.

    Returns
    -------
    int
        Number of positions in every array of mutation rates.
    """
    if not mus:
        raise ValueError("Finding the number of positions in every array "
                         "requires at least 1 array, but got 0 arrays")
    counts = list(set(counts_pos(*mus)))
    if len(counts) > 1:
        raise ValueError(f"Finding the number of positions in every array "
                         f"requires all arrays to have the same number of "
                         f"positions, but got {sorted(counts)}")
    return counts[0]

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
