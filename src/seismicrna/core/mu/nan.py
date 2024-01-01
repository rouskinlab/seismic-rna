"""
Comparisons of arbitrary numbers of mutation rates.
"""

from functools import wraps
from typing import Callable

import numpy as np
import pandas as pd

from .dim import count_pos, counts_pos_consensus


def any_nan(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Boolean array of positions where any mutation rate is NaN.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray | pandas.Series
        Boolean array of positions where any mutation rate is NaN.
    """
    # Reduce np.isnan over all axes but the first axis (i.e. axis 0),
    # which is -- by convention -- the position, and thus return an
    # array that has the same length as the first axis of mus.
    if mus.ndim <= 1:
        # If there are 1 or fewer axes, then no non-positional axes
        # exist to reduce with np.any(). Compute isnan without reducing.
        return np.isnan(mus)
    # Otherwise, reduce over the non-positional axes with np.any().
    return np.any(np.isnan(mus),
                  axis=(1 if mus.ndim == 2 else tuple(range(1, mus.ndim))))


def no_nan(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Boolean array of positions where no mutation rate is NaN.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray | pandas.Series
        Boolean array of positions where no mutation rate is NaN.
    """
    return np.logical_not(any_nan(mus))


def remove_nan(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Remove positions at which any mutation rate is NaN.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    tuple[numpy.ndarray | pandas.Series | pandas.DataFrame, ...]
        Mutation rates without NaN values.
    """
    # List the 0-indexed positions.
    positions = np.arange(count_pos(mus))
    # Find the positions with no NaN values.
    pos_no_nan = positions[no_nan(mus)]
    # Return only those positions (taking from the positional axis, 0).
    return np.take(mus, pos_no_nan, axis=0)


def removes_nan(*mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Remove positions at which any mutation rate in any group is NaN.

    Parameters
    ----------
    *mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Groups of mutation rates; each can contain multiple sets as the
        columns of a multidimensional array.

    Returns
    -------
    tuple[numpy.ndarray | pandas.Series | pandas.DataFrame, ...]
        Mutation rates without NaN values.
    """
    # List the 0-indexed positions.
    positions = np.arange(counts_pos_consensus(*mus))
    # Find positions with no NaN values in any group of mutation rates.
    pos_no_nan = positions[np.logical_and.reduce(list(map(no_nan, mus)))]
    # Return only those positions (along axis 0) from each group.
    return tuple(np.take(mu, pos_no_nan, axis=0) for mu in mus)


def auto_remove_nan(func: Callable):
    """ Decorate a function with one positional argument of mutation
    rates so that it automatically removes positions with NaNs from the
    input argument (but, if while using the NaN-less input, the function
    produces any new NaNs, then those NaNs will be returned).

    Note that if @auto_reframe and @auto_remove_nan are used to decorate
    the same function, then auto_reframe should be the inner decorator.
    If auto_remove_nan is the inner decorator and removes any NaNs, then
    auto_reframe will attempt to broadcast the NaN-less axis 0 over the
    original (longer) axis 0. This operation would raise a ValueError
    or, worse, if the NaN-less axis 0 happened to have length 1, would
    still broadcast to the original axis, causing a silent bug.
    """

    @wraps(func)
    def wrapper(mus: np.ndarray | pd.Series | pd.DataFrame, *args, **kwargs):
        return func(remove_nan(mus), *args, **kwargs)

    return wrapper


def auto_removes_nan(func: Callable):
    """ Decorate a function with positional argument(s) of mutation
    rates so that it automatically removes positions with NaNs from the
    input argument (but, if while using the NaN-less input, the function
    produces any new NaNs, then those NaNs will be returned). """

    @wraps(func)
    def wrapper(*mus: np.ndarray | pd.Series | pd.DataFrame, **kwargs):
        return func(*removes_nan(*mus), **kwargs)

    return wrapper

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
