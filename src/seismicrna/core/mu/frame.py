"""

"""

from functools import wraps
from numbers import Number
from typing import Callable

import numpy as np
import pandas as pd


def reframe(values: Number | np.ndarray | pd.Series | pd.DataFrame,
            axes: None | tuple[int | np.ndarray | pd.Index, ...] = None):
    """ Place the values in an array object with the given axes.

    Parameters
    ----------
    values: Number | numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) to place in a new array-like object.
    axes: None | tuple[int | numpy.ndarray | pandas.Index, ...] = None
        Axes of the new array-like object, specified as follows:

        - If None, then return just the values as a NumPy array.
        - If a tuple, then each element creates an axis as follows:

          - If an integer, then force the corresponding axis to be of
            that length.
          - If an array-like, then assign the axis a Pandas Index from
            the values in the element.

          Then, the array and index types are determined as follows:

          - If all elements are integers, then return a NumPy array in
            which the values are broadcast to the shape given by axes.
          - If at least one element is array-like, then return a Pandas
            object (a Series if axes has one item, a DataFrame if two).
          - If integers and array-like items are mixed, then replace
            each integer with a Pandas RangeIndex.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) in their new array-like object.
    """
    if axes is None:
        # No axes specified: return just the values as a NumPy array.
        return np.asarray(values)
    # Determine whether to interpret all axes as dimensions or indexes.
    lengths = list()
    indexes = list()
    for axis in axes:
        if isinstance(axis, int):
            # For the current axis, only the dimension is specified.
            if indexes:
                # If any previous axes were indexes, then all axes must
                # be indexes, so promote this axis to a RangeIndex.
                indexes.append(pd.RangeIndex(axis))
            else:
                # Otherwise, just specify the length of the axis.
                lengths.append(axis)
        elif isinstance(axis, (np.ndarray, pd.Index)):
            # For the current axis, an index is explicitly specified.
            if lengths and not indexes:
                # If this is the first axis with an index, but at least
                # one axis was already given (by its dimension only),
                # then promote all previously given axes to indexes.
                indexes.extend(map(pd.RangeIndex, lengths))
            indexes.append(axis)
        else:
            raise TypeError("Expected each axis to be int, ndarray, or Index, "
                            f"but got {type(axis).__name__}")
    # Determine the shape of the output array and broadcast the values.
    shape = tuple(idx.size for idx in indexes) if indexes else tuple(lengths)
    broadcast = np.broadcast_to(values, shape)
    num_indexes = len(indexes)
    if num_indexes == 0:
        # No indexes were specified, so just return the broadcast array.
        return np.array(broadcast)
    if num_indexes == 1:
        # Exactly one index was specified, so return a Series.
        return pd.Series(broadcast, index=indexes[0])
    if num_indexes == 2:
        # Exactly two indexes were specified, so return a DataFrame.
        return pd.DataFrame(broadcast, index=indexes[0], columns=indexes[1])
    raise ValueError("A Pandas object must have 1 or 2 axes, "
                     f"but got {num_indexes}")


def reframe_like(values: Number | np.ndarray | pd.Series | pd.DataFrame,
                 target: np.ndarray | pd.Series | pd.DataFrame):
    """ Place the values in an array object with the same type and axes
    as target.

    Parameters
    ----------
    values: Number | numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) to place in a new array-like object.
    target: numpy.ndarray | pandas.Series | pandas.DataFrame
        Array object whose type and axes are to be used for constructing
        the returned array.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) in their new array-like object.
    """
    if isinstance(target, np.ndarray):
        axes = target.shape
    elif isinstance(target, pd.Series):
        axes = target.index,
    elif isinstance(target, pd.DataFrame):
        axes = target.index, target.columns
    else:
        raise TypeError("Expected target to be ndarray, Series, or Dataframe, "
                        f"but got {type(target).__name__}")
    return reframe(values, axes)


def auto_reframe(func: Callable):
    """ Decorate a function with one positional argument of mutation
    rates so that it automatically reframes the return value to the
    input value. """

    @wraps(func)
    def wrapper(mus: np.ndarray | pd.Series | pd.DataFrame, *args, **kwargs):
        # First, compute the result of the function as a NumPy array.
        result = np.asarray(func(mus, *args, **kwargs))
        # Then, determine how to reframe the result.
        if result.ndim == mus.ndim:
            # If the result has the same number of dimensions as the
            # input argument, then reframe based on the input argument.
            return reframe_like(result, mus)



########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
