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
                 target: np.ndarray | pd.Series | pd.DataFrame,
                 drop: int = 0):
    """ Place the values in an array object with the same type and axes
    as target.

    Parameters
    ----------
    values: Number | numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) to place in a new array-like object.
    target: numpy.ndarray | pandas.Series | pandas.DataFrame
        Array object whose type and axes are to be used for constructing
        the returned array.
    drop: int = 0
        Reduce the dimensionality of the target by dropping this number
        of axes, starting from axis 0 and continuing upwards.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Value(s) in their new array-like object.
    """
    # Determine the axes to pass to reframe based on the type of target.
    if isinstance(target, np.ndarray):
        axes = target.shape
    elif isinstance(target, pd.Series):
        axes = target.index,
    elif isinstance(target, pd.DataFrame):
        axes = target.index, target.columns
    else:
        raise TypeError("Expected target to be ndarray, Series, or Dataframe, "
                        f"but got {type(target).__name__}")
    # Optionally, drop axes, starting from axis 0.
    if drop < 0:
        raise ValueError(f"Cannot drop a negative number ({drop}) of axes")
    if drop > len(axes):
        raise ValueError(f"Cannot drop {drop} axes from a {len(axes)}-D array")
    # Reframe the values using the axes from the target.
    return reframe(values, axes[drop:])


def auto_reframe(func: Callable):
    """ Decorate a function with one positional argument of data so that
    it converts the input data to a NumPy array, runs, and then reframes
    the return value using the original argument as the target.

    Note that if @auto_reframe and @auto_remove_nan are used to decorate
    the same function, then auto_reframe should be the inner decorator.
    If auto_remove_nan is the inner decorator and removes any NaNs, then
    auto_reframe will attempt to broadcast the NaN-less axis 0 over the
    original (longer) axis 0. This operation would raise a ValueError
    or, worse, if the NaN-less axis 0 happened to have length 1, would
    still broadcast to the original axis, causing a silent bug.
    """

    @wraps(func)
    def wrapper(data: np.ndarray | pd.Series | pd.DataFrame, *args, **kwargs):
        # Compute the result of the function as a NumPy array.
        result = np.asarray(func(np.asarray(data), *args, **kwargs))
        # Reframe the result like the input argument, dropping any axes
        # that were eliminated by a reducing operation (e.g. summation).
        return reframe_like(result, data, data.ndim - result.ndim)

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
