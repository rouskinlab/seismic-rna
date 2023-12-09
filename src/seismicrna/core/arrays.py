from functools import wraps
from typing import Callable, Iterable

import numpy as np
import pandas as pd


def get_priority(arrays: Iterable[np.ndarray | pd.Series | pd.DataFrame]):
    """ Determine the highest-priority type among the given arrays:
    pandas.DataFrame > pandas.Series > np.ndarray

    Parameters
    ----------
    *arrays: np.ndarray | pd.Series | pd.DataFrame
        Arrays among which to find the highest priority.

    Returns
    -------
    type[np.ndarray | pd.Series | pd.DataFrame]
        Highest-priority type.
    """
    # Ensure arrays is a list-like object.
    arrays = list(arrays)
    # Check if there are any invalid types.
    unknown = [array for array in arrays
               if not isinstance(array, (np.ndarray,
                                         pd.Series,
                                         pd.DataFrame))]
    if unknown:
        raise TypeError("Cannot determine priority of array types "
                        f"{[type(array).__name__ for array in unknown]}")
    # Check if there are any arrays of each type in order of decreasing
    # priority.
    for array_type in (pd.DataFrame, pd.Series, np.ndarray):
        if any(isinstance(array, array_type) for array in arrays):
            # If so, then that type is the highest priority type.
            return array_type
    raise ValueError("Cannot determine priority of zero arrays")


def get_shared_index(indexes: Iterable[pd.Index | pd.MultiIndex]):
    """ Verify that all given indexes are equal, and return the first.

    Parameters
    ----------
    indexes: Iterable[pandas.Index | pandas.MultiIndex]
        Indexes to compare.

    Returns
    -------
    pandas.Index | pandas.MultiIndex
        The shared index.
    """
    # Ensure indexes is a list-like object.
    indexes = list(indexes)
    try:
        # Get the first index.
        index = indexes[0]
    except IndexError:
        raise ValueError("No indexes were given")
    if not isinstance(index, pd.Index):
        raise TypeError(f"Expected Index, but got {type(index).__name__}")
    # Ensure all indexes are identical.
    for i, other in enumerate(indexes[1:], start=1):
        if not isinstance(other, pd.Index):
            raise TypeError(f"Expected Index, but got {type(other).__name__}")
        if isinstance(index, pd.MultiIndex):
            if not isinstance(other, pd.MultiIndex):
                raise TypeError(
                    f"Expected MultiIndex, but got {type(other).__name__}"
                )
            if other.names != index.names or not other.equal_levels(index):
                raise ValueError(f"Indexes 0 and {i} differ: {index} ≠ {other}")
        else:
            if isinstance(other, pd.MultiIndex):
                raise TypeError(
                    f"Expected Index, but got {type(other).__name__}"
                )
        if other.name != index.name or not other.equals(index):
            raise ValueError(f"Indexes 0 and {i} differ: {index} ≠ {other}")
    return index


def get_shared_indexes(arrays: Iterable[np.ndarray | pd.Series | pd.DataFrame]):
    """ Verify that all given arrays have the same indexes, and return
    the indexes of the first array.

    Parameters
    ----------
    arrays: Iterable[numpy.ndarray | pandas.Series | pandas.DataFrame]
        Arrays to compare.

    Returns
    -------
    dict[str, pd.Index]
        The shared indexes.
    """
    # Ensure arrays is a list-like object.
    arrays = list(arrays)
    # Determine the highest-priority type of the arrays.
    array_type = get_priority(arrays)
    # Check if any arrays are not of that type.
    not_type = list(filter(lambda a: not isinstance(a, array_type), arrays))
    if not_type:
        raise TypeError(
            f"Expected every array to be {type(array_type).__name__}, "
            f"but got {not_type}"
        )
    # Determine the indexes of the array.
    if array_type is pd.DataFrame:
        return dict(index=get_shared_index(array.index for array in arrays),
                    columns=get_shared_index(array.columns for array in arrays))
    if array_type is pd.Series:
        return dict(index=get_shared_index(array.index for array in arrays))
    if array_type is np.ndarray:
        return dict()
    raise TypeError(f"Invalid type for array: {type(array_type).__name__}")


def get_shared_shape(arrays: Iterable[np.ndarray | pd.Series | pd.DataFrame]):
    """ Verify that all given arrays have the same shape, and return it.

    Parameters
    ----------
    arrays: Iterable[numpy.ndarray | pandas.Series | pandas.DataFrame]
        Arrays to compare.

    Returns
    -------
    tuple[int, ...]
        The shared shape.
    """
    # Determine the shape of each array.
    shapes = [tuple(array.shape) for array in arrays]
    # Determine the number of unique shapes.
    nuniq = len(set(shapes))
    if nuniq == 0:
        raise ValueError("Cannot determine shape of zero arrays")
    if nuniq > 1:
        raise ValueError(f"Got more than one shape for arrays: {shapes}")
    return shapes[0]


def promote_arrays(*arrays: np.ndarray | pd.Series | pd.DataFrame):
    """ Promote all given arrays to the highest priority.

    Parameters
    ----------
    *arrays: numpy.ndarray | pandas.Series | pandas.DataFrame
        Arrays to promote.

    Returns
    -------
    tuple[numpy.ndarray | pandas.Series | pandas.DataFrame, ...]
        Promoted arrays.
    """
    # Verify that all arrays have the same shape.
    _ = get_shared_shape(arrays)
    # Determine the type to which to promote the arrays.
    return_type = get_priority(arrays)
    # Determine the indexes of the arrays of the highest type.
    indexes = get_shared_indexes(filter(lambda a: isinstance(a, return_type),
                                        arrays))
    # Promote the arrays.
    if return_type is np.ndarray:
        # NumPy arrays cannot be created by giving data to np.ndarray(),
        # but rather to np.array() or np.asarray().
        constructor = np.asarray
    else:
        # Other types (Series, DataFrame) can be created by calling the
        # type itself.
        constructor = return_type
    return tuple(constructor(array, **indexes) for array in arrays)


def rearray(values: np.ndarray | pd.Series | pd.DataFrame,
            astype: np.ndarray | pd.Series | pd.DataFrame,
            broadcast: bool = False):
    """ Return the given values in a container of the given type and,
    if applicable, with the given indexes and columns.

    Parameters
    ----------
    values: numpy.ndarray | pandas.Series | pandas.DataFrame
        Values to copy into the new container.
    astype: numpy.ndarray | pandas.Series | pandas.DataFrame
        Type of the container to return; indexes and columns will be
        returned in the output as well.
    broadcast: bool = False
        Broadcast values to the shape of astype.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Values in the new container.
    """
    # Ensure that values is a NumPy array.
    values = np.asarray(values)
    if broadcast:
        # Broadcast values to the shape of astype.
        values = np.array(np.broadcast_to(values, astype.shape))
    if isinstance(astype, pd.DataFrame):
        if values.ndim == 2:
            # If the given values are 2D, then return a DataFrame.
            return pd.DataFrame(values, astype.index, astype.columns)
        if values.ndim == 1:
            # If the given values are 1D, then return a Series whose
            # index is the columns of the DataFrame.
            return pd.Series(values, astype.columns)
        if values.ndim == 0:
            # If the given values are 0D, then return a scalar array.
            return values
        raise ValueError(f"Cannot return a {values.ndim}-dimensional array "
                         f"as type {type(astype).__name__}")
    if isinstance(astype, pd.Series):
        if values.ndim == 1:
            # If the given values are 1D, then return a Series.
            return pd.Series(values, astype.index)
        if values.ndim == 0:
            # If the given values are 0D, then return a scalar array.
            return values
        raise ValueError(f"Cannot return a {values.ndim}-dimensional array "
                         f"as type {type(astype).__name__}")
    if isinstance(astype, np.ndarray):
        # Return the values as a NumPy array.
        return values
    raise TypeError(f"Invalid type for array: {type(astype).__name__}")


def _np_internal(func: Callable, broadcast: bool):
    """ Decorate a function that accepts one NumPy or Pandas array-like
    argument so that it processes the values (NumPy array) and rearrays
    them to the original type, optionally after broadcasting. """

    @wraps(func)
    def wrapper(array: np.ndarray | pd.Series | pd.DataFrame, *args, **kwargs):
        # Call the function with a NumPy array argument, then rearray
        # the output to the original type.
        return rearray(func(np.asarray(array), *args, **kwargs),
                       array,
                       broadcast=broadcast)

    return wrapper


def np_internal(func: Callable):
    """ Decorate a function that accepts one NumPy or Pandas array-like
    argument so that it processes the values (NumPy array) and rearrays
    them to the original type. """

    return _np_internal(func, False)


def np_internal_broadcast(func: Callable):
    """ Decorate a function that accepts one NumPy or Pandas array-like
    argument so that it processes the values (NumPy array) and rearrays
    them to the original type, after broadcasting. """

    return _np_internal(func, True)


def np2_internal(func: Callable):
    """ Decorate a function that accepts two NumPy or Pandas array-like
    arguments so that it processes the values (NumPy arrays) and returns
    an object with the same type as the input argument. """

    @wraps(func)
    def wrapper(array1: np.ndarray | pd.Series | pd.DataFrame,
                array2: np.ndarray | pd.Series | pd.DataFrame,
                *args,
                **kwargs):
        # Call the function with two NumPy array arguments.
        result = func(np.asarray(array1), np.asarray(array2), *args, **kwargs)
        # Promote the result to the highest-priority array type.
        result, _, _ = promote_arrays(result, array1, array2)
        return result

    return wrapper
