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
    input argument (but not from the return value). """

    @wraps(func)
    def wrapper(mus: np.ndarray | pd.Series | pd.DataFrame, *args, **kwargs):
        return func(remove_nan(mus), *args, **kwargs)

    return wrapper


def auto_removes_nan(func: Callable):
    """ Decorate a function with positional argument(s) of mutation
    rates so that it automatically removes positions with NaNs from the
    input argument (but not from the return value). """

    @wraps(func)
    def wrapper(*mus: np.ndarray | pd.Series | pd.DataFrame, **kwargs):
        return func(*removes_nan(*mus), **kwargs)

    return wrapper
