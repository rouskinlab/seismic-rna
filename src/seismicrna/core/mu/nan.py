"""
Comparisons of arbitrary numbers of mutation rates.
"""

import numpy as np
import pandas as pd

from .dim import get_common_npos


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
    # Otherwise, reduce the non-positional axes with np.any().
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


def without_nans(*mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Remove positions from axis 0 at which any mutation rate is NaN.

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
    # Generate an array of the positions.
    positions = np.arange(get_common_npos(*mus))
    # Find positions with no NaN values in any group of mutation rates.
    pos_no_nan = positions[np.logical_and.reduce(list(map(no_nan, mus)))]
    # Return only those positions from each group.
    return tuple(np.take(mu, pos_no_nan, axis=0) for mu in mus)
