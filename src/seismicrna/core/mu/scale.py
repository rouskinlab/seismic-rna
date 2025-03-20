import numpy as np
import pandas as pd

from .dim import count_pos
from .frame import auto_reframe
from .nan import auto_remove_nan


@auto_remove_nan
@auto_reframe
def calc_quantile(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Calculate the mutation rate at a quantile, ignoring NaNs.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile to return from the mutation rates; must be in [0, 1].

    Returns
    -------
    float | numpy.ndarray | pandas.Series
        Value of the quantile from the mutation rates.
    """
    if not isinstance(quantile, float):
        # Although numpy.quantile supports array-like values for the
        # quantile argument, get_quantile does not because the result
        # would have one or more extra dimensions.
        raise TypeError("Expected quantile to be float, "
                        f"but got {type(quantile).__name__}")
    if count_pos(mus) == 0:
        # If there are no positions, then return an all-NaN array with
        # the same dimensions as the input but without axis 0, instead
        # of raising an error, which np.quantile would do.
        return np.full(mus.shape[1:], np.nan)
    # Return the quantile by reducing along axis 0.
    return np.quantile(mus, quantile, axis=0)


def normalize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize the mutation rates to a quantile, so that the value of
    the quantile is scaled to 1 and all other mutation rates are scaled
    by the same factor. If quantile is 0, then do not normalize.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile for normalizing the mutation rates; must be in [0, 1].

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Normalized mutation rates.
    """
    return mus / calc_quantile(mus, quantile) if quantile > 0. else mus


@auto_reframe
def winsorize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize and winsorize the mutation rates to a quantile so that
    all mutation rates greater than or equal to the mutation rate at the
    quantile are set to 1, and lesser mutation rates are normalized.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile for normalizing the mutation rates; must be in [0, 1].

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Normalized and winsorized mutation rates.
    """
    return np.clip(normalize(mus, quantile), 0., 1.)


@auto_remove_nan
@auto_reframe
def calc_ranks(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Rank the mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Ranks of the mutation rates.
    """
    # Initialize all ranks to -1.
    ranks = np.full_like(mus, -1, dtype=int)
    # Make coordinate arrays to index every element of ranks in order.
    coords = np.unravel_index(np.arange(ranks.size), ranks.shape)
    # Replace the coordinate array of the first axis with the indexes
    # in order of their ranks.
    coords[0][:] = np.argsort(mus, axis=0).reshape(coords[0].shape)
    # Fill the ranks.
    ranks[coords] = np.ravel(
        np.broadcast_to(np.expand_dims(np.arange(mus.shape[0]),
                                       axis=tuple(range(1, mus.ndim))),
                        ranks.shape)
    )
    # Confirm all ranks were filled.
    if np.any(ranks < 0):
        # This error should be impossible, but checking just in case.
        raise ValueError(f"Ranks cannot be negative, but got {ranks}")
    return ranks
