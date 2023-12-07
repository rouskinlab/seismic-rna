import numpy as np
import pandas as pd

from .scale import standardize


def calc_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
              mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Compute the root-mean-square deviation between two groups of
    mutation rates.

    Parameters
    ----------
    mus1: np.ndarray | pd.Series | pd.DataFrame
        First group of mutation rates; can contain multiple sets as the
        columns of a multidimensional array or DataFrame.
    mus2: np.ndarray | pd.Series | pd.DataFrame
        Second group of mutation rates; can contain multiple sets as the
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    np.ndarray | pd.Series | pd.DataFrame
        Standardized mutation rates.
    """
    # Compute the difference between the standardized mutation rates.
    diff = standardize(mus2) - standardize(mus1)
    # Return the root-mean-square distance from the line of best fit
    # (which equals the difference divided by the square root of 2),
    # ignoring NaN values.
    return np.sqrt(np.nanmean(diff * diff, axis=0) / 2.)
