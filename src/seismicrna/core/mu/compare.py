import numpy as np
import pandas as pd

from .scale import calc_rms


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
    # Normalize the mutation rates so the maximum of each group is 1.
    norm1 = mus1 / np.max(mus1, axis=0)
    norm2 = mus2 / np.max(mus2, axis=0)
    # Compute the root-mean-square mutation rate for each group.
    rms1 = calc_rms(norm1)
    rms2 = calc_rms(norm2)
    # Standardize the mutation rates so that the root-mean-square of
    # each group is 1, and then compute the difference.
    diff = norm1 / rms1 - norm2 / rms2
    # Compute the root-mean-square difference and restore the original
    # scale by dividing by the geometric mean of the standardization
    # coefficients.
    return np.sqrt(np.nanmean(diff * diff, axis=0) * (rms1 * rms2))
