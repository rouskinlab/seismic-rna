"""
Pairwise comparisons of mutation rates.
"""

from typing import Callable

import numpy as np
import pandas as pd

from .one import calc_rms, normalize
from ..arg import MUCOMP_DETERM, MUCOMP_PEARSON, MUCOMP_RMSD, MUCOMP_SPEARMAN
from ..seq import get_shared_index, get_windows


def calc_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
              mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square deviation between two groups of
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
    norm1 = normalize(mus1, 1.)
    norm2 = normalize(mus2, 1.)
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


def calc_pearson(mus1: np.ndarray | pd.Series | pd.DataFrame,
                 mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the Pearson correlation coefficient between two groups
    of mutation rates.

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
    float | np.ndarray | pd.Series
        Pearson correlation coefficient.
    """
    # Calculate the mean of each input over the first axis.
    mean1 = np.nanmean(mus1, axis=0)
    mean2 = np.nanmean(mus2, axis=0)
    # Calculate the difference of each element from the mean.
    diff1 = mus1 - mean1
    diff2 = mus2 - mean2
    # Calculate the variance of each dataset over the first axis.
    var1 = np.nansum(diff1 * diff1, axis=0)
    var2 = np.nansum(diff2 * diff2, axis=0)
    # Calculate the covariance of the datasets over the first axis.
    cov = np.nansum(diff1 * diff2, axis=0)
    # Calculate the Pearson correlation coefficient.
    return cov / np.sqrt(var1 * var2)


def calc_coeff_determination(mus1: np.ndarray | pd.Series | pd.DataFrame,
                             mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the coefficient of determination (a.k.a. R-squared)
    between two groups of mutation rates.

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
    float | np.ndarray | pd.Series
        Coefficient of determination.
    """
    # The coefficient of determination equals the Pearson r, squared.
    r = calc_pearson(mus1, mus2)
    return r * r


def get_mucomp(name: str):
    if name == MUCOMP_RMSD:
        return calc_rmsd
    if name == MUCOMP_PEARSON:
        return calc_pearson
    if name == MUCOMP_DETERM:
        return calc_coeff_determination
    if name == MUCOMP_SPEARMAN:
        return
    raise ValueError(f"Invalid comparison method: {repr(name)}")


def compare_sliding(mus1: pd.Series,
                    mus2: pd.Series,
                    mucomp: Callable,
                    window: int,
                    winmin: int):
    """ Compare two Series via sliding windows.
    """
    # Initialize an empty Series for the sliding comparison.
    values = pd.Series(np.nan, index=get_shared_index([mus1.index, mus2.index]))
    # Calculate the value of the comparison for each window.
    for center, (win1, win2) in get_windows(mus1,
                                            mus2,
                                            size=window,
                                            min_count=winmin):
        values.loc[center] = mucomp(win1, win2)
    return values
