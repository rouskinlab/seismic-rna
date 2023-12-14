"""
Pairwise comparisons of mutation rates.
"""

from typing import Callable

import numpy as np
import pandas as pd

from .nan import auto_removes_nan
from .scale import calc_rms, calc_ranks, normalize
from ..arg import MUCOMP_DETERM, MUCOMP_PEARSON, MUCOMP_RMSD, MUCOMP_SPEARMAN
from ..seq import get_shared_index, get_windows


@auto_removes_nan
def calc_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
              mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square deviation (RMSD) between two
    groups of mutation rates, ignoring NaNs.

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
    mus1 = normalize(mus1, 1.)
    mus2 = normalize(mus2, 1.)
    # Compute the root-mean-square mutation rate for each group.
    rms1 = calc_rms(mus1)
    rms2 = calc_rms(mus2)
    # Standardize the mutation rates so that the root-mean-square of
    # each group is 1, and then compute the difference.
    diff = mus1 / rms1 - mus2 / rms2
    # Compute the root-mean-square difference and restore the original
    # scale by multiplying by the geometric mean of the root-mean-square
    # mutation rates.
    return np.sqrt(np.mean(np.square(diff), axis=0) * (rms1 * rms2))


@auto_removes_nan
def calc_pearson(mus1: np.ndarray | pd.Series | pd.DataFrame,
                 mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the Pearson correlation coefficient between two groups
    of mutation rates, ignoring NaNs.

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
    mean1 = np.mean(mus1, axis=0)
    mean2 = np.mean(mus2, axis=0)
    # Calculate the difference of each element from the mean.
    diff1 = mus1 - mean1
    diff2 = mus2 - mean2
    # Calculate the variance of each dataset over the first axis.
    var1 = np.sum(diff1 * diff1, axis=0)
    var2 = np.sum(diff2 * diff2, axis=0)
    # Calculate the covariance of the datasets over the first axis.
    cov = np.sum(diff1 * diff2, axis=0)
    # Calculate the Pearson correlation coefficient.
    return cov / np.sqrt(var1 * var2)


def calc_coeff_determ(mus1: np.ndarray | pd.Series | pd.DataFrame,
                      mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the coefficient of determination (a.k.a. R-squared)
    between two groups of mutation rates, ignoring NaNs.

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
    # The coefficient of determination equals the squared Pearson r.
    return np.square(calc_pearson(mus1, mus2))


@auto_removes_nan
def calc_spearman(mus1: np.ndarray | pd.Series | pd.DataFrame,
                  mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the Spearman rank correlation coefficient between two
    groups of mutation rates, ignoring NaNs.

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
        Spearman rank correlation coefficient.
    """
    # The Spearman correlation is the Pearson correlation of the ranks.
    return calc_pearson(calc_ranks(mus1), calc_ranks(mus2))


def _get_comp_method(key: str):
    """ Get a comparison method based on its key. """
    if key == MUCOMP_RMSD:
        return calc_rmsd, "Root-Mean-Square Deviation", "RMSD"
    if key == MUCOMP_PEARSON:
        return calc_pearson, "Pearson Correlation Coefficient", "PCC"
    if key == MUCOMP_DETERM:
        return calc_coeff_determ, "Coefficient of Determination", "R-Squared"
    if key == MUCOMP_SPEARMAN:
        return calc_spearman, "Spearman Correlation Coefficient", "SCC"
    raise ValueError(f"Invalid method of comparison: {repr(key)}")


def get_comp_func(key: str) -> Callable:
    """ Get the function of a comparison method based on its key.

    Parameters
    ----------
    key: str
        Key with which to retrieve the comparison function.

    Returns
    -------
    Callable
        Function to compare mutation rates.
    """
    func, _, __ = _get_comp_method(key)
    return func


def get_comp_name(key: str) -> str:
    """ Get the name of a comparison method based on its key.

    Parameters
    ----------
    key: str
        Key with which to retrieve the comparison method name.

    Returns
    -------
    str
        Name of the comparison method.
    """
    _, name, __ = _get_comp_method(key)
    return name


def get_comp_abbr(key: str) -> str:
    """ Get the abbreviation of a comparison method based on its key.

    Parameters
    ----------
    key: str
        Key with which to retrieve the comparison method abbreviation.

    Returns
    -------
    str
        Abbreviation of the comparison method.
    """
    _, __, abbr = _get_comp_method(key)
    return abbr


def compare_windows(mus1: pd.Series,
                    mus2: pd.Series,
                    method: str | Callable,
                    size: int,
                    min_count: int = 2):
    """ Compare two Series via sliding windows.
    """
    if isinstance(method, str):
        # If the comparison method is given a string, then fetch the
        # function itself.
        method = get_comp_func(method)
    # Initialize an empty Series for the sliding comparison.
    values = pd.Series(np.nan, index=get_shared_index([mus1.index, mus2.index]))
    # Calculate the value of the comparison for each window.
    for center, (win1, win2) in get_windows(mus1,
                                            mus2,
                                            size=size,
                                            min_count=min_count):
        values.loc[center] = method(win1, win2)
    return values
