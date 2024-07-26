"""
Pairwise comparisons of mutation rates.
"""

from typing import Callable

import numpy as np
import pandas as pd

from .nan import auto_removes_nan
from .scale import calc_rms, calc_ranks, normalize_max
from ..arg import KEY_DETERM, KEY_PEARSON, KEY_NRMSD, KEY_SPEARMAN
from ..seq import get_shared_index, iter_windows


@auto_removes_nan
def calc_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
              mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square difference (RMSD) of two groups
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
    np.ndarray | pd.Series | pd.DataFrame
        Root-mean-square deviation (RMSD)
    """
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
def calc_nrmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
               mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the normalized root-mean-square difference (NRMSD)
    of two groups of mutation rates, ignoring NaNs.

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
        Normalized root-mean-square deviation (NRMSD)
    """
    # Normalize the mutation rates so the maximum of each group is 1.
    return calc_rmsd(normalize_max(mus1),
                     normalize_max(mus2))


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


def get_comp_method(key: str):
    """ Get a comparison method based on its key. """
    lowerkey = key.lower()
    if lowerkey == KEY_NRMSD:
        return calc_nrmsd, "Normalized Root-Mean-Square Deviation"
    if lowerkey == KEY_PEARSON:
        return calc_pearson, "Pearson Correlation Coefficient"
    if lowerkey == KEY_SPEARMAN:
        return calc_spearman, "Spearman Correlation Coefficient"
    if lowerkey == KEY_DETERM:
        return calc_coeff_determ, "Coefficient of Determination"
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
    func, _ = get_comp_method(key)
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
    _, name = get_comp_method(key)
    return name


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
    for center, (win1, win2) in iter_windows(mus1,
                                             mus2,
                                             size=size,
                                             min_count=min_count):
        values.loc[center] = method(win1, win2)
    return values

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
