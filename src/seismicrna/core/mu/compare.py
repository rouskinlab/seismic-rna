from typing import Callable

import numpy as np
import pandas as pd

from .nan import auto_removes_nan
from .frame import find_highest_type
from .scale import calc_rms, calc_ranks, normalize_max
from ..arg import KEY_DETERM, KEY_PEARSON, KEY_NRMSD, KEY_SPEARMAN
from ..seq import get_shared_index, iter_windows

DEFAULT_CLIP_LOG_ODDS = 1.e-3


def _calc_diff_log_odds(mus1: float | np.ndarray | pd.Series | pd.DataFrame,
                        mus2: float | np.ndarray | pd.Series | pd.DataFrame,
                        p_min: float,
                        p_max: float):
    assert 0. <= p_min <= p_max <= 1.
    mus1 = np.clip(mus1, p_min, p_max)
    mus2 = np.clip(mus2, p_min, p_max)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log((mus1 * (1. - mus2)) / ((1. - mus1) * mus2))


def calc_diff_log_odds(mus1: float | np.ndarray | pd.Series | pd.DataFrame,
                       mus2: float | np.ndarray | pd.Series | pd.DataFrame,
                       p_min: float = DEFAULT_CLIP_LOG_ODDS,
                       p_max: float = (1. - DEFAULT_CLIP_LOG_ODDS)):
    """ Calculate the difference in log odds between mus1 and mus2.
    Assume that mus1 and mus2 are on the same scale (e.g. two clusters
    from the same sample), so perform no scaling or normalization.
    If mus1 and mus2 are equal, they will always return 0, even if they
    are both 0 or both 1.

    Parameters
    ----------
    mus1: np.ndarray | pd.Series | pd.DataFrame
        First group of mutation rates; can contain multiple sets as the
        columns of a multidimensional array or DataFrame.
    mus2: np.ndarray | pd.Series | pd.DataFrame
        Second group of mutation rates; can contain multiple sets as the
        columns of a multidimensional array or DataFrame.
    p_min: float
        Ensure that all probabilities are at least this value to prevent
        dividing by or taking the log of 0.
    p_max: float
        Ensure that all probabilities are at most this value to prevent
        dividing by or taking the log of 0.

    Returns
    -------
    float | np.ndarray | pd.Series | pd.DataFrame
        Difference in log odds:
        log(mus1 / (1 - mus1)) - log(mus2 / (1 - mus2))
    """
    if not 0. <= p_min <= p_max <= 1.:
        raise ValueError(f"Must have 0 ≤ p_min ≤ p_max ≤ 1, "
                         f"but got p_min={p_min} and p_max={p_max}")
    diff_log_odds = _calc_diff_log_odds(mus1, mus2, p_min, p_max)
    if 0. < p_min <= p_max < 1.:
        # If the probabilities are clipped on both sides, then two 0s
        # will always give a log odds difference of 0, as will two 1s.
        return diff_log_odds
    # If p_min == 0, then two 0s will give a log odds difference of NaN;
    # if p_max == 1, then two 1s will give a log odds difference of NaN;
    # so equal input values must explicitly return 0.
    indexes = dict()
    result_type = find_highest_type(mus1, mus2, creatable=True)
    if result_type is pd.DataFrame or result_type is pd.Series:
        if not mus1.index.equals(mus2.index):
            raise ValueError("Indexes of mus1 and mus2 are different")
        indexes["index"] = mus1.index
        if result_type is pd.DataFrame:
            if not mus1.columns.equals(mus2.columns):
                raise ValueError("Columns of mus1 and mus2 are different")
            indexes["columns"] = mus1.columns
    return result_type(np.where(mus1 == mus2, 0., diff_log_odds), **indexes)


@auto_removes_nan
def calc_sum_abs_diff_log_odds(mus1: np.ndarray | pd.Series | pd.DataFrame,
                               mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the sum of absolute differences in log odds between
    mus1 and mus2. See calc_diff_log_odds for details.

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
        Sum of absolute differences in log odds
    """
    return np.sum(np.abs(calc_diff_log_odds(mus1, mus2)), axis=0)


@auto_removes_nan
def calc_raw_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
                  mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square difference (RMSD) of two groups
    of mutation rates, ignoring NaNs. Assume that mus1 and mus2 are on
    the same scale (e.g. two clusters from the same sample), so perform
    no scaling or normalization.

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
        Root-mean-square deviation (RMSD)
    """
    return np.sqrt(np.mean(np.square(mus1 - mus2), axis=0))


@auto_removes_nan
def calc_std_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
                  mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the standardized root-mean-square difference (RMSD)
    of two groups of mutation rates, ignoring NaNs. Assume that mus1
    and mus2 may be on different scales (e.g. different experiments),
    so scale each to the same root-mean-square value before calculating
    RMSD, then adjust the RMSD back to the original scale.

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
        Standardized root-mean-square deviation (SRMSD)
    """
    # Compute the root-mean-square mutation rate for each group.
    rms1 = calc_rms(mus1)
    rms2 = calc_rms(mus2)
    # Standardize the mutation rates so that the root-mean-square of
    # each group is 1, compute the root-mean-square difference, and
    # restore the original scale by multiplying by the geometric mean
    # of the root-mean-square mutation rates.
    return calc_raw_rmsd(mus1 / rms1, mus2 / rms2) * np.sqrt(rms1 * rms2)


@auto_removes_nan
def calc_norm_rmsd(mus1: np.ndarray | pd.Series | pd.DataFrame,
                   mus2: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the normalized root-mean-square difference (NRMSD)
    of two groups of mutation rates, ignoring NaNs. Like calc_std_rmsd,
    except both groups are initially scaled so that their maxima are 1,
    which makes it possible to compare the normalized RMSD between two
    datasets with high mutation rates to that betweeen two datasets with
    low mutation rates.

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
        Normalized root-mean-square deviation (NRMSD)
    """
    # Normalize the mutation rates so the maximum of each group is 1.
    return calc_std_rmsd(normalize_max(mus1), normalize_max(mus2))


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
        return calc_norm_rmsd, "Normalized Root-Mean-Square Deviation"
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
