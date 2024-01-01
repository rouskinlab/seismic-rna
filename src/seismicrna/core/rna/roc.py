import numpy as np
import pandas as pd

from ..seq import POS_NAME, iter_windows


def _sort_paired(paired: pd.Series, profile: pd.Series):
    """ Ensure the structure and mutation data have the same indexes,
    filter out positions without mutation data, and sort the structure
    data in ascending order by mutation rate.

    Parameters
    ----------
    paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.

    Returns
    -------
    pandas.Series
        Values in `paired` sorted by values in `profile`.
    """
    if not paired.index.equals(profile.index):
        raise ValueError(f"Indexes of paired bases and mutational profile "
                         f"differ: {paired.index} ≠ {profile.index}")
    # Use only positions with non-missing mutation data.
    profile_not_nan = profile.loc[np.logical_not(np.isnan(profile.values))]
    # Sort the positions in ascending order of the mutation data.
    sorted_indexes = profile_not_nan.sort_values().index.get_level_values(POS_NAME)
    return paired.loc[sorted_indexes].astype(bool, copy=False)


def _compute_fpr_tpr(paired: pd.Series):
    """ Compute the receiver operating characteristic (ROC) curve to
    indicate how well chemical reactivities agree with a structure.

    Parameters
    ----------
    paired: pandas.Series
        Boolean series with one index per position, sorted in ascending
        order of the mutation data.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        FPR and TPR axes, respectively.
    """
    # Count the positions.
    n_total = paired.size
    # Compute the cumulative number of paired bases at each threshold.
    paired_cumsum = np.hstack([[0], np.cumsum(paired)])
    # Count the paired and unpaired positions.
    n_paired = int(paired_cumsum[-1])
    n_unpaired = n_total - n_paired
    # Get the false positive rate: (paired and reactive) / paired.
    fpr = (1. - paired_cumsum / n_paired
           if n_paired > 0
           else np.full(n_total + 1, np.nan))
    # Get the true positive rate: (unpaired and reactive) / unpaired.
    tpr = (1. - (np.arange(n_total + 1) - paired_cumsum) / n_unpaired
           if n_unpaired > 0
           else np.full(n_total + 1, np.nan))
    # Traditionally, false positive rate is plotted on the x-axis and
    # true positive rate on the y-axis of an ROC curve.
    return fpr, tpr


def compute_roc_curve(paired: pd.Series, profile: pd.Series):
    """ Compute the receiver operating characteristic (ROC) curve to
    indicate how well mutation data agree with a structure.

    Parameters
    ----------
    paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        FPR and TPR axes, respectively, of the ROC curve.
    """
    return _compute_fpr_tpr(_sort_paired(paired, profile))


def compute_auc(fpr: np.ndarray, tpr: np.ndarray):
    """ Compute the area under the curve (AUC) of the receiver operating
    characteristic (ROC).

    Parameters
    ----------
    fpr: numpy.ndarray
        False positive rate (FPR) of the ROC curve.
    tpr: numpy.ndarray
        True positive rate (TPR) of the ROC curve.

    Returns
    -------
    float
        AUC-ROC
    """
    return float(-np.vdot(np.diff(fpr), tpr[1:]))


def compute_auc_roc(paired: pd.Series, profile: pd.Series):
    """ Compute the receiver operating characteristic (ROC) and the area
    under the curve (AUC) to indicate how well mutation data agree with
    a structure.

    Parameters
    ----------
    paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.

    Returns
    -------
    float
        AUC-ROC
    """
    return compute_auc(*compute_roc_curve(paired, profile))


def compute_rolling_auc(paired: pd.Series,
                        profile: pd.Series,
                        size: int,
                        min_data: int = 2):
    """ Compute the area under the curve (AUC) of the receiver operating
    characteristic (ROC) at each position using a sliding window.

    Parameters
    ----------
    paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.
    size: int
        Size of the window.
    min_data: int = 2
        Minimum number of data in a window to use it (otherwise NaN).

    Returns
    -------
    pandas.Series
        AUC-ROC at each position.
    """
    # Initialize an empty series to hold the AUC-ROC values.
    aucrocs = pd.Series(np.nan, index=paired.index)
    # Compute the AUC-ROC for each sliding window.
    for center, window in iter_windows(paired,
                                       profile,
                                       size=size,
                                       min_count=min_data,
                                       include_nan=False):
        aucrocs.loc[center] = compute_auc_roc(*window)
    return aucrocs

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
