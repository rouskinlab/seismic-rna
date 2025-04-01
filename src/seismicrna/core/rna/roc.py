import numpy as np
import pandas as pd

from ..array import get_length
from ..seq import POS_NAME, iter_windows
from ..validate import require_isinstance, require_equal


def _sort_paired(profile: pd.Series, *structs: pd.Series | None):
    """ Ensure the structure and mutation data have the same indexes,
    filter out positions without mutation data, and sort the structure
    data in ascending order by mutation rate.

    Parameters
    ----------
    profile: pandas.Series
        Mutational profile with one index per position; each element is
        the mutation rate at the position.
    *structs: pandas.Series | None
        Boolean series with one index per position; each element is the
        pairing status of a base at that position. Ignore values that
        are None.

    Returns
    -------
    tuple[pandas.Series, ...]
        Values in `paired` sorted by values in `profile`.
    """
    assert isinstance(profile, pd.Series)
    # Use only positions with non-missing mutation data.
    profile_not_nan = profile.loc[np.logical_not(np.isnan(profile.values))]
    # Sort the positions in ascending order of the mutation data.
    sorted_indexes = profile_not_nan.sort_values().index.get_level_values(
        POS_NAME
    )
    sorted_structs = list()
    for struct in structs:
        if struct is not None:
            assert isinstance(struct, pd.Series)
            assert struct.index.equals(profile.index)
            sorted_structs.append(struct.loc[sorted_indexes])
    return tuple(sorted_structs)


def _compute_fpr_tpr(sorted_paired: pd.Series,
                     sorted_unpaired: pd.Series | None = None):
    """ Compute the receiver operating characteristic (ROC) curve to
    indicate how well chemical reactivities agree with a structure.

    Parameters
    ----------
    sorted_paired: pandas.Series
        Whether each position is paired, sorted in ascending order of
        mutation rate.
    sorted_unpaired: pandas.Series | None
        Whether each position is unpaired, sorted in ascending order of
        mutation rate.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        FPR and TPR axes, respectively.
    """
    assert isinstance(sorted_paired, pd.Series)
    assert sorted_paired.dtype == bool
    if sorted_unpaired is not None:
        assert isinstance(sorted_unpaired, pd.Series)
        assert sorted_unpaired.dtype == bool
        # The paired and unpaired series must have the same positions,
        # and no position can be both paired and unpaired, although a
        # position can be considered neither paired nor unpaired.
        assert sorted_unpaired.index.equals(sorted_paired.index)
        assert np.all(~(sorted_paired & sorted_unpaired))
    else:
        # Assume every position is either paired or unpaired.
        sorted_unpaired = ~sorted_paired
    # Count the positions.
    n_total = sorted_paired.size
    # Compute the cumulative number of paired and unpaired bases at each
    # threshold.
    paired_cumsum = np.hstack([[0], np.cumsum(sorted_paired)])
    unpaired_cumsum = np.hstack([[0], np.cumsum(sorted_unpaired)])
    # Count the paired and unpaired positions.
    n_paired = int(paired_cumsum[-1])
    n_unpaired = int(unpaired_cumsum[-1])
    # The sum of paired and unpaired positions equal to the total number
    # of positions, or be less if some positions are considered neither
    # paired nor unpaired.
    assert n_paired + n_unpaired <= n_total
    # Get the false positive rate: (paired and reactive) / paired.
    fpr = (1. - paired_cumsum / n_paired
           if n_paired > 0
           else np.full(n_total + 1, np.nan))
    # Get the true positive rate: (unpaired and reactive) / unpaired.
    tpr = (1. - unpaired_cumsum / n_unpaired
           if n_unpaired > 0
           else np.full(n_total + 1, np.nan))
    # Traditionally, false positive rate is plotted on the x-axis and
    # true positive rate on the y-axis of an ROC curve.
    return fpr, tpr


def compute_roc_curve(profile: pd.Series,
                      is_paired: pd.Series,
                      is_unpaired: pd.Series | None = None):
    """ Compute the receiver operating characteristic (ROC) curve to
    indicate how well mutation data agree with a structure.

    Parameters
    ----------
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.
    is_paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    is_unpaired: pandas.Series | None
        Boolean series with one index per position, where each value is
        True if the base at the position is unpaired, otherwise False.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        FPR and TPR axes, respectively, of the ROC curve.
    """
    return _compute_fpr_tpr(*_sort_paired(profile, is_paired, is_unpaired))


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
    require_isinstance("fpr", fpr, np.ndarray)
    require_isinstance("tpr", tpr, np.ndarray)
    n = get_length(fpr)
    require_equal("get_length(fpr)",
                  n,
                  get_length(tpr),
                  "get_length(tpr)")
    return float(-np.vdot(np.diff(fpr), tpr[1:])) if n > 1 else np.nan


def compute_auc_roc(profile: pd.Series, 
                    is_paired: pd.Series,
                    is_unpaired: pd.Series | None = None):
    """ Compute the receiver operating characteristic (ROC) and the area
    under the curve (AUC) to indicate how well mutation data agree with
    a structure.

    Parameters
    ----------
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.
    is_paired: pandas.Series
        Boolean series with one index per position, where each value is
        True if the base at the position is paired, otherwise False.
    is_unpaired: pandas.Series | None
        Boolean series with one index per position, where each value is
        True if the base at the position is unpaired, otherwise False.

    Returns
    -------
    float
        AUC-ROC
    """
    return compute_auc(*compute_roc_curve(profile, is_paired, is_unpaired))


def compute_rolling_auc(profile: pd.Series,
                        *structs: pd.Series,
                        size: int,
                        min_data: int = 2):
    """ Compute the area under the curve (AUC) of the receiver operating
    characteristic (ROC) at each position using a sliding window.

    Parameters
    ----------
    profile: pandas.Series
        Mutational profile with one index per position, where each value
        is the mutation rate at the position.
    *structs: pandas.Series
        Boolean series with one index per position; each element is the
        pairing status of a base at that position. Ignore values that
        are None.
    size: int
        Size of the window.
    min_data: int = 2
        Minimum number of data in a window to use it, otherwise NaN.

    Returns
    -------
    pandas.Series
        AUC-ROC at each position.
    """
    # Initialize an empty series to hold the AUC-ROC values.
    aucrocs = pd.Series(np.nan, index=profile.index)
    # Compute the AUC-ROC for each sliding window.
    for center, window in iter_windows(profile,
                                       *structs,
                                       size=size,
                                       min_count=min_data,
                                       include_nan=False):
        aucrocs.loc[center] = compute_auc_roc(*window)
    return aucrocs
