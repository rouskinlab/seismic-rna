from __future__ import annotations

from .dim import count_pos
from .frame import auto_reframe
from .nan import auto_remove_nan

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd


@auto_remove_nan
@auto_reframe
def calc_gini(mus: np.ndarray | pd.Series | pd.DataFrame):
    """Calculate the Gini coefficient of mutation rates, ignoring NaNs.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    float | numpy.ndarray | pandas.Series
        Value of the Gini coefficient.
    """
    import numpy as np

    if (npos := count_pos(mus)) == 0:
        # If there are no positions, then return an all-NaN array with
        # the same dimensions as the input but without axis 0.
        return np.full(mus.shape[1:], np.nan)
    # Imported here (not at module level) so importing this module does not
    # import numba via measure_jit.
    from .measure_jit import calc_sum_abs_diff

    with np.errstate(divide="ignore", invalid="ignore"):
        return calc_sum_abs_diff(mus) / (npos * npos * mus.mean(axis=0))


@auto_remove_nan
@auto_reframe
def calc_signal_noise(
    mus: np.ndarray | pd.Series | pd.DataFrame, is_signal: np.ndarray | pd.Series
):
    """Calculate the signal-to-noise ratio of mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a DataFrame.
    is_signal: np.ndarray | pd.Series
        Whether to count each position as signal.

    Returns
    -------
    float | numpy.ndarray | pandas.Series
        Signal-to-noise ratio.
    """
    import numpy as np

    signal = mus[is_signal]
    noise = mus[~is_signal]
    if count_pos(signal) == 0 or count_pos(noise) == 0:
        # If there is not at least one signal and at least one noise,
        # then return an all-NaN array with the same dimensions as the
        # input but without axis 0.
        return np.full(mus.shape[1:], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return signal.mean(axis=0) / noise.mean(axis=0)
