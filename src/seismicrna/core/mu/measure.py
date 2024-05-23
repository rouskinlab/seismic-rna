"""
Calculate trends in mutation rates.
"""

import numpy as np
import pandas as pd
from numba import jit

from .dim import count_pos
from .frame import auto_reframe
from .nan import auto_remove_nan


@jit()
def _calc_sum_abs_diff(x: np.ndarray):
    """ Sum the absolute difference along axis 0. """
    n = x.shape[0]
    sum_abs_diff = np.zeros(x.shape[1:])
    for i in range(n):
        for j in range(i + 1, n):
            sum_abs_diff += np.abs(x[i] - x[j])
    return sum_abs_diff


@auto_remove_nan
@auto_reframe
def calc_gini(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the Gini coefficient of mutation rates, ignoring NaNs.

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
    if (npos := count_pos(mus)) == 0:
        # If there are no positions, then return an all-NaN array with
        # the same dimensions as the input but without axis 0.
        return np.full(mus.shape[1:], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _calc_sum_abs_diff(mus) / (npos * npos * mus.mean(axis=0))


@auto_remove_nan
@auto_reframe
def calc_signal_noise(mus: np.ndarray | pd.Series | pd.DataFrame,
                      is_signal: np.ndarray | pd.Series):
    """ Calculate the signal-to-noise ratio of mutation rates.

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
    signal = mus[is_signal]
    noise = mus[~is_signal]
    if count_pos(signal) == 0 or count_pos(noise) == 0:
        # If there is not at least one signal and at least one noise,
        # then return an all-NaN array with the same dimensions as the
        # input but without axis 0.
        return np.full(mus.shape[1:], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        return signal.mean(axis=0) / noise.mean(axis=0)

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
