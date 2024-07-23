"""
Scale mutation rates.
"""

import numpy as np
import pandas as pd

from .dim import count_pos
from .frame import auto_reframe
from .nan import auto_remove_nan


@auto_remove_nan
@auto_reframe
def calc_quantile(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Calculate the mutation rate at a quantile, ignoring NaNs.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile to return from the mutation rates; must be in [0, 1].

    Returns
    -------
    float | numpy.ndarray | pandas.Series
        Value of the quantile from the mutation rates.
    """
    if not isinstance(quantile, float):
        # Although numpy.quantile supports array-like values for the
        # quantile argument, get_quantile does not because the result
        # would have one or more extra dimensions.
        raise TypeError("Expected quantile to be float, "
                        f"but got {type(quantile).__name__}")
    if count_pos(mus) == 0:
        # If there are no positions, then return an all-NaN array with
        # the same dimensions as the input but without axis 0, instead
        # of raising an error, which np.quantile would do.
        return np.full(mus.shape[1:], np.nan)
    # Return the quantile by reducing along axis 0.
    return np.quantile(mus, quantile, axis=0)


def normalize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize the mutation rates to a quantile, so that the value of
    the quantile is scaled to 1 and all other mutation rates are scaled
    by the same factor. If quantile is 0, then do not normalize.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile for normalizing the mutation rates; must be in [0, 1].

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Normalized mutation rates.
    """
    return mus / calc_quantile(mus, quantile) if quantile > 0. else mus


def normalize_max(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Normalize the mutation rates so their maximum becomes 1. """
    return normalize(mus, 1.0)


def normalize_med(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Normalize the mutation rates so their median becomes 1. """
    return normalize(mus, 0.5)


@auto_reframe
def winsorize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize and winsorize the mutation rates to a quantile so that
    all mutation rates greater than or equal to the mutation rate at the
    quantile are set to 1, and lesser mutation rates are normalized.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.
    quantile: float
        Quantile for normalizing the mutation rates; must be in [0, 1].

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Normalized and winsorized mutation rates.
    """
    return np.clip(normalize(mus, quantile), 0., 1.)


@auto_remove_nan
def calc_rms(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square mutation rate, ignoring NaNs.

    Parameters
    ----------
    mus: np.ndarray | pd.Series | pd.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    float | numpy.ndarray | pandas.Series
        Root-mean-square mutation rate.
    """
    return np.sqrt(np.mean(np.square(mus), axis=0))


def standardize(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Standardize mutation rates so that the root-mean-square mutation
    rate equals 1. Note that approximately half of the standardized
    mutation rates will be greater than 1.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Standardized mutation rates.
    """
    return mus / calc_rms(mus)


@auto_remove_nan
@auto_reframe
def calc_ranks(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Rank the mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Ranks of the mutation rates.
    """
    # Initialize all ranks to -1.
    ranks = np.full_like(mus, -1, dtype=int)
    # Make coordinate arrays to index every element of ranks in order.
    coords = np.unravel_index(np.arange(ranks.size), ranks.shape)
    # Replace the coordinate array of the first axis with the indexes
    # in order of their ranks.
    coords[0][:] = np.argsort(mus, axis=0).reshape(coords[0].shape)
    # Fill the ranks.
    ranks[coords] = np.ravel(
        np.broadcast_to(np.expand_dims(np.arange(mus.shape[0]),
                                       axis=tuple(range(1, mus.ndim))),
                        ranks.shape)
    )
    # Confirm all ranks were filled.
    if np.any(ranks < 0):
        # This error should be impossible, but checking just in case.
        raise ValueError(f"Ranks cannot be negative, but got {ranks}")
    return ranks

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
