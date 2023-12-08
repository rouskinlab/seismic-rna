"""
Calculations on single series of mutation rates.
"""

import warnings
from logging import getLogger

import numpy as np
import pandas as pd

from ..arrays import recontainerize

logger = getLogger(__name__)


def get_ranks(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Compute the ranks of the mutation rates.

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
    if mus.ndim == 0:
        raise ValueError(f"mus must have at least 1 dimension, but got {mus}")
    if np.any(np.isnan(mus)):
        raise ValueError(f"Cannot rank mutation rates with NaN values: {mus}")
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
    return recontainerize(ranks, mus)


def get_quantile(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Compute the mutation rate at a quantile.

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
    if mus.size == 0:
        # If there are no values, then return NaN instead of raising an
        # error, as np.nanquantile would.
        value = np.nan
    else:
        with warnings.catch_warnings():
            # Temporarily ignore warnings about all-NaN arrays.
            warnings.simplefilter("ignore")
            # Determine the quantile or, if the array is all NaN, then
            # set value to NaN.
            value = np.nanquantile(mus, quantile, axis=0)
    if np.any(np.isnan(np.atleast_1d(value))):
        # If a NaN value was returned for either of the above reasons,
        # then issue a warning.
        logger.warning(f"Got NaN quantile {quantile} of {mus}")
    return value


def normalize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize the mutation rates to a quantile.

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
    if quantile == 0.:
        # Do not normalize the mutation rates if quantile == 0.
        return mus.copy()
    return mus / get_quantile(mus, quantile)


def winsorize(mus: np.ndarray | pd.Series | pd.DataFrame, quantile: float):
    """ Normalize and winsorize the mutation rates to a quantile.

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
    winsorized = np.clip(normalize(mus, quantile), 0., 1.)
    # Return the same data type as was given for mus.
    if isinstance(mus, pd.DataFrame):
        return pd.DataFrame(winsorized,
                            index=mus.index,
                            columns=mus.columns,
                            copy=False)
    if isinstance(mus, pd.Series):
        return pd.Series(winsorized,
                         index=mus.index,
                         copy=False)
    return winsorized


def calc_rms(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the root-mean-square mutation rate.

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
    return np.sqrt(np.nanmean(mus * mus, axis=0))


def standardize(mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Standardize mutation rates so that the root-mean-square mutation
    rate equals 1.

    Parameters
    ----------
    mus: np.ndarray | pd.Series | pd.DataFrame
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    np.ndarray | pd.Series | pd.DataFrame
        Standardized mutation rates.
    """
    return mus / calc_rms(mus)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
