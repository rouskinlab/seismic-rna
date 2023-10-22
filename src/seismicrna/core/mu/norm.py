from logging import getLogger
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov

from ..seq import Section

logger = getLogger(__name__)



def get_mu_quantile(mus: np.ndarray | pd.Series, quantile: float) -> float:
    """ Compute the mutation rate at a quantile. """
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
            value = np.nanquantile(mus, quantile)
    if np.isnan(value):
        # If a NaN value was returned for either of the above reasons,
        # then issue a warning.
        logger.warning(f"Got NaN quantile {quantile} of {mus}")
    return value


def normalize(mus: np.ndarray | pd.Series, quantile: float):
    """ Normalize the mutation rates to a quantile. """
    if quantile == 0.:
        # Do not normalize the mutation rates if quantile == 0.
        return mus.copy()
    return mus / get_mu_quantile(mus, quantile)


def winsorize(mus: np.ndarray | pd.Series, quantile: float):
    """ Normalize and winsorize the mutation rates to a quantile. """
    winsorized = np.clip(normalize(mus, quantile), 0., 1.)
    # Return the same data type as was given for mus.
    if isinstance(mus, pd.Series):
        return pd.Series(winsorized, index=mus.index)
    return winsorized
