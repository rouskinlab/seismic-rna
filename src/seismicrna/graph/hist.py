from abc import ABC

import numpy as np
import pandas as pd

from .base import GraphBase, GraphRunner
from ..core.arg import opt_hist_bins, opt_hist_margin


COUNT_NAME = "Count"
LOWER_NAME = "Lower"
UPPER_NAME = "Upper"


def get_edges_index(edges: np.ndarray, use_ratio: bool):
    """ Generate an index for the edges of histogram bins.

    Parameters
    ----------
    edges: numpy.ndarray
        Edges of histogram bins.
    use_ratio: bool
        Assume the edges represent ratios rather than counts.

    Returns
    -------
    pandas.Index | pandas.MultiIndex
        Index for the edges.
    """
    # Determine the lower and upper edge of each bin.
    lower = edges[:-1]
    upper = edges[1:]
    if use_ratio:
        # Make a MultiIndex of the lower and upper edge of each bin.
        return pd.MultiIndex.from_arrays([lower, upper],
                                         names=[LOWER_NAME, UPPER_NAME])
    # Make an Index of the count for each bin.
    if lower.size > 0:
        min_count = round(lower[0] + 0.5)
        counts = np.arange(min_count, min_count + lower.size)
        if not np.allclose(lower + 0.5, counts):
            raise ValueError("If the edges represent counts, then every lower "
                             "edge must be a half-integer and 1 more than the "
                             f"previous lower edge, but got {lower}")
    else:
        counts = np.arange(0)
    return pd.Index(counts, name=COUNT_NAME)


class HistogramGraph(GraphBase, ABC):
    """ Generic histogram. """

    def __init__(self, *, hist_bins: int, hist_margin: float, **kwargs):
        super().__init__(**kwargs)
        if hist_bins <= 0:
            raise ValueError(
                f"hist_bins must be a positive integer, but got {hist_bins}"
            )
        self.num_bins = hist_bins
        if hist_margin < 0.:
            raise ValueError(f"hist_margin must be ≥ 0, but got {hist_margin}")
        self.margin = hist_margin

    def get_bounds(self, data: pd.DataFrame):
        """ Get the lower and upper bounds of the histogram. """
        try:
            lo = float(np.nanmin(data))
        except ValueError:
            # If data contains no non-NaN values, then default to 0.
            lo = 0.
        else:
            # If the minimum is slightly greater than 0, then make it 0.
            if self.use_ratio and 0. < lo < self.margin:
                lo = 0.
        try:
            up = float(np.nanmax(data))
        except ValueError:
            # If data contains no non-NaN values, then default to 1.
            up = 1.
        else:
            # If the maximum is slightly less than 1, then make it 1.
            if self.use_ratio and 0. < 1. - up < self.margin:
                up = 1.
        return lo, up

    def get_edges(self, data: pd.DataFrame):
        """ Get the edges of the histogram bins. """
        # Make bins of floating-point values.
        lower, upper = self.get_bounds(data)
        if self.use_ratio:
            # Use the specified number of bins.
            num_bins = self.num_bins
        else:
            # Make one bin per integer, with half-integer bounds.
            lower = round(lower) - 0.5
            upper = round(upper) + 0.5
            num_bins = round(upper - lower)
        return np.linspace(lower, upper, num_bins + 1)


class HistogramRunner(GraphRunner, ABC):

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_hist_bins, opt_hist_margin]

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
