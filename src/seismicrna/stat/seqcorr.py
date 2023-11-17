import os
from logging import getLogger

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairGraphWriter, SeqPairOneAxisGraph
from .traces import iter_seq_line_traces
from ..core.seq import POS_NAME

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class SeqCorrGraphRunner(SeqPairGraphRunner):

    @classmethod
    def writer_type(cls):
        return SeqCorrGraphWriter


@command(COMMAND, params=SeqPairGraphRunner.params)
def cli(*args, **kwargs):
    """ Create line graphs of rolling correlations between datasets. """
    return SeqCorrGraphRunner.run(*args, **kwargs)


class SeqCorrGraph(SeqPairOneAxisGraph):

    def __init__(self, *args,
                 margin: int = 20,
                 correl_min_n: int = 12,
                 correl_method: str = "spearman",
                 **kwargs):
        # Import scipy here instead of at the top of this module because
        # its import is slow enough to impact global startup time.
        from scipy.stats import pearsonr, spearmanr
        super().__init__(*args, **kwargs)
        self._margin = margin
        self._correl_min_n = correl_min_n
        if correl_method == "pearson":
            self._correl_method = pearsonr
        elif correl_method == "spearman":
            self._correl_method = spearmanr
        else:
            raise ValueError(f"Invalid correlation method: '{correl_method}'")

    @classmethod
    def graph_type(cls):
        return COMMAND

    @property
    def y_title(self):
        return f"Correlation of {self.quantity}-1 and {self.quantity}-2"

    @classmethod
    def _trace_function(cls):
        return iter_seq_line_traces

    @property
    def _merge_data(self):

        def get_rolling_corr(vals1: pd.Series, vals2: pd.Series):
            """ Compute the rolling correlation between the Series. """
            if not (index := vals1.index).equals(vals2.index):
                raise ValueError("Got different indexes for series 1 "
                                 f"({vals1.index}) and 2 ({vals2.index})")
            # Get the numeric positions in the sequence.
            pos = index.get_level_values(POS_NAME)
            if pos.has_duplicates:
                raise ValueError(f"Duplicate positions: {pos.duplicated()}")
            if not pos.is_monotonic_increasing:
                raise ValueError(f"Positions are not sorted: {pos}")
            # Initialize an empty Series for the rolling correlation.
            roll_corr = pd.Series(np.nan, index=index)
            # Compute the correlation for each rolling window.
            for win_pos in pos:
                # Determine the 5' and 3' ends of the window.
                end5 = win_pos - self._margin
                end3 = win_pos + self._margin
                # Get the both sets of values in the rolling window.
                win1 = vals1.loc[end5: end3]
                win2 = vals2.loc[end5: end3]
                # Count the positions at which both windows have values.
                not_isnan = np.logical_not(np.logical_or(np.isnan(win1),
                                                         np.isnan(win2)))
                if np.count_nonzero(not_isnan) < self._correl_min_n:
                    # If there are not enough positions where both series
                    # have values, then skip this position.
                    continue
                # Compute the correlation between the two series over
                # the window, ignoring missing data.
                win_corr = self._correl_method(win1, win2, nan_policy="omit")
                roll_corr.loc[win_pos] = win_corr.statistic
            return roll_corr

        return get_rolling_corr

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqCorrGraphWriter(SeqPairGraphWriter):

    @property
    def graph_type(self):
        return SeqCorrGraph

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
