import os
from logging import getLogger

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairGraphWriter, SeqPairOneAxisGraph
from .traces import iter_seq_line_traces
from ..core.arg import (CORREL_DETERM,
                        CORREL_PEARSON,
                        CORREL_SPEARMAN,
                        opt_correl,
                        opt_window,
                        opt_winmin)
from ..core.seq import get_shared_index, get_windows

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class SeqCorrGraph(SeqPairOneAxisGraph):

    def __init__(self, *, window: int, winmin: int, correl: str, **kwargs):
        # Import scipy here instead of at the top of this module because
        # its import is slow enough to impact global startup time.
        from scipy.stats import pearsonr, spearmanr
        super().__init__(**kwargs)
        self._window = window
        self._winmin = winmin
        if correl == CORREL_DETERM:
            self._correl = pearsonr
            self._square = True
        elif correl == CORREL_PEARSON:
            self._correl = pearsonr
            self._square = False
        elif correl == CORREL_SPEARMAN:
            self._correl = spearmanr
            self._square = False
        else:
            raise ValueError(f"Invalid correlation coefficient: {repr(correl)}")

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
            # Initialize an empty Series for the rolling correlation.
            roll_corr = pd.Series(np.nan, index=get_shared_index([vals1.index,
                                                                  vals2.index]))
            # Compute the correlation for each rolling window.
            for center, (win1, win2) in get_windows(vals1,
                                                    vals2,
                                                    size=self._window,
                                                    min_count=self._winmin):
                # Compute the correlation between the two series over
                # the window, ignoring missing data.
                corr = self._correl(win1, win2, nan_policy="omit").statistic
                if self._square:
                    # Square the correlation.
                    corr = corr * corr
                roll_corr.loc[center] = corr
            return roll_corr

        return get_rolling_corr

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqCorrGraphWriter(SeqPairGraphWriter):

    @classmethod
    def graph_type(cls):
        return SeqCorrGraph


class SeqCorrGraphRunner(SeqPairGraphRunner):

    @classmethod
    def var_params(cls):
        return [opt_correl, opt_window, opt_winmin]

    @classmethod
    def writer_type(cls):
        return SeqCorrGraphWriter


@command(COMMAND, params=SeqCorrGraphRunner.params())
def cli(*args, **kwargs):
    """ Create line graphs of rolling correlations between datasets. """
    return SeqCorrGraphRunner.run(*args, **kwargs)

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
