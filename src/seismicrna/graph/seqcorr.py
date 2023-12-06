import os
from logging import getLogger
from typing import Callable

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairGraphWriter, SeqPairOneAxisGraph
from .traces import iter_seq_line_traces
from ..core.arg import (METHOD_DETERM,
                        METHOD_PEARSON,
                        METHOD_SPEARMAN,
                        opt_correl,
                        opt_window,
                        opt_winmin)
from ..core.seq import get_shared_index, get_windows

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


def _get_correl(func: Callable, square: bool = False):
    """ Correlation method for comparing datasets. """

    def correl(data1: pd.Series, data2: pd.Series):
        # Compute the correlation between the two series over the
        # window, ignoring missing data.
        statistic = func(data1, data2, nan_policy="omit").statistic
        return statistic * statistic if square else statistic

    return correl


def _get_method(method: str):
    """ Method for comparing datasets. """
    if method == METHOD_SOMETHING:
        pass
    if method == METHOD_DETERM:
        # Use the coefficient of determination.
        from scipy.stats import pearsonr
        return _get_correl(pearsonr, square=True)
    if method == METHOD_PEARSON:
        # Use Pearson's correlation coefficient.
        from scipy.stats import pearsonr
        return _get_correl(pearsonr)
    if method == METHOD_SPEARMAN:
        # Use Spearman's correlation coefficient.
        from scipy.stats import spearmanr
        return _get_correl(spearmanr)
    raise ValueError(f"Invalid comparison method: {repr(method)}")


class SeqCorrGraph(SeqPairOneAxisGraph):

    def __init__(self, *, window: int, winmin: int, method: str, **kwargs):
        # Import scipy here instead of at the top of this module because
        # its import is slow enough to impact global startup time.
        super().__init__(**kwargs)
        self._window = window
        self._winmin = winmin
        self._method = _get_method(method)

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
        def get_rolling(vals1: pd.Series, vals2: pd.Series):
            """ Compute the rolling comparison between the Series. """
            # Initialize an empty Series for the rolling comparison.
            rolling = pd.Series(np.nan, index=get_shared_index([vals1.index,
                                                                vals2.index]))
            # Compare each window.
            for center, (win1, win2) in get_windows(vals1,
                                                    vals2,
                                                    size=self._window,
                                                    min_count=self._winmin):
                rolling.loc[center] = self._method(win1, win2)
            return rolling

        return get_rolling

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
