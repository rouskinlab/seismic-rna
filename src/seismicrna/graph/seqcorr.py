import os
from functools import partial
from logging import getLogger

from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairGraphWriter, SeqPairOneAxisGraph
from .traces import iter_seq_line_traces
from ..core.arg import opt_mucomp, opt_window, opt_winmin
from ..core.mu import compare_windows, get_comp_abbr

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class SeqCorrGraph(SeqPairOneAxisGraph):

    @classmethod
    def graph_type(cls):
        return COMMAND

    @classmethod
    def _trace_function(cls):
        return iter_seq_line_traces

    def __init__(self, *, mucomp: str, window: int, winmin: int, **kwargs):
        super().__init__(**kwargs)
        self._method = mucomp
        self._size = window
        self._min_count = winmin

    @property
    def y_title(self):
        return f"{get_comp_abbr(self._method)} of {self.quantity}s"

    @property
    def _merge_data(self):
        return partial(compare_windows,
                       method=self._method,
                       size=self._size,
                       min_count=self._min_count)

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
        return [opt_mucomp, opt_window, opt_winmin]

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
