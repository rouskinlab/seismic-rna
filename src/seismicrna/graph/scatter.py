import os
from itertools import product
from logging import getLogger

from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairTwoAxisGraph, SeqPairGraphWriter
from .traces import iter_seq_base_scatter_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class SeqScatterGraphRunner(SeqPairGraphRunner):

    @classmethod
    def writer_type(cls):
        return SeqScatterGraphWriter


@command(COMMAND, params=SeqScatterGraphRunner.params)
def cli(*args, **kwargs):
    """ Create scatter plots between pairs of samples at each position
    in a sequence. """
    return SeqScatterGraphRunner.run(*args, **kwargs)


class SeqScatterGraph(SeqPairTwoAxisGraph):

    @classmethod
    def graph_type(cls):
        return COMMAND

    def get_traces(self):
        for (col, (_, vals1)), (row, (_, vals2)) in product(
                enumerate(self.data1.items(), start=1),
                enumerate(self.data2.items(), start=1)
        ):
            for trace in iter_seq_base_scatter_traces(vals1, vals2, self.cmap):
                yield (row, col), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_xaxes(gridcolor="#d0d0d0")
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqScatterGraphWriter(SeqPairGraphWriter):

    @property
    def graph_type(self):
        return SeqScatterGraph

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
