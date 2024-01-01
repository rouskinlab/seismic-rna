import os
from functools import cached_property
from itertools import product
from logging import getLogger

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .color import ColorMapGraph, SeqColorMap
from .trace import iter_seq_base_scatter_traces
from .twotable import SAMPLE_NAME, TwoTableGraph, TwoTableRunner, TwoTableWriter

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class ScatterPlotGraph(TwoTableGraph, ColorMapGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Scatter plot"

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @property
    def x_title(self):
        return self.sample1

    @property
    def y_title(self):
        return self.sample2

    @cached_property
    def data(self):
        # Join data tables 1 and 2 horizontally.
        data = pd.concat([self.data1, self.data2], axis=1, join="inner")
        # Add the sample names as the first level of the columns.
        samples = np.hstack([np.repeat(self.sample1, self.data1.columns.size),
                             np.repeat(self.sample2, self.data2.columns.size)])
        names = [SAMPLE_NAME] + list(data.columns.names)
        data.columns = pd.MultiIndex.from_arrays(
            [(samples if name == SAMPLE_NAME
              else data.columns.get_level_values(name).values)
             for name in names],
            names=names
        )
        return data

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


class ScatterPlotWriter(TwoTableWriter, PosGraphWriter):

    @classmethod
    def get_graph_type(cls):
        return ScatterPlotGraph


class ScatterPlotRunner(TwoTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return ScatterPlotWriter


@command(COMMAND, params=ScatterPlotRunner.params())
def cli(*args, **kwargs):
    """ Scatter plot comparing two profiles. """
    return ScatterPlotRunner.run(*args, **kwargs)

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
