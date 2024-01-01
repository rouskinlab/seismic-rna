import os
from logging import getLogger

import numpy as np
from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .color import ColorMapGraph, SeqColorMap
from .trace import iter_seq_base_bar_traces
from .twotable import TwoTableRunner, TwoTableWriter, TwoTableMergedGraph
from ..core.seq import POS_NAME

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class DeltaProfileGraph(TwoTableMergedGraph, ColorMapGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Delta profile graph"

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @classmethod
    def _trace_function(cls):
        return iter_seq_base_bar_traces

    @property
    def x_title(self) -> str:
        return POS_NAME

    @property
    def _trace_kwargs(self):
        return super()._trace_kwargs | dict(cmap=self.cmap)

    @property
    def y_title(self):
        return f"Difference between {self.data_kind}s 1 and 2"

    @property
    def _merge_data(self):
        """ Compute the difference between the profiles. """
        return np.subtract

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class DeltaProfileWriter(TwoTableWriter, PosGraphWriter):

    @classmethod
    def get_graph_type(cls):
        return DeltaProfileGraph


class DeltaProfileRunner(TwoTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return DeltaProfileWriter


@command(COMMAND, params=DeltaProfileRunner.params())
def cli(*args, **kwargs):
    """ Bar graph of differences between two profiles per position. """
    return DeltaProfileRunner.run(*args, **kwargs)

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
