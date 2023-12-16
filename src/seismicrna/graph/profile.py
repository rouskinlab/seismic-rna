import os
from abc import ABC
from logging import getLogger

from click import command
from plotly import graph_objects as go

from .base import arrange_table
from .color import ColorMapGraph, RelColorMap, SeqColorMap
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .traces import iter_seq_base_bar_traces, iter_seqbar_stack_traces
from ..core.seq import POS_NAME

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class ProfileGraph(OneTableGraph, ColorMapGraph, ABC):
    """ Graph of a mutational profile. """

    @property
    def x_title(self):
        return POS_NAME

    @property
    def y_title(self):
        return self.data_kind.capitalize()


class SingleRelProfileGraph(ProfileGraph, ABC):
    """ Bar graph where each bar shows one relationship of the base. """

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_seq_base_bar_traces(values, self.cmap):
                yield (row, 1), trace


class MultiRelProfileGraph(ProfileGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    relationships for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        # Stack the bars at each position.
        figure.update_layout(barmode="stack")

    def get_traces(self):
        for row, index in enumerate(self.row_index, start=1):
            for trace in iter_seqbar_stack_traces(self.data.loc[:, index],
                                                  self.cmap):
                yield (row, 1), trace


class ProfileWriter(OneTableWriter):

    @classmethod
    def graph_type(cls):
        return ProfileGraph

    def iter(self,
             rels_sets: tuple[str, ...],
             use_ratio: bool,
             quantile: float,
             arrange: str):
        for cparams in arrange_table(self.table, arrange):
            for rels in rels_sets:
                graph_type = (SingleRelProfileGraph
                              if len(rels) == 1
                              else MultiRelProfileGraph)
                yield graph_type(table=self.table,
                                 rels=rels,
                                 use_ratio=use_ratio,
                                 quantile=quantile,
                                 **cparams)


class ProfileRunner(OneTableRunner):

    @classmethod
    def writer_type(cls):
        return ProfileWriter


@command(COMMAND, params=ProfileRunner.params())
def cli(*args, **kwargs):
    """ Create bar graphs of positions in a sequence. """
    return ProfileRunner.run(*args, **kwargs)

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
