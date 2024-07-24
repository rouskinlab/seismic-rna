import os
from abc import ABC
from functools import cached_property
from logging import getLogger

from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .color import ColorMapGraph, RelColorMap, SeqColorMap
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import MultiRelsGraph, OneRelGraph
from .trace import iter_seq_base_bar_traces, iter_seqbar_stack_traces
from ..core.header import parse_header
from ..core.seq import POS_NAME

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class ProfileGraph(OneTableGraph, ColorMapGraph, ABC):
    """ Bar graph of a mutational profile for one table. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @property
    def x_title(self):
        return POS_NAME

    @property
    def y_title(self):
        return self.data_kind.capitalize()

    @cached_property
    def data(self):
        return self._fetch_data(self.table,
                                k=self.k,
                                clust=self.clust)

    @cached_property
    def data_header(self):
        """ Header of the selected data (not of the entire table). """
        return parse_header(self.data.columns)


class OneRelProfileGraph(OneRelGraph, ProfileGraph):
    """ Bar graph with one relationship per position. """

    @classmethod
    def what(cls):
        return "Profile"

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_seq_base_bar_traces(values, self.cmap):
                yield (row, 1), trace


class MultiRelsProfileGraph(MultiRelsGraph, ProfileGraph):
    """ Stacked bar graph with multiple relationships per position. """

    @classmethod
    def what(cls):
        return "Typed profile"

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        # Stack the bars at each position.
        figure.update_layout(barmode="stack")

    def get_traces(self):
        for row, index in zip(range(1, self.nrows + 1),
                              self.data_header.iter_clust_indexes(),
                              strict=True):
            for trace in iter_seqbar_stack_traces(self.data.loc[:, index],
                                                  self.cmap):
                yield (row, 1), trace


class ProfileWriter(OneTableWriter, PosGraphWriter):

    def get_graph(self, rels_group: str, **kwargs):
        return (OneRelProfileGraph(table=self.table,
                                   rel=rels_group,
                                   **kwargs)
                if len(rels_group) == 1
                else MultiRelsProfileGraph(table=self.table,
                                           rels=rels_group,
                                           **kwargs))


class ProfileRunner(OneTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return ProfileWriter


@command(COMMAND, params=ProfileRunner.params())
def cli(*args, **kwargs):
    """ Bar graph of relationships(s) per position. """
    return ProfileRunner.run(*args, **kwargs)

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
