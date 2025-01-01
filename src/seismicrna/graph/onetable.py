from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain

from .base import (cgroup_table,
                   get_action_name,
                   make_tracks)
from .onedata import OneDataGraph
from .table import TableGraph, TableGraphRunner, TableGraphWriter
from ..core.table import Table, PositionTable
from ..core.task import dispatch


class OneTableGraph(TableGraph, OneDataGraph, ABC):
    """ Graph of data from one Table. """

    def __init__(self, *,
                 table: Table | PositionTable,
                 **kwargs):
        super().__init__(**kwargs)
        self.table = table

    @property
    def top(self):
        return self.table.top

    @property
    def sample(self):
        return self.table.sample

    @property
    def ref(self):
        return self.table.ref

    @property
    def reg(self):
        return self.table.reg

    @property
    def seq(self):
        return self.table.seq

    @cached_property
    def action(self):
        return get_action_name(self.table)

    @cached_property
    def row_tracks(self):
        return make_tracks(self.table, self.k, self.clust)


class OneTableWriter(TableGraphWriter, ABC):

    def __init__(self, table: Table):
        super().__init__(table)

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return self.tables[0]

    @abstractmethod
    def get_graph(self, *args, **kwargs) -> OneTableGraph:
        """ Return a graph instance. """

    def iter_graphs(self,
                    rels: tuple[str, ...],
                    cgroup: str,
                    **kwargs):
        for cparams in cgroup_table(self.table, cgroup):
            for rels_group in rels:
                yield self.get_graph(rels_group, **kwargs | cparams)


class OneTableRunner(TableGraphRunner, ABC):

    @classmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            max_procs: int,
            **kwargs):
        # Generate a table writer for each table.
        writer_type = cls.get_writer_type()
        writers = [writer_type(table_file)
                   for table_file in cls.list_table_files(input_path)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    pass_n_procs=False,
                                    kwargs=kwargs)))

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
