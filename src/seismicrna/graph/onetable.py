from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from pathlib import Path
from typing import Iterable

from .base import get_action_name
from .cgroup import (ClusterGroupRunner,
                     cgroup_table,
                     make_tracks)
from .onesource import OneSourceGraph, OneSourceClusterGroupGraph
from .table import (TableGraph,
                    TableRunner,
                    TableWriter,
                    RelTableGraph,
                    RelTableRunner)
from ..core.table import Table, PositionTable, AbundanceTable
from ..core.task import dispatch


class OneTableGraph(TableGraph, OneSourceGraph, ABC):
    """ Graph of data from one Table. """

    def __init__(self, *,
                 table: Table | PositionTable | AbundanceTable,
                 **kwargs):
        super().__init__(**kwargs)
        self.table = table

    @property
    def top(self):
        return self.table.top

    @property
    def branches(self):
        return self.table.branches

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


class OneTableRelClusterGroupGraph(OneTableGraph,
                                   RelTableGraph,
                                   OneSourceClusterGroupGraph,
                                   ABC):

    @cached_property
    def row_tracks(self):
        return make_tracks(self.table,
                           self.k,
                           self.clust,
                           k_clust_list=self.k_clust_list)


class OneTableWriter(TableWriter, ABC):

    def __init__(self, table: Table, **kwargs):
        super().__init__(table, **kwargs)

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        assert len(self.tables) == 1
        return self.tables[0]

    @abstractmethod
    def get_graph(self, *args, **kwargs) -> OneTableGraph:
        """ Return a graph instance. """


class OneTableRelClusterGroupWriter(OneTableWriter, ABC):

    def iter_graphs(self, *, rels: list[str], cgroup: str, **kwargs):
        for cparams in cgroup_table(self.table, cgroup):
            for rels_group in rels:
                yield self.get_graph(rels_group, **kwargs | cparams)


class OneTableRunner(TableRunner, ABC):

    @classmethod
    def run(cls,
            input_path: Iterable[str | Path], *,
            verify_times: bool,
            num_cpus: int,
            **kwargs):
        # Generate a table writer for each table.
        writer_type = cls.get_writer_type()
        writers = [writer_type(table_file)
                   for table_file
                   in cls.load_input_files(input_path,
                                           verify_times=verify_times)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    num_cpus=num_cpus,
                                    pass_num_cpus=False,
                                    as_list=False,
                                    ordered=False,
                                    raise_on_error=False,
                                    kwargs=kwargs)))


class OneTableRelClusterGroupRunner(OneTableRunner,
                                    RelTableRunner,
                                    ClusterGroupRunner,
                                    ABC):
    pass
