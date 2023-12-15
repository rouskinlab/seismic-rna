from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path

from .base import GraphBase, GraphWriter, get_source_name
from ..table.base import Table, PosTable
from ..table.load import load


class OneTableGraph(GraphBase, ABC):
    """ Graph of one Table. """

    def __init__(self, *, table: Table | PosTable, **kwargs):
        super().__init__(**kwargs)
        self.table = table

    @property
    def sample_source(self):
        return get_source_name(self.table)

    @property
    def sample(self):
        return self.table.sample

    @property
    def ref(self):
        return self.table.ref

    @property
    def sect(self):
        return self.table.sect

    @property
    def seq(self):
        return self.table.seq

    @property
    def col_index(self):
        return None


class OneTableWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for a given table. """

    @classmethod
    @abstractmethod
    def graph_type(cls) -> type[OneTableGraph]:
        """ Type of the graph to write. """

    def __init__(self, table_file: Path):
        super().__init__(table_file)

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return load(self.table_files[0])
