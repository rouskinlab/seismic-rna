from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from pathlib import Path

from .base import (GraphBase,
                   GraphRunner,
                   GraphWriter,
                   get_source_name,
                   make_source_sample,
                   make_subject)
from ..core.parallel import dispatch
from ..table.base import Table, PosTable
from ..table.load import find_table_files, load


class OneTableGraph(GraphBase, ABC):
    """ Graph of one Table. """

    def __init__(self, *,
                 table: Table | PosTable,
                 order: int | None,
                 clust: int | None,
                 **kwargs):
        super().__init__(**kwargs)
        self.table = table
        self.order = order
        self.clust = clust

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

    @cached_property
    def source(self):
        return get_source_name(self.table)

    @cached_property
    def source_sample(self):
        return make_source_sample(self.source, self.sample)

    @cached_property
    def row_index(self):
        return self.table.header.select(order=self.order, clust=self.clust)

    @property
    def col_index(self):
        return None

    @cached_property
    def subject(self):
        return make_subject(self.source_sample, self.order, self.clust)

    @cached_property
    def data(self):
        return self._fetch_data(self.table,
                                order=self.order,
                                clust=self.clust)


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


class OneTableRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def writer_type(cls) -> type[OneTableWriter]:
        """ Type of GraphWriter. """

    @classmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            max_procs: int,
            parallel: bool,
            **kwargs):
        # List all table files.
        table_files = list(find_table_files(input_path))
        # Generate a table writer for each table.
        writers = [cls.writer_type()(table_file)
                   for table_file in table_files]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    parallel,
                                    pass_n_procs=False,
                                    kwargs=kwargs)))
