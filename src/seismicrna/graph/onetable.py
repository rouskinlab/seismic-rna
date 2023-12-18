from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from pathlib import Path

from .base import (GraphBase,
                   GraphRunner,
                   GraphWriter,
                   arrange_table,
                   get_source_name,
                   make_index,
                   make_source_sample,
                   make_subject)
from ..core.header import parse_header
from ..core.parallel import dispatch
from ..core.rna import RNAProfile, RNAState, from_ct
from ..table.base import Table, PosTable, get_rel_name
from ..table.load import find_table_files, load


class OneTableGraph(GraphBase, ABC):
    """ Graph of data from one Table. """

    def __init__(self, *,
                 table: Table | PosTable,
                 order: int | None,
                 clust: int | None,
                 **kwargs):
        super().__init__(**kwargs)
        self.table = table
        self.order = order
        self.clust = clust

    @cached_property
    def rel_names(self):
        return list(map(get_rel_name, self.rel_codes))

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
        return make_index(self.table.header, self.order, self.clust)

    @property
    def col_index(self):
        return None

    @cached_property
    def subject(self):
        return make_subject(self.source, self.order, self.clust)

    @cached_property
    def data(self):
        return self._fetch_data(self.table,
                                order=self.order,
                                clust=self.clust)

    @cached_property
    def data_header(self):
        """ Header of the selected data (not of the entire table). """
        return parse_header(self.data.columns)


class OneTableStructureGraph(OneTableGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *, ct_file: Path | None = None, **kwargs):
        super().__init__(**kwargs)
        self.ct_file = ct_file

    @cached_property
    def rel_name(self):
        """ Name of the relationship to graph. """
        return get_rel_name(self.rel_codes)

    def iter_profiles(self):
        """ Yield each RNAProfile from the table. """
        yield from self.table.iter_profiles(quantile=self.quantile,
                                            rel=self.rel_name,
                                            order=self.order,
                                            clust=self.clust)

    def iter_states(self):
        """ Yield each RNAState. """
        for profile in self.iter_profiles():
            ct_file = (self.ct_file
                       if self.ct_file is not None
                       else profile.get_ct_file(self.top))
            for struct in from_ct(ct_file):
                yield RNAState.from_struct_profile(struct, profile)


class OneTableWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for a given table. """

    @classmethod
    @abstractmethod
    def get_graph_type(cls, *args, **kwargs) -> type[OneTableGraph]:
        """ Type of the graph to write. """

    def __init__(self, table_file: Path):
        super().__init__(table_file)

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return load(self.table_files[0])

    def iter(self,
             rels: tuple[str, ...],
             arrange: str,
             **kwargs):
        for cparams in arrange_table(self.table, arrange):
            for rel in rels:
                graph_type = self.get_graph_type(rel)
                yield graph_type(rels=rel,
                                 table=self.table,
                                 order=cparams.get("order"),
                                 clust=cparams.get("clust"),
                                 **kwargs)


class OneTableRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def writer_type(cls, *args, **kwargs) -> type[OneTableWriter]:
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
