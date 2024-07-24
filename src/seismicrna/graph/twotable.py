from abc import ABC, abstractmethod
from collections import defaultdict
from functools import cached_property
from itertools import chain, combinations, product
from logging import getLogger
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from .base import (LINKER,
                   GraphRunner,
                   GraphWriter,
                   cgroup_table,
                   get_action_name,
                   make_tracks,
                   make_title_action_sample,
                   make_path_subject)
from .rel import OneRelGraph
from ..core import path
from ..core.arg import opt_comppair, opt_compself, opt_out_dir
from ..core.task import dispatch
from ..table.base import ClustTable, PosTable, Table
from ..table.load import find_pos_tables, load_pos_table

logger = getLogger(__name__)

# Index level names.
SAMPLE_NAME = "Sample"
ROW_NAME = "Row"
COL_NAME = "Column"


class TwoTableGraph(OneRelGraph, ABC):
    """ Graph of two Tables. """

    def __init__(self, *,
                 out_dir: str | Path,
                 table1: Table | PosTable,
                 k1: int | None,
                 clust1: int | None,
                 table2: Table | PosTable,
                 k2: int | None,
                 clust2: int | None,
                 **kwargs):
        super().__init__(**kwargs)
        self._top = Path(out_dir)
        self.table1 = table1
        self.k1 = k1
        self.clust1 = clust1
        self.table2 = table2
        self.k2 = k2
        self.clust2 = clust2

    def _get_common_attribute(self, name: str):
        """ Get the common attribute for tables 1 and 2. """
        attr1 = getattr(self.table1, name)
        attr2 = getattr(self.table2, name)
        if attr1 != attr2:
            raise ValueError(f"Attribute {repr(name)} differs between "
                             f"tables 1 ({repr(attr1)}) and 2 ({repr(attr2)})")
        return attr1

    @property
    def top(self):
        return self._top

    @property
    def sample1(self):
        """ Name of sample 1. """
        return self.table1.sample

    @property
    def sample2(self):
        """ Name of sample 2. """
        return self.table2.sample

    @cached_property
    def sample(self):
        return (self.sample1 if self.sample1 == self.sample2
                else LINKER.join([self.sample1, self.sample2]))

    @cached_property
    def ref(self):
        return self._get_common_attribute("ref")

    @cached_property
    def sect(self):
        return self._get_common_attribute("sect")

    @cached_property
    def seq(self):
        return self._get_common_attribute("seq")

    @cached_property
    def action1(self):
        """ Action that generated dataset 1. """
        return get_action_name(self.table1)

    @cached_property
    def action2(self):
        """ Action that generated dataset 2. """
        return get_action_name(self.table2)

    @cached_property
    def action_sample1(self):
        """ Action and sample of dataset 1. """
        return make_title_action_sample(self.action1, self.sample1)

    @cached_property
    def action_sample2(self):
        """ Action and sample of dataset 2. """
        return make_title_action_sample(self.action2, self.sample2)

    @cached_property
    def title_action_sample(self):
        return (self.action_sample1
                if self.action_sample1 == self.action_sample2
                else " vs. ".join([self.action_sample1, self.action_sample2]))

    @cached_property
    def path_subject1(self):
        """ Name of subject 1. """
        return (make_path_subject(self.action1, self.k1, self.clust1)
                if isinstance(self.table1, ClustTable)
                else self.action1)

    @cached_property
    def path_subject2(self):
        """ Name of subject 2. """
        return (make_path_subject(self.action2, self.k2, self.clust2)
                if isinstance(self.table2, ClustTable)
                else self.action2)

    @cached_property
    def path_subject(self):
        return (self.path_subject1
                if self.path_subject1 == self.path_subject2
                else LINKER.join([self.path_subject1, self.path_subject2]))

    @cached_property
    def data1(self):
        """ Data from table 1. """
        return self._fetch_data(self.table1,
                                k=self.k1,
                                clust=self.clust1)

    @cached_property
    def data2(self):
        """ Data from table 2. """
        return self._fetch_data(self.table2,
                                k=self.k2,
                                clust=self.clust2)

    @cached_property
    def row_tracks(self):
        return make_tracks(self.table2.header, self.k2, self.clust2)

    @cached_property
    def col_tracks(self):
        return make_tracks(self.table1.header, self.k1, self.clust1)


class TwoTableMergedGraph(TwoTableGraph, ABC):
    """ Graph of a pair of datasets over the same sequence in which the
    data series are merged in some fashion into another series, and the
    original data are not graphed directly. """

    @classmethod
    @abstractmethod
    def _trace_function(cls) -> Callable:
        """ Function to generate the graph's traces. """

    @property
    def _trace_kwargs(self) -> dict[str, Any]:
        """ Keyword arguments for self._trace_function. """
        return dict()

    @property
    @abstractmethod
    def _merge_data(self) -> Callable:
        """ Function to merge the two datasets into one. """

    @cached_property
    def data(self):
        # Merge each pair in the Cartesian product of datasets 1 and 2.
        data = pd.DataFrame.from_dict(
            {(row, col): self._merge_data(vals1, vals2)
             for (col, (key1, vals1)), (row, (key2, vals2))
             in product(enumerate(self.data1.items(), start=1),
                        enumerate(self.data2.items(), start=1))}
        )
        # Indicate that the column index levels now correspond to the
        # rows and columns of the graph.
        data.columns.rename([ROW_NAME, COL_NAME], inplace=True)
        return data

    def get_traces(self):
        trace_func = self._trace_function()
        for (row, col), series in self.data.items():
            for trace in trace_func(series, **self._trace_kwargs):
                yield (row, col), trace


class TwoTableWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for two given tables. """

    @classmethod
    @abstractmethod
    def get_graph_type(cls, *args, **kwargs) -> type[TwoTableGraph]:
        """ Type of graph. """

    def __init__(self, table1_file: Path, table2_file: Path):
        super().__init__(table1_file, table2_file)

    @cached_property
    def table1(self):
        """ The first table providing the data for the graph(s). """
        return load_pos_table(self.table_files[0])

    @cached_property
    def table2(self):
        """ The second table providing the data for the graph(s). """
        return load_pos_table(self.table_files[1])

    def iter_graphs(self,
                    rels: tuple[str, ...],
                    cgroup: str,
                    **kwargs):
        for cparams1, cparams2 in product(cgroup_table(self.table1, cgroup),
                                          cgroup_table(self.table2, cgroup)):
            for rel in rels:
                graph_type = self.get_graph_type()
                yield graph_type(rel=rel,
                                 table1=self.table1,
                                 k1=cparams1["k"],
                                 clust1=cparams1["clust"],
                                 table2=self.table2,
                                 k2=cparams2["k"],
                                 clust2=cparams2["clust"],
                                 **kwargs)


def _iter_table_pairs(table_files: Iterable[Path],
                      table_segs: tuple[path.Segment, ...]):
    """ Yield every pair of files whose reference and section match. """
    logger.debug("Seeking all pairs of table files with identical references "
                 f"and sections matching segments {list(map(str, table_segs))}")
    # Determine the reference and section of each table.
    table_fields = defaultdict(set)
    for file in table_files:
        fields = path.parse(file, *table_segs)
        key = fields[path.REF], fields[path.SECT]
        if file in table_fields[key]:
            logger.warning(f"Duplicate table file: {file}")
        else:
            table_fields[key].add(file)
    # Yield every pair of table files.
    for (ref, sect), tables in table_fields.items():
        n_files = len(tables)
        n_pairs = n_files * (n_files - 1) // 2
        logger.debug(f"Found {n_files} table files ({n_pairs} pairs) with "
                     f"reference {repr(ref)} and section {repr(sect)}")
        yield from combinations(sorted(tables), 2)


def iter_pos_table_pairs(table_files: Iterable[Path]):
    """ Yield every pair of files of positional tables whose reference
    and section match. """
    yield from _iter_table_pairs(table_files, path.POS_TABLE_SEGS)


def iter_read_table_pairs(table_files: Iterable[Path]):
    """ Yield every pair of files of per-read tables whose reference and
    section match. """
    yield from _iter_table_pairs(table_files, path.READ_TABLE_SEGS)


class TwoTableRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[TwoTableWriter]:
        """ Type of GraphWriter. """

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_comppair, opt_compself, opt_out_dir]

    @classmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            compself: bool,
            comppair: bool,
            max_procs: int,
            parallel: bool,
            **kwargs):
        # List all table files.
        table_files = list(find_pos_tables(input_path))
        # Determine all pairs of tables to compare.
        table_pairs = list()
        if compself:
            # Compare every table with itself.
            table_pairs.extend((file, file) for file in table_files)
        if comppair:
            # Compare every pair of two different tables.
            table_pairs.extend(iter_pos_table_pairs(table_files))
        # Generate a table writer for each pair of tables.
        writers = [cls.get_writer_type()(table1_file, table2_file)
                   for table1_file, table2_file in table_pairs]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    parallel,
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
