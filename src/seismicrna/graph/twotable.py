from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain, combinations, product
from logging import getLogger
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from .base import (LINKER,
                   GraphBase,
                   GraphRunner,
                   GraphWriter,
                   get_source_name,
                   make_subject)
from .seq import get_table_params
from ..core.arg import opt_comppair, opt_compself
from ..core.header import make_header
from ..core.parallel import dispatch
from ..core.seq import POS_NAME
from ..table.base import ClustTable, PosTable, Table
from ..table.load import find_table_files, load

logger = getLogger(__name__)

# Sample index level name
SAMPLE_NAME = "Sample"
ROW_NAME = "Row"
COL_NAME = "Column"


def _get_clusts(max_order: int,
                min_order: int,
                order: int | None,
                clust: int | None):
    return (make_header(max_order=max_order,
                        min_order=min_order).select(order=order,
                                                    clust=clust)
            if max_order > 0
            else None)


# Base Sequence Pair Graph #############################################

class TwoTableGraph(GraphBase, ABC):
    """ Graph of two Tables. """

    def __init__(self, *,
                 table1: Table | PosTable,
                 table2: Table | PosTable,
                 order1: int | None,
                 clust1: int | None,
                 order2: int | None,
                 clust2: int | None,
                 **kwargs):
        super().__init__(**kwargs)
        self.table1 = table1
        self.table2 = table2
        self.order1 = order1
        self.clust1 = clust1
        self.order2 = order2
        self.clust2 = clust2

    @property
    def rel_code(self):
        """ Code of the relationship. """
        if len(self.rel_codes) != 1:
            raise ValueError("Expected exactly one relationship, "
                             f"but got {list(self.rel_codes)}")
        return self.rel_codes

    def _get_common_attribute(self, name: str):
        """ Get the common attribute for tables 1 and 2. """
        attr1 = getattr(self.table1, name)
        attr2 = getattr(self.table2, name)
        if attr1 != attr2:
            raise ValueError(f"Attribute {repr(name)} differs between "
                             f"tables 1 ({repr(attr1)}) and 2 ({repr(attr2)})")
        return attr1

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
    def source1(self):
        """ Source of dataset 1. """
        return get_source_name(self.table1)

    @cached_property
    def source2(self):
        """ Source of dataset 2. """
        return get_source_name(self.table2)

    @cached_property
    def sample_source1(self):
        """ Sample and source of dataset 1. """
        return f"{self.source1} reads from sample {repr(self.sample1)}"

    @cached_property
    def sample_source2(self):
        """ Sample and source of dataset 2. """
        return f"{self.source2} reads from sample {repr(self.sample2)}"

    @cached_property
    def sample_source(self):
        return (self.sample_source1
                if self.sample_source1 == self.sample_source2
                else " vs. ".join([self.sample_source1, self.sample_source2]))

    @cached_property
    def subject1(self):
        """ Name of subject 1. """
        return (make_subject(self.source1, self.order1, self.clust1)
                if isinstance(self.table1, ClustTable)
                else self.source1)

    @cached_property
    def subject2(self):
        """ Name of subject 2. """
        return (make_subject(self.source2, self.order2, self.clust2)
                if isinstance(self.table2, ClustTable)
                else self.source2)

    @cached_property
    def subject(self):
        return (self.subject1 if self.subject1 == self.subject2
                else LINKER.join([self.subject1, self.subject2]))

    @cached_property
    def data1(self):
        """ Data from table 1. """
        return self._fetch_data(self.table1,
                                order=self.order1,
                                clust=self.clust1)

    @cached_property
    def data2(self):
        """ Data from table 2. """
        return self._fetch_data(self.table2,
                                order=self.order2,
                                clust=self.clust2)

    @cached_property
    def row_index(self):
        return _get_clusts(self.table2.header.max_order,
                           self.table2.header.min_order,
                           self.order2,
                           self.clust2)

    @cached_property
    def col_index(self):
        return _get_clusts(self.table1.header.max_order,
                           self.table1.header.min_order,
                           self.order1,
                           self.clust1)


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
    def x_title(self) -> str:
        return POS_NAME

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
        for (row, col), series in self.data.items():
            for trace in self._trace_function()(series, **self._trace_kwargs):
                yield (row, col), trace


class TwoTableWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for two given tables. """

    @classmethod
    @abstractmethod
    def graph_type(cls) -> type[TwoTableGraph]:
        """ Type of the graph to write. """

    def __init__(self, table1_file: Path, table2_file: Path):
        super().__init__(table1_file, table2_file)

    @cached_property
    def table1(self):
        """ The first table providing the data for the graph(s). """
        return load(self.table_files[0])

    @cached_property
    def table2(self):
        """ The second table providing the data for the graph(s). """
        return load(self.table_files[1])

    def iter(self,
             rels: tuple[str, ...],
             arrange: str,
             **kwargs):
        _, _, csparams1 = get_table_params(self.table1, arrange)
        _, _, csparams2 = get_table_params(self.table2, arrange)
        for cparams1, cparams2 in product(csparams1, csparams2):
            for rel in rels:
                yield self.graph_type()(rels=rel,
                                        table1=self.table1,
                                        table2=self.table2,
                                        order1=cparams1.get("order"),
                                        clust1=cparams1.get("clust"),
                                        order2=cparams2.get("order"),
                                        clust2=cparams2.get("clust"),
                                        **kwargs)


class TwoTableRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def writer_type(cls) -> type[TwoTableWriter]:
        """ Type of GraphWriter. """

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_comppair, opt_compself]

    @classmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            compself: bool,
            comppair: bool,
            max_procs: int,
            parallel: bool,
            **kwargs):
        # List all table files.
        table_files = list(find_table_files(input_path))
        # Determine all pairs of tables to compare.
        table_pairs = list()
        if compself:
            # Compare every table with itself.
            table_pairs.extend((file, file) for file in table_files)
        if comppair:
            # Compare every pair of two different tables.
            table_pairs.extend(combinations(table_files, 2))
        # Generate a table writer for each pair of tables.
        writers = [cls.writer_type()(table1_file=t1, table2_file=t2)
                   for t1, t2 in table_pairs]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    parallel,
                                    pass_n_procs=False,
                                    kwargs=kwargs)))

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
