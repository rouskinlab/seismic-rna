from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from inspect import Parameter, signature
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from .base import (COVER_REL,
                   DELET_REL,
                   INSRT_REL,
                   MATCH_REL,
                   MUTAT_REL,
                   SUBST_REL,
                   SUB_A_REL,
                   SUB_C_REL,
                   SUB_G_REL,
                   SUB_T_REL,
                   INFOR_REL,
                   TABLE_RELS,
                   Table,
                   PositionTable,
                   ReadTable,
                   AbundanceTable)
from ..batch import accumulate_batches, accumulate_counts
from ..data import MutsDataset
from ..header import Header, make_header
from ..logs import logger
from ..rel import RelPattern, HalfRelPattern
from ..seq import Region
from ..write import need_write

# These relationships are of all subtypes of mutations.
SUBMUTS = [SUBST_REL,
           SUB_A_REL,
           SUB_C_REL,
           SUB_G_REL,
           SUB_T_REL,
           DELET_REL,
           INSRT_REL]

PRECISION = 1


class Tabulator(ABC):
    """ Base class for tabulating data for one or more tables. """

    @classmethod
    @abstractmethod
    def table_types(cls) -> list[type[TableWriter]]:
        """ Types of tables that this tabulator can write. """

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    @classmethod
    def _format_table(cls, counts: pd.DataFrame, header: Header):
        """ Format the count data with the proper header. """
        # Initialize an empty table.
        table = pd.DataFrame(cls.get_null_value(), counts.index, header.index)
        # Count reads with each relationship at each position.
        for rel, rel_counts in counts.items():
            table.loc[:, rel] = rel_counts
        # Add a column for informative relationships.
        table.loc[:, INFOR_REL] = (table.loc[:, MATCH_REL].values
                                   +
                                   table.loc[:, MUTAT_REL].values)
        return table

    def __init__(self, *,
                 top: Path,
                 sample: str,
                 region: Region,
                 count_pos: bool,
                 count_read: bool,
                 validate: bool = True):
        self.top = top
        self.sample = sample
        self.refseq = region.seq
        self.region = region
        self.pattern = None
        self.ks = None
        self.count_ends = False
        self.count_pos = count_pos
        self.count_read = count_read
        self.validate = validate

    @property
    def ref(self):
        """ Name of the reference. """
        return self.region.ref

    @cached_property
    def pos_header(self):
        """ Header of the per-position data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def read_header(self):
        """ Header of the per-read data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def _accum_kwargs(self):
        """ Keyword arguments for the accumulate function. """
        return dict(refseq=self.refseq,
                    pos_nums=self.region.unmasked_int,
                    patterns=all_patterns(self.pattern),
                    ks=self.ks,
                    count_ends=self.count_ends,
                    count_pos=self.count_pos,
                    count_read=self.count_read,
                    validate=self.validate)

    @cached_property
    @abstractmethod
    def _counts(self):
        """ All counts for all table(s). """

    @property
    def num_reads(self):
        """ Raw number of reads. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return num_reads

    @property
    def counts_per_pos(self):
        """ Raw counts per position. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return per_pos

    @property
    def counts_per_read(self):
        """ Raw counts per read. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return per_read

    @property
    def end_counts(self):
        """ Raw counts for each pair of end coordinates. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return end_counts

    @cached_property
    def data_per_pos(self):
        """ DataFrame of per-position data. """
        return self._format_table(self.counts_per_pos,
                                  self.pos_header).reindex(self.region.range)

    @cached_property
    def data_per_read(self):
        """ DataFrame of per-read data. """
        return self._format_table(self.counts_per_read,
                                  self.read_header)

    @cached_property
    @abstractmethod
    def data_per_clust(self) -> pd.Series | None:
        """ Series of per-cluster data (or None if no clusters). """

    def generate_tables(self, *,
                        pos: bool = True,
                        read: bool = True,
                        clust: bool = True):
        """ Generate the tables from this data. """
        for table_type in self.table_types():
            if issubclass(table_type, PositionTableWriter):
                if pos:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, ReadTableWriter):
                if read:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, AbundanceTableWriter):
                if clust:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            else:
                # This should never happen; checking just in case.
                raise TypeError(table_type)

    def write_tables(self, *, force: bool = False, **kwargs):
        files = list()
        for table in self.generate_tables(**kwargs):
            files.append(table.write(force))
        return files


class CountTabulator(Tabulator, ABC):
    """ Tabulator that accepts pre-counted data from batches. """

    def __init__(self, *,
                 batch_counts: Iterable[tuple[Any, Any, Any, Any]],
                 **kwargs):
        super().__init__(**kwargs)
        self._batch_counts = iter(batch_counts)

    @cached_property
    def _counts(self):
        return accumulate_counts(self._batch_counts, **self._accum_kwargs)


class BatchTabulator(Tabulator, ABC):
    """ Tabulator that accepts batches as the input data. """

    def __init__(self, *,
                 get_batch_count_all: Callable,
                 num_batches: int,
                 max_procs: int = 1,
                 **kwargs):
        super().__init__(**kwargs)
        self._get_batch_count_all = get_batch_count_all
        self.num_batches = num_batches
        self.max_procs = max_procs

    @cached_property
    def _counts(self):
        return accumulate_batches(self._get_batch_count_all,
                                  self.num_batches,
                                  max_procs=self.max_procs,
                                  **self._accum_kwargs)


class DatasetTabulator(BatchTabulator, ABC):
    """ Tabulator made from one dataset. """

    @classmethod
    def _list_args(cls, func: Callable):
        """ List the positional arguments of a function. """
        return [name for name, param in signature(func).parameters.items()
                if param.kind is Parameter.KEYWORD_ONLY]

    @classmethod
    @cache
    def _init_data(cls):
        """ Attributes of the dataset to use as keyword arguments in
        super().__init__(). """
        return ["top", "sample", "get_batch_count_all", "num_batches"]

    def __init__(self, *,
                 dataset: MutsDataset,
                 validate: bool = False,
                 **kwargs):
        # Since the batches come from a Dataset, they do not need to be
        # validated, so make validate False by default.
        super().__init__(**{attr: getattr(dataset, attr)
                            for attr in self._init_data()},
                         validate=validate,
                         **kwargs)


class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: Tabulator):
        self.tabulator = tabulator

    @property
    def top(self):
        return self.tabulator.top

    @property
    def sample(self):
        return self.tabulator.sample

    @property
    def ref(self):
        return self.tabulator.ref

    @property
    def reg(self):
        return self.tabulator.region.name

    @property
    def refseq(self):
        return self.tabulator.refseq

    @property
    def columns(self):
        return self.header.index

    def write(self, force: bool):
        """ Write the table's rounded data to the table's CSV file. """
        if need_write(self.path, force):
            self.data.round(decimals=PRECISION).to_csv(self.path)
            logger.routine(f"Wrote {self} to {self.path}")
        return self.path


class PositionTableWriter(TableWriter, PositionTable, ABC):

    @cached_property
    def data(self):
        return self.tabulator.data_per_pos


class ReadTableWriter(TableWriter, ReadTable, ABC):

    @cached_property
    def data(self):
        return self.tabulator.data_per_read


class AbundanceTableWriter(TableWriter, AbundanceTable, ABC):

    @cached_property
    def data(self):
        return self.tabulator.data_per_clust


def _iter_mut_patterns():
    """ Yield a HalfRelPattern for each type of mutation. """
    yield SUBST_REL, HalfRelPattern.from_counts(count_sub=True)
    yield SUB_A_REL, HalfRelPattern("ca", "ga", "ta")
    yield SUB_C_REL, HalfRelPattern("ac", "gc", "tc")
    yield SUB_G_REL, HalfRelPattern("ag", "cg", "tg")
    yield SUB_T_REL, HalfRelPattern("at", "ct", "gt")
    yield DELET_REL, HalfRelPattern.from_counts(count_del=True)
    yield INSRT_REL, HalfRelPattern.from_counts(count_ins=True)


def _iter_patterns(mask: RelPattern | None = None):
    """ Yield a RelPattern for every type of relationship. """
    # Count everything except for no coverage.
    yield COVER_REL, RelPattern.allc()
    # Count matches to the reference sequence.
    yield MATCH_REL, RelPattern.muts().intersect(mask, invert=True)
    # Count all types of mutations, relative to reference matches.
    yield MUTAT_REL, RelPattern.muts().intersect(mask)
    # Count each type of mutation, relative to reference matches.
    for mut, pattern in _iter_mut_patterns():
        yield mut, RelPattern(pattern, HalfRelPattern.refs()).intersect(mask)


def all_patterns(mask: RelPattern | None = None):
    """ Every RelPattern, keyed by its name. """
    return dict(_iter_patterns(mask))

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
