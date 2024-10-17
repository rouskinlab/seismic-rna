from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
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
                   UNAMB_REL,
                   TABLE_RELS,
                   Table,
                   PosTable,
                   ReadTable,
                   AbundanceTable)
from ..batch import END5_COORD, END3_COORD, accumulate_batches
from ..data import MutsDataset
from ..header import Header, make_header
from ..logs import logger
from ..rel import RelPattern, HalfRelPattern
from ..unbias import calc_p_ends_observed
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
    """ Base class for tabulating data for multiple tables from a report
    loader. """

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
        table.loc[:, UNAMB_REL] = (table.loc[:, MATCH_REL].values
                                   +
                                   table.loc[:, MUTAT_REL].values)
        return table

    def __init__(self, dataset: MutsDataset):
        self.dataset = dataset

    @property
    def top(self):
        """ Top-level directory. """
        return self.dataset.top

    @property
    def sample(self):
        """ Name of the sample. """
        return self.dataset.sample

    @property
    def ref(self):
        """ Name of the reference. """
        return self.dataset.ref

    @property
    def section(self):
        """ Section of the dataset. """
        return self.dataset.section

    @property
    def refseq(self):
        """ Reference sequence. """
        return self.dataset.refseq

    @property
    @abstractmethod
    def ks(self) -> list[int] | None:
        """ Numbers of clusters. """

    @cached_property
    def pos_header(self):
        """ Header of the per-position data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def read_header(self):
        """ Header of the per-read data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def _counts(self):
        return accumulate_batches(self.dataset.iter_batches(),
                                  self.dataset.refseq,
                                  self.dataset.section.unmasked_int,
                                  all_patterns(self.dataset.pattern),
                                  ks=self.ks,
                                  validate=False)

    @property
    def _num_reads(self):
        """ Raw number of reads. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return num_reads

    @property
    def _counts_per_pos(self):
        """ Raw counts per position. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return per_pos

    @property
    def _counts_per_read(self):
        """ Raw counts per read. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return per_read

    @property
    def _end_counts(self):
        """ Raw counts for each pair of end coordinates. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return end_counts

    @cached_property
    def p_ends_given_clust_noclose(self):
        """ Probability of each end coordinate. """
        # Ensure end_counts has 2 dimensions.
        if self._end_counts.ndim == 1:
            end_counts = self._end_counts.values[:, np.newaxis]
        else:
            end_counts = self._end_counts.values
        end5s = (self._end_counts.index.get_level_values(END5_COORD).values
                 - self.section.end5)
        end3s = (self._end_counts.index.get_level_values(END3_COORD).values
                 - self.section.end5)
        return calc_p_ends_observed(self.section.length,
                                    end5s,
                                    end3s,
                                    end_counts)

    @cached_property
    def data_per_pos(self):
        """ DataFrame of per-position data. """
        return self._format_table(self._counts_per_pos,
                                  self.pos_header).reindex(self.section.range)

    @cached_property
    def data_per_read(self):
        """ DataFrame of per-read data. """
        return self._format_table(self._counts_per_read,
                                  self.read_header)

    def generate_tables(self, *,
                        table_pos: bool = True,
                        table_read: bool = True,
                        table_clust: bool = True):
        """ Generate the tables from this data. """
        for table_type in self.table_types():
            if issubclass(table_type, PosTableWriter):
                if table_pos:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, ReadTableWriter):
                if table_read:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, AbundanceTableWriter):
                if table_clust:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            else:
                # This should never happen; checking just in case.
                raise TypeError(table_type)

    def write_tables(self, *, force: bool = False, **kwargs):
        for table in self.generate_tables(**kwargs):
            table.write(force)


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
    def sect(self):
        return self.tabulator.section.name

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


class PosTableWriter(TableWriter, PosTable, ABC):

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
