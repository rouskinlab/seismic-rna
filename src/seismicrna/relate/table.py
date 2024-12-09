from abc import ABC
from functools import cache, cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .data import load_relate_dataset
from .report import RelateReport
from ..core import path
from ..core.header import RelHeader, parse_header
from ..core.logs import logger
from ..core.seq import FULL_NAME, DNA, Region
from ..core.table import (Table,
                          Tabulator,
                          DatasetTabulator,
                          CountTabulator,
                          PositionTable,
                          PositionTableWriter,
                          ReadTable,
                          ReadTableWriter,
                          RelTypeTable)


class AvgTable(RelTypeTable, ABC):
    """ Average over an ensemble of RNA structures. """

    @classmethod
    def header_type(cls):
        return RelHeader


class FullTable(Table, ABC):

    @property
    def path_fields(self):
        return {path.TOP: self.top,
                path.SAMP: self.sample,
                path.REF: self.ref}


class FullPositionTable(FullTable, PositionTable, ABC):

    @classmethod
    def path_segs(cls):
        return path.REF_DIR_SEGS + (path.PositionTableSeg,)


class FullReadTable(FullTable, ReadTable, ABC):

    @classmethod
    def path_segs(cls):
        return path.REF_DIR_SEGS + (path.ReadTableSeg,)


class RelateTable(AvgTable, ABC):

    @classmethod
    def kind(cls):
        return path.CMD_REL_DIR


class RelatePositionTable(RelateTable, FullPositionTable, ABC):

    def _iter_profiles(self, *,
                       regions: Iterable[Region] | None,
                       quantile: float,
                       rel: str,
                       k: int | None,
                       clust: int | None):
        # Relation table loaders have unmasked, unfiltered reads and are
        # thus unsuitable for making RNA profiles. Yield no profiles.
        yield from ()


class RelateReadTable(RelateTable, FullReadTable, ABC):
    pass


class RelatePositionTableWriter(PositionTableWriter, RelatePositionTable):
    pass


class RelateReadTableWriter(ReadTableWriter, RelateReadTable):
    pass


class FullTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return 0

    def __init__(self, *, ref: str, refseq: DNA, **kwargs):
        # For a full tabulator, the full reference sequence must be used
        # as the region.
        super().__init__(region=Region(ref, refseq), **kwargs)


class AverageTabulator(Tabulator, ABC):

    @cached_property
    def data_per_clust(self):
        # An ensemble average tabulator has no per-cluster data.
        return None


class RelateTabulator(FullTabulator, AverageTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [RelatePositionTableWriter, RelateReadTableWriter]


class RelateCountTabulator(CountTabulator, RelateTabulator):
    pass


class RelateDatasetTabulator(DatasetTabulator, RelateTabulator):

    @classmethod
    @cache
    def _init_data(cls):
        return super()._init_data() + cls._list_args(FullTabulator.__init__)


class TableLoader(Table, ABC):
    """ Load a table from a file. """

    @classmethod
    def find_tables(cls, paths: Iterable[str | Path]):
        """ Yield files of the tables within the given paths. """
        for file in path.find_files_chain(paths, cls.path_segs()):
            if file.name.startswith(cls.kind()):
                yield file

    @classmethod
    def load_tables(cls, paths: Iterable[str | Path]):
        """ Yield tables within the given paths. """
        for file in cls.find_tables(paths):
            try:
                yield cls(file)
            except Exception as error:
                logger.error(error)

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        self._reg = fields.get(path.REG, FULL_NAME)
        if not self.path.with_suffix(table_file.suffix).samefile(table_file):
            raise ValueError(f"{type(self).__name__} got path {table_file}, "
                             f"but expected {self.path}")

    @property
    def top(self) -> Path:
        return self._out_dir

    @property
    def sample(self) -> str:
        return self._sample

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def reg(self) -> str:
        return self._reg

    @cached_property
    def refseq(self):
        dataset = load_relate_dataset(RelateReport.build_path(
            top=self.top, sample=self.sample, ref=self.ref)
        )
        return dataset.refseq


class RelTypeTableLoader(TableLoader, RelTypeTable, ABC):
    """ Load a table of relationship types. """

    @cached_property
    def data(self) -> pd.DataFrame:
        data = pd.read_csv(self.path,
                           index_col=self.index_cols(),
                           header=self.header_rows())
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.columns)
        # The columns must be replaced with the header index for the
        # type casting to take effect.
        data.columns = header.index
        return data


class PositionTableLoader(RelTypeTableLoader, PositionTable, ABC):
    """ Load data indexed by position. """


class ReadTableLoader(RelTypeTableLoader, ReadTable, ABC):
    """ Load data indexed by read. """


class RelatePositionTableLoader(PositionTableLoader, RelatePositionTable):
    """ Load relate data indexed by position. """


class RelateReadTableLoader(ReadTableLoader, RelateReadTable):
    """ Load relate data indexed by read. """

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
