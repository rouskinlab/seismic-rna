from abc import ABC
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .base import (Table,
                   RelTypeTable,
                   PosTable,
                   ReadTable,
                   RelPosTable,
                   RelReadTable,
                   MaskPosTable,
                   MaskReadTable,
                   ClustPosTable,
                   ClustFreqTable)
from ..core import path
from ..core.header import parse_header

logger = getLogger(__name__)


# Table Loader Base Classes ############################################

class TableLoader(Table, ABC):
    """ Load a table from a file. """

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        self._sect = fields[path.SECT]
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
    def sect(self) -> str:
        return self._sect


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


# Load by Index (position/read/frequency) ##############################

class PosTableLoader(RelTypeTableLoader, PosTable, ABC):
    """ Load data indexed by position. """


class ReadTableLoader(RelTypeTableLoader, ReadTable, ABC):
    """ Load data indexed by read. """


# Instantiable Table Loaders ###########################################

class RelPosTableLoader(PosTableLoader, RelPosTable):
    """ Load relation data indexed by position. """


class RelReadTableLoader(ReadTableLoader, RelReadTable):
    """ Load relation data indexed by read. """


class MaskPosTableLoader(PosTableLoader, MaskPosTable):
    """ Load masked bit vector data indexed by position. """


class MaskReadTableLoader(ReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClustPosTableLoader(PosTableLoader, ClustPosTable):
    """ Load cluster data indexed by position. """


class ClustFreqTableLoader(TableLoader, ClustFreqTable):
    """ Load cluster data indexed by cluster. """

    @cached_property
    def data(self) -> pd.Series:
        data = pd.read_csv(self.path,
                           index_col=self.index_cols()).squeeze(axis=1)
        if not isinstance(data, pd.Series):
            raise ValueError(f"{self} must have one column, but got\n{data}")
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.index)
        # The index must be replaced with the header index for the
        # type casting to take effect.
        data.index = header.index
        return data


# Helper Functions #####################################################


POS_LOADS = RelPosTableLoader, MaskPosTableLoader, ClustPosTableLoader
READ_LOADS = RelReadTableLoader, MaskReadTableLoader
FREQ_LOADS = ClustFreqTableLoader,


def load_table(types: Iterable[type[PosTableLoader
                                    | ReadTableLoader
                                    | ClustFreqTableLoader]],
               table_file: Path):
    """ Load a Table of one of several types from a file. """
    for table_type in types:
        try:
            # Try to load the table file using the type of TableLoader.
            return table_type(table_file)
        except (FileNotFoundError, ValueError):
            # The TableLoader was not the correct type for the file.
            pass
    raise TypeError(f"Failed to open {table_file} as one of {types}")


def load_any_table(table_file: Path):
    return load_table([*POS_LOADS, *READ_LOADS, *FREQ_LOADS], table_file)


def load_pos_table(table_file: Path) -> PosTableLoader:
    return load_table(POS_LOADS, table_file)


def load_read_table(table_file: Path) -> ReadTableLoader:
    return load_table(READ_LOADS, table_file)


def find_tables(segments: Iterable[path.Segment], files: Iterable[str | Path]):
    """ Yield every table file with the given type of segment from among
    the given paths. """
    for segment in segments:
        yield from path.find_files_chain(files, [segment])


def find_all_tables(files: Iterable[str | Path]):
    yield from find_tables([path.PosTableSeg,
                            path.ReadTableSeg,
                            path.FreqTableSeg],
                           files)


def find_pos_tables(files: Iterable[str | Path]):
    yield from find_tables([path.PosTableSeg], files)


def find_read_tables(files: Iterable[str | Path]):
    yield from find_tables([path.ReadTableSeg], files)


def load_tables(finder: Callable[[Iterable[str | Path]], Iterable[Path]],
                loader: Callable[[Path], TableLoader],
                files: Iterable[str | Path]):
    """ Yield every table with the given type of segment from among the
    given paths. """
    for file in finder(files):
        try:
            yield loader(file)
        except TypeError as error:
            logger.error(error)


def load_all_tables(files: Iterable[str | Path]):
    """ Yield every table among the given paths. """
    yield from load_tables(find_all_tables, load_any_table, files)


def load_pos_tables(files: Iterable[str | Path]):
    """ Yield every positional table among the given paths. """
    yield from load_tables(find_pos_tables, load_pos_table, files)


def load_read_tables(files: Iterable[str | Path]):
    """ Yield every per-read table among the given paths. """
    yield from load_tables(find_read_tables, load_read_table, files)

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
