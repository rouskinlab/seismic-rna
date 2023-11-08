from abc import ABC
from logging import getLogger
from pathlib import Path

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
                   ClustReadTable,
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
            raise ValueError(f"{type(self).__name__} got path {table_file},"
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

    @property
    def _data(self):
        data = (pd.read_csv(self.path,
                            index_col=self.header_rows(),
                            header=self.index_cols()).T
                if self.transposed() else
                pd.read_csv(self.path,
                            index_col=self.index_cols(),
                            header=self.header_rows()))
        # Any numeric data in the header will be read as strings and
        # must be cast to integers, which can be done with parse_header.
        header = parse_header(data.columns)
        # The columns must be replaced with the header index explicitly
        # in order for the type casting to take effect.
        data.columns = header.index
        return data


class RelTypeTableLoader(TableLoader, RelTypeTable, ABC):
    """ Load a table of relationship types. """


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


class ClustReadTableLoader(ReadTableLoader, ClustReadTable):
    """ Load cluster data indexed by read. """


class ClustFreqTableLoader(TableLoader, ClustFreqTable):
    """ Load cluster data indexed by cluster. """


# Helper Functions #####################################################

def load(table_file: Path):
    """ Helper function to load a TableLoader from a table file. """
    for loader_type in (RelPosTableLoader,
                        RelReadTableLoader,
                        MaskPosTableLoader,
                        MaskReadTableLoader,
                        ClustPosTableLoader,
                        ClustReadTableLoader,
                        ClustFreqTableLoader):
        try:
            # Try to load the table file using the type of TableLoader.
            return loader_type(table_file)
        except (FileNotFoundError, ValueError):
            # The TableLoader was not the correct type for the file.
            pass
    # None of the TableLoader types were able to open the file.
    raise ValueError(f"Failed to open table: {table_file}")


def find_tables(tables: tuple[str, ...]):
    """ Yield a file for each given file/directory of a table. """
    yield from path.find_files_chain(map(Path, tables), [path.TableSeg])


def load_tables(tables: tuple[str, ...]):
    """ Yield a table for each given file/directory of a table. """
    for file in find_tables(tables):
        try:
            yield load(file)
        except Exception as error:
            logger.error(f"Failed to load table from {file}: {error}")

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
