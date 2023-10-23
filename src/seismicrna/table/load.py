from abc import ABC, abstractmethod
from functools import cached_property
from itertools import product
from pathlib import Path
from typing import Iterable

import pandas as pd

from .base import (MUTAT_REL,
                   Table,
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
from ..cluster.names import AVERAGE_NAME, fmt_clust_name
from ..core import path
from ..core.batch import OC_INDEX_NAMES, POC_INDEX_NAMES, REL_NAME
from ..core.rna.profile import RnaProfile
from ..core.seq import INDEX_NAMES, Section


# Table Loader Base Classes ############################################

class TableLoader(Table, ABC):
    """ Load a table from a file. """

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        self._sect = fields[path.SECT]
        if not self.path.samefile(table_file):
            raise ValueError(f"Invalid path for '{type(self).__name__}': "
                             f"{table_file} (expected {self.path})")

    @property
    def top(self):
        return self._out_dir

    @property
    def sample(self):
        return self._sample

    @property
    def ref(self):
        return self._ref

    @property
    def sect(self):
        return self._sect

    @classmethod
    @abstractmethod
    def index_col(cls) -> list[int]:
        """ Column(s) of the file to use as the index. """

    @classmethod
    @abstractmethod
    def header_row(cls) -> list[int]:
        """ Row(s) of the file to use as the columns. """

    @cached_property
    @abstractmethod
    def data(self):
        return pd.read_csv(self.path,
                           index_col=self.index_col(),
                           header=self.header_row())


class RelTypeTableLoader(TableLoader, RelTypeTable, ABC):
    """ Load a table of relationship types. """


# Load by Index (position/read/frequency) ##############################

class PosTableLoader(RelTypeTableLoader, PosTable, ABC):
    """ Load data indexed by position. """

    @classmethod
    def index_col(cls):
        return list(range(len(INDEX_NAMES)))

    @abstractmethod
    def iter_profiles(self, sections: Iterable[Section], quantile: float):
        """ Yield RNA mutational profiles from the table. """
        for section in sections:
            yield RnaProfile(section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=pd.Series())


class ReadTableLoader(RelTypeTableLoader, ReadTable, ABC):
    """ Load data indexed by read. """

    @classmethod
    def index_col(cls):
        return 0


# Load by Source (relate/mask/cluster) #################################


class AvgTableLoader(RelTypeTableLoader, ABC):

    @classmethod
    def header_row(cls):
        return 0

    @cached_property
    def data(self):
        # Load the data in the same manner as the superclass.
        data = super().data
        # Rename the column level of the relationships.
        data.columns.rename(REL_NAME, inplace=True)
        return data


class RelTableLoader(AvgTableLoader, ABC):
    pass


class MaskTableLoader(AvgTableLoader, ABC):
    pass


class ClustTableLoader(RelTypeTableLoader, ABC):

    @classmethod
    def header_row(cls):
        return list(range(len(POC_INDEX_NAMES)))

    @cached_property
    def data(self):
        # Load the data in the same manner as the superclass.
        data = super().data
        # Reformat the columns of clusters.
        data.columns = reformat_cluster_index(data.columns)
        return data


# Instantiable Table Loaders ###########################################

class RelPosTableLoader(RelTableLoader, PosTableLoader, RelPosTable):
    """ Load relation data indexed by position. """

    def iter_profiles(self, *args, **kwargs):
        # Relation table loaders have unmasked, unfiltered reads and are
        # thus unsuitable for making RNA profiles. Yield no profiles.
        yield from ()


class RelReadTableLoader(RelTableLoader, ReadTableLoader, RelReadTable):
    """ Load relation data indexed by read. """


class MaskPosTableLoader(MaskTableLoader, PosTableLoader, MaskPosTable):
    """ Load masked bit vector data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section], quantile: float):
        for section in sections:
            yield RnaProfile(path.fill_whitespace(AVERAGE_NAME),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self.process(ratio=True,
                                                 quantile=quantile,
                                                 rels=[MUTAT_REL])[MUTAT_REL])


class MaskReadTableLoader(MaskTableLoader, ReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClustPosTableLoader(ClustTableLoader, PosTableLoader, ClustPosTable):
    """ Load cluster data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section], quantile: float):
        """ Yield RNA mutational profiles from a table. """
        for section, (order, clust) in product(sections, self.ord_clust):
            yield RnaProfile(path.fill_whitespace(fmt_clust_name(order, clust)),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self.process(ratio=True,
                                                 quantile=quantile,
                                                 rels=[MUTAT_REL],
                                                 order=[order],
                                                 cluster=[clust])[MUTAT_REL])


class ClustReadTableLoader(ClustTableLoader, ReadTableLoader, ClustReadTable):
    """ Load cluster data indexed by read. """


class ClustFreqTableLoader(TableLoader, ClustFreqTable):
    """ Load cluster data indexed by cluster. """

    @classmethod
    def index_col(cls):
        return list(range(len(OC_INDEX_NAMES)))

    @classmethod
    def header_row(cls):
        return 0

    @cached_property
    def data(self):
        # Load the data in the same manner as the superclass.
        data = super().data
        # Reformat the index of clusters.
        data.index = reformat_cluster_index(data.index)
        return data


# Helper Functions #####################################################

def find_tables(tables: tuple[str, ...]):
    """ Return a file for each given file/directory of a table. """
    yield from path.find_files_chain(map(Path, tables), [path.TableSeg])


def load(table_file: Path):
    """ Helper function to load a TableLoader from a table file. """
    for loader_type in (RelPosTableLoader, RelReadTableLoader,
                        MaskPosTableLoader, MaskReadTableLoader,
                        ClustPosTableLoader, ClustReadTableLoader,
                        ClustFreqTableLoader):
        try:
            # Try to load the table file using the type of TableLoader.
            return loader_type(table_file)
        except (FileNotFoundError, ValueError):
            # The TableLoader was not the correct type for the file.
            pass
    # None of the TableLoader types were able to open the file.
    raise ValueError(f"Failed to open table: {table_file}")


def reformat_cluster_index(index: pd.MultiIndex):
    """ Ensure that the columns are a MultiIndex and that the order and
    cluster numbers are integers, not strings. """
    return pd.MultiIndex.from_arrays(
        [index.get_level_values(n).astype(int if n in OC_INDEX_NAMES else str)
         for n in index.names],
        names=index.names
    )


def get_clusters(columns: pd.Index | pd.MultiIndex, allow_zero: bool = False):
    """ Return a MultiIndex of non-redundant orders and cluster numbers
    from columns with order and cluster numbers as levels. """
    try:
        return pd.MultiIndex.from_arrays([columns.get_level_values(level)
                                          for level in OC_INDEX_NAMES],
                                         names=OC_INDEX_NAMES).drop_duplicates()
    except KeyError:
        # The index did not contain levels named "order" and "cluster".
        if allow_zero:
            # Default to an index of zero for each level.
            return pd.MultiIndex.from_tuples([(0,) * len(OC_INDEX_NAMES)],
                                             names=OC_INDEX_NAMES)
        # Re-raise the error.
        raise

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
