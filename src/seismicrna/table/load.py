from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .base import (MUTAT_REL, CLUST_INDEX_NAMES, REL_NAME,
                   Table, RelTypeTable, PosTable, ReadTable,
                   RelPosTable, RelReadTable,
                   MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustFreqTable)
from ..cluster.names import ENSEMBLE_NAME, ORD_CLS_NAME, fmt_clust_name
from ..core import path
from ..core.rna import RnaProfile
from ..core.sect import Section, INDEX_NAMES


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
            raise ValueError(f"Invalid path for '{self.__class__.__name__}': "
                             f"{table_file} (expected {self.path})")

    @property
    def out_dir(self):
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
        return list(range(len(CLUST_INDEX_NAMES)))

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
            yield RnaProfile(path.fill_whitespace(ENSEMBLE_NAME),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self._ratio_col(MUTAT_REL, quantile))


class MaskReadTableLoader(MaskTableLoader, ReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClustPosTableLoader(ClustTableLoader, PosTableLoader, ClustPosTable):
    """ Load cluster data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section], quantile: float):
        """ Yield RNA mutational profiles from a table. """
        for section in sections:
            for ok, fmut in self._ratio_col(MUTAT_REL, quantile).items():
                yield RnaProfile(path.fill_whitespace(fmt_clust_name(*ok)),
                                 section=section,
                                 sample=self.sample,
                                 data_sect=self.sect,
                                 reacts=fmut)


class ClustReadTableLoader(ClustTableLoader, ReadTableLoader, ClustReadTable):
    """ Load cluster data indexed by read. """


class ClustFreqTableLoader(TableLoader, ClustFreqTable):
    """ Load cluster data indexed by cluster. """

    @classmethod
    def index_col(cls):
        return list(range(len(ORD_CLS_NAME)))

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
    return path.find_files_chain(map(Path, tables), [path.MutTabSeg])


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
        [index.get_level_values(n).astype(int if n in ORD_CLS_NAME else str)
         for n in index.names],
        names=index.names
    )


def get_clusters(columns: pd.Index | pd.MultiIndex, allow_zero: bool = False):
    """ Return a MultiIndex of non-redundant orders and cluster numbers
    from columns with order and cluster numbers as levels. """
    try:
        return pd.MultiIndex.from_arrays([columns.get_level_values(level)
                                          for level in ORD_CLS_NAME],
                                         names=ORD_CLS_NAME).drop_duplicates()
    except KeyError:
        # The index did not contain levels named "order" and "cluster".
        if allow_zero:
            # Default to an index of zero for each level.
            return pd.MultiIndex.from_tuples([(0,) * len(ORD_CLS_NAME)],
                                             names=ORD_CLS_NAME)
        # Re-raise the error.
        raise
