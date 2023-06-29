from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .base import (MUTAT_REL, POPAVG_TITLE, CLUST_INDEX_NAMES,
                   Table, PosTable, ReadTable,
                   RelPosTable, RelReadTable,
                   MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustFreqTable)
from ..cluster.names import ORD_CLS_NAME
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
    def data(self):
        return pd.read_csv(self.path,
                           index_col=self.index_col(),
                           header=self.header_row())


# Load by Index (position/read/frequency) ##############################

class PosTableLoader(TableLoader, PosTable, ABC):
    """ Load data indexed by position. """

    @classmethod
    def index_col(cls):
        return tuple(range(len(INDEX_NAMES)))

    @abstractmethod
    def iter_profiles(self, sections: Iterable[Section]):
        """ Yield RNA mutational profiles from the table. """
        for section in sections:
            yield RnaProfile("", section, "", "", pd.Series())


class ReadTableLoader(TableLoader, ReadTable, ABC):
    """ Load data indexed by read. """

    @classmethod
    def index_col(cls):
        return [0]


# Load by Source (relate/mask/cluster) #################################

class RelTableLoader(TableLoader, ABC):

    @classmethod
    def header_row(cls):
        return [0]


class MaskTableLoader(TableLoader, ABC):

    @classmethod
    def header_row(cls):
        return [0]


class ClustTableLoader(TableLoader, ABC):

    @classmethod
    def header_row(cls):
        return list(range(len(CLUST_INDEX_NAMES)))


# Instantiable Table Loaders ###########################################

class RelPosTableLoader(RelTableLoader, PosTableLoader, RelPosTable):
    """ Load relation data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        # Relation table loaders have unmasked, unfiltered reads and are
        # thus unsuitable for making RNA profiles. Yield no profiles.
        yield from ()


class RelReadTableLoader(RelTableLoader, ReadTableLoader, RelReadTable):
    """ Load relation data indexed by read. """


class MaskPosTableLoader(MaskTableLoader, PosTableLoader, MaskPosTable):
    """ Load masked bit vector data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        for section in sections:
            yield RnaProfile(title=path.fill_whitespace(POPAVG_TITLE),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self._fract_rel(MUTAT_REL))


class MaskReadTableLoader(MaskTableLoader, ReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClustPosTableLoader(ClustTableLoader, PosTableLoader, ClustPosTable):
    """ Load cluster data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        """ Yield RNA mutational profiles from a table. """
        for section in sections:
            for (order, k), fmut in self._fract_rel(MUTAT_REL).items():
                yield RnaProfile(f"Cluster_{order}-{k}",
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
        return [0]


# Helper Functions #####################################################

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
