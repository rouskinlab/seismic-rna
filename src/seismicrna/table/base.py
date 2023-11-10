from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any, Generator, Iterable

import pandas as pd

from ..core import path
from ..core.batch import RB_INDEX_NAMES
from ..core.header import (REL_NAME,
                           Header,
                           RelHeader,
                           ClustHeader,
                           RelClustHeader,
                           format_clust_name,
                           parse_header)
from ..core.mu import winsorize
from ..core.rna import RnaProfile
from ..core.seq import SEQ_INDEX_NAMES, Section, index_to_pos, index_to_seq

# General fields
READ_TITLE = "Read Name"
R_OBS_TITLE = "Reads Observed"
R_ADJ_TITLE = "Reads Adjusted"

# Count relationships
COVER_REL = "Covered"
INFOR_REL = "Informed"
MATCH_REL = "Matched"
MUTAT_REL = "Mutated"
DELET_REL = "Deleted"
INSRT_REL = "Inserted"
SUBST_REL = "Subbed"
SUB_A_REL = "Subbed-A"
SUB_C_REL = "Subbed-C"
SUB_G_REL = "Subbed-G"
SUB_T_REL = "Subbed-T"

# One-letter codes for each type of relationship
REL_CODES = {
    'v': COVER_REL,
    'n': INFOR_REL,
    'r': MATCH_REL,
    'm': MUTAT_REL,
    's': SUBST_REL,
    'a': SUB_A_REL,
    'c': SUB_C_REL,
    'g': SUB_G_REL,
    't': SUB_T_REL,
    'd': DELET_REL,
    'i': INSRT_REL,
}

# Columns of each relation-based table
TABLE_RELS = list(REL_CODES.values())


def get_rel_name(rel_code: str):
    """ Get the name of a relationship from its code. """
    try:
        return REL_CODES[rel_code]
    except KeyError:
        raise ValueError(f"Relationship code must be one of {list(REL_CODES)}, "
                         f"but got {repr(rel_code)}")


def _get_denom_rel(rel: str):
    """ Get the relationship that serves as the denominator. """
    return COVER_REL if rel == COVER_REL or rel == INFOR_REL else INFOR_REL


def _get_denom_cols(numer_cols: pd.Index):
    """ Get the denominator columns based on the numerator columns. """
    if isinstance(numer_cols, pd.MultiIndex):
        return pd.MultiIndex.from_arrays(
            [list(map(_get_denom_rel, numer_cols.get_level_values(name)))
             if name == REL_NAME
             else numer_cols.get_level_values(name)
             for name in numer_cols.names],
            names=numer_cols.names
        )
    if numer_cols.name is not None and numer_cols.name != REL_NAME:
        raise ValueError(f"Expected index to be named {repr(REL_NAME)}, "
                         f"but got {repr(numer_cols.name)}")
    return pd.Index(list(map(_get_denom_rel, numer_cols)),
                    name=numer_cols.name)


# Table Base Classes ###################################################

class Table(ABC):
    """ Table base class. """

    @classmethod
    @abstractmethod
    def kind(cls) -> str:
        """ Kind of table. """

    @classmethod
    @abstractmethod
    def by_read(cls) -> bool:
        """ Whether the table contains data for each read. """

    @classmethod
    @abstractmethod
    def transposed(cls) -> bool:
        """ Whether the data is saved as its transpose. """

    @classmethod
    @abstractmethod
    def header_type(cls) -> type[Header]:
        """ Type of the header for the table. """

    @classmethod
    def header_depth(cls):
        return cls.header_type().num_levels()

    @classmethod
    def header_rows(cls) -> list[int]:
        """ Row(s) of the file to use as the columns. """
        return list(range(cls.header_depth()))

    @classmethod
    @abstractmethod
    def index_depth(cls) -> int:
        """ Number of columns in the index. """

    @classmethod
    def index_cols(cls) -> list[int]:
        """ Column(s) of the file to use as the index. """
        return list(range(cls.index_depth()))

    @classmethod
    def path_segs(cls):
        """ Table's path segments. """
        return path.TABLE_SEGS

    @classmethod
    def gzipped(cls):
        """ Whether the table's file is compressed with gzip. """
        return cls.by_read()

    @classmethod
    def ext(cls):
        """ Table's file extension: either '.csv' or '.csv.gz'. """
        return path.CSVZIP_EXT if cls.gzipped() else path.CSV_EXT

    @property
    @abstractmethod
    def top(self) -> Path:
        """ Path of the table's output directory. """

    @property
    @abstractmethod
    def sample(self) -> str:
        """ Name of the table's sample. """

    @property
    @abstractmethod
    def ref(self) -> str:
        """ Name of the table's reference. """

    @property
    @abstractmethod
    def sect(self) -> str:
        """ Name of the table's section. """

    @property
    @abstractmethod
    def _data(self) -> pd.DataFrame:
        """ Table's raw data frame. """

    @cached_property
    def data(self):
        """ Table's data frame. """
        return self._data

    @cached_property
    def header(self):
        """ Header for the table's data. """
        header = parse_header(self.data.columns)
        if not isinstance(header, self.header_type()):
            raise TypeError(f"Expected {self.header_type().__name__}, "
                            f"but got {type(header).__name__}")
        return header

    @property
    def path_fields(self) -> dict[str, Any]:
        """ Table's path fields. """
        return {path.TOP: self.top,
                path.CMD: path.CMD_TBL_DIR,
                path.SAMP: self.sample,
                path.REF: self.ref,
                path.SECT: self.sect,
                path.TABLE: self.kind(),
                path.EXT: self.ext()}

    @cached_property
    def path(self):
        """ Path of the table's CSV file (possibly gzipped). """
        return path.buildpar(*self.path_segs(), **self.path_fields)

    def __str__(self):
        return f"{type(self).__name__} at {self.path}"


class RelTypeTable(Table, ABC):
    """ Table with multiple types of relationships. """

    @classmethod
    def transposed(cls):
        return False

    def fetch_count(self, **kwargs) -> pd.Series | pd.DataFrame:
        """ Fetch counts of one or more columns. """
        return self.data.loc[:, self.header.select(**kwargs)]

    def fetch_ratio(self, *,
                    quantile: float = 0.,
                    precision: int | None = None,
                    **kwargs) -> pd.Series | pd.DataFrame:
        """ Fetch ratios of one or more columns. """
        # Fetch the data for the numerator.
        numer = self.fetch_count(**kwargs)
        # Fetch the data for the denominator.
        denom = self.data.loc[:, _get_denom_cols(numer.columns)]
        # Compute the ratio of the numerator and the denominator.
        ratio = winsorize(numer / denom.values, quantile)
        # Round the ratio to the desired precision.
        return ratio.round(precision) if precision is not None else ratio


# Table by Source (relate/mask/cluster) ################################

class AvgTable(RelTypeTable, ABC):
    """ Average over an ensemble of RNA structures. """

    @classmethod
    def header_type(cls):
        return RelHeader


class RelTable(AvgTable, ABC):
    pass


class MaskTable(AvgTable, ABC):
    pass


class ClustTable(RelTypeTable, ABC):
    """ Cluster for each RNA structure in an ensemble. """

    @classmethod
    def header_type(cls):
        return RelClustHeader


# Table by Index (position/read/frequency) #############################

class PosTable(RelTypeTable, ABC):
    """ Table indexed by position. """

    MASK = "pos-mask"

    @classmethod
    def by_read(cls):
        return False

    @classmethod
    def index_depth(cls):
        return len(SEQ_INDEX_NAMES)

    @cached_property
    def range(self):
        return self._data.index

    @cached_property
    def range_int(self):
        return index_to_pos(self.range)

    @cached_property
    def seq(self):
        return index_to_seq(self.range)

    @property
    def end5(self):
        return int(self.range_int[0])

    @property
    def end3(self):
        return int(self.range_int[-1])

    @cached_property
    def masked_bool(self):
        return self._data.isna().all(axis=1)

    @cached_property
    def unmasked_bool(self):
        return ~self.masked_bool

    @cached_property
    def unmasked(self):
        return self._data.index[self.unmasked_bool]

    @cached_property
    def unmasked_int(self):
        return index_to_pos(self.unmasked)

    @cached_property
    def section(self):
        """ Section covered by the table. """
        section = Section(self.ref,
                          self.seq,
                          seq5=self.end5,
                          end5=self.end5,
                          end3=self.end3,
                          name=self.sect)
        section.add_mask(self.MASK, self.unmasked_int, invert=True)
        return section

    @cached_property
    def data(self):
        return self._data.loc[self.unmasked]

    @abstractmethod
    def _iter_profiles(self,
                       sections: Iterable[Section],
                       quantile: float,
                       order: int | None,
                       clust: int | None) -> Generator[RnaProfile, Any, Any]:
        """ Yield RNA mutational profiles from the table. """

    def iter_profiles(self,
                      sections: Iterable[Section] | None = None,
                      quantile: float = 0.,
                      order: int | None = None,
                      clust: int | None = None):
        """ Yield RNA mutational profiles from the table. """
        yield from self._iter_profiles((sections if sections is not None
                                        else [self.section]),
                                       quantile,
                                       order,
                                       clust)


class ReadTable(RelTypeTable, ABC):
    """ Table indexed by read. """

    @classmethod
    def by_read(cls):
        return True

    @classmethod
    def index_depth(cls):
        return len(RB_INDEX_NAMES)

    @property
    def reads(self):
        return self._data.index.values


# Table by Source and Index ############################################

class RelPosTable(RelTable, PosTable, ABC):

    @classmethod
    def kind(cls):
        return path.RELATE_POS_TAB

    def _iter_profiles(self,
                       sections: Iterable[Section],
                       quantile: float,
                       order: int | None,
                       clust: int | None):
        # Relation table loaders have unmasked, unfiltered reads and are
        # thus unsuitable for making RNA profiles. Yield no profiles.
        yield from ()


class ProfilePosTable(PosTable, ABC):

    def _iter_profiles(self,
                       sections: Iterable[Section],
                       quantile: float,
                       order: int | None,
                       clust: int | None):
        """ Yield RNA mutational profiles from a table. """
        for section in sections if sections is not None else [self.section]:
            for ho, hc in self.header.clusts:
                if (not order or order == ho) and (not clust or clust == hc):
                    name = format_clust_name(ho, hc, allow_zero=True)
                    yield RnaProfile(title=path.fill_whitespace(name),
                                     section=section,
                                     sample=self.sample,
                                     data_sect=self.sect,
                                     data=self.fetch_ratio(quantile=quantile,
                                                           rel=MUTAT_REL,
                                                           order=ho,
                                                           clust=hc))


class MaskPosTable(MaskTable, ProfilePosTable, ABC):

    @classmethod
    def kind(cls):
        return path.MASKED_POS_TAB


class ClustPosTable(ClustTable, ProfilePosTable, ABC):

    @classmethod
    def kind(cls):
        return path.CLUST_POS_TAB


class RelReadTable(RelTable, ReadTable, ABC):

    @classmethod
    def kind(cls):
        return path.RELATE_READ_TAB


class MaskReadTable(MaskTable, ReadTable, ABC):

    @classmethod
    def kind(cls):
        return path.MASKED_READ_TAB


class ClustReadTable(ClustTable, ReadTable, ABC):

    @classmethod
    def kind(cls):
        return path.CLUST_READ_TAB


class ClustFreqTable(Table, ABC):

    @classmethod
    def kind(cls):
        return path.CLUST_FREQ_TAB

    @classmethod
    def by_read(cls):
        return False

    @classmethod
    def transposed(cls):
        return True

    @classmethod
    def header_type(cls):
        return ClustHeader

    @classmethod
    def index_depth(cls):
        return 1

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
