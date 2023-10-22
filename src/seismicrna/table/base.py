from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..core import path
from ..core.batch import CLUST_NAME, ORDER_NAME, REL_NAME
from ..core.mu import winsorize
from ..core.seq import index_to_pos, index_to_seq

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

# Keyword argument for each level of the column MultiIndex.
LEVEL_KEYS = {REL_NAME: "rels", ORDER_NAME: "order", CLUST_NAME: "cluster"}


# Table Base Classes ###################################################

class Table(ABC):
    """ Table base class. """

    @property
    @abstractmethod
    def top(self):
        """ Path of the table's output directory. """
        return Path()

    @property
    @abstractmethod
    def sample(self):
        """ Name of the table's sample. """
        return ""

    @property
    @abstractmethod
    def ref(self):
        """ Name of the table's reference. """
        return ""

    @property
    @abstractmethod
    def sect(self):
        """ Name of the table's section. """
        return ""

    @cached_property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        """ Table's data frame. """
        return pd.DataFrame()

    @classmethod
    @abstractmethod
    def kind(cls) -> str:
        """ Kind of table. """
        return ""

    @classmethod
    @abstractmethod
    def by_read(cls):
        """ Whether the table contains data for each read. """
        return False

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

    @cached_property
    def multiindexed(self):
        """ Whether the columns are a MultiIndex. """
        return isinstance(self.data.columns, pd.MultiIndex)

    @cached_property
    def rel_level(self):
        """ Index of the column level indicating the relationship. """
        return self.data.columns.names.index(REL_NAME)

    def _switch_rel(self, column: str | tuple, new_rel: str):
        """ Switch the relationship in a column label. """
        if isinstance(column, str):
            column = column,
        return (
            column[: self.rel_level] + (new_rel,) + column[self.rel_level + 1:]
            if self.multiindexed else new_rel
        )

    @cache
    def get_ratio(self,
                  column: tuple,
                  quantile: float,
                  precision: int | None = None):
        """ Compute the ratio for a column. """
        # Determine the relationship to use as the numerator.
        numer_rel = column[self.rel_level]
        # Determine the relationship to use as the denominator.
        denom_rel = COVER_REL if numer_rel == INFOR_REL else INFOR_REL
        # Determine the column to use as the denominator.
        denom_col = self._switch_rel(column, denom_rel)
        # Compute the ratio of the numerator and the denominator.
        ratio = self.data[column] / self.data[denom_col]
        # If a quantile was given, then winsorize to it.
        if quantile is not None:
            ratio = winsorize(ratio, quantile)
        # Round the ratio to the desired precision.
        if precision is not None:
            ratio = ratio.round(precision)
        return ratio

    def _get_indexer(self, **kwargs):
        """ Format a column selection into a column indexer. """
        return (
            tuple(kwargs.get(LEVEL_KEYS[level_name], slice(None))
                  for level_name in self.data.columns.names)
            if self.multiindexed
            else kwargs.get(LEVEL_KEYS[self.data.columns.name], slice(None))
        )

    def process(self, *,
                ratio: bool,
                quantile: float = 0.,
                precision: int | None = None,
                **kwargs: list):
        """ Select, process, and return data from the table. """
        # Instantiate an empty DataFrame with the index and columns.
        columns = self.data.loc[:, self._get_indexer(**kwargs)].columns
        data = pd.DataFrame(index=self.data.index, columns=columns, dtype=float)
        # Fill in the DataFrame one column at a time.
        for column in columns:
            data[column] = (self.get_ratio(column, quantile, precision) if ratio
                            else self.data[column]).values
        return data


# Table by Source (relate/mask/cluster) ################################

class RelTable(RelTypeTable, ABC):
    pass


class MaskTable(RelTypeTable, ABC):
    pass


class ClustTable(RelTypeTable, ABC):

    @cached_property
    def ord_clust(self):
        """ MultiIndex of all order-cluster pairs. """
        return self.data.columns.droplevel(REL_NAME).drop_duplicates()

    @cached_property
    def orders(self):
        """ Index of all order numbers. """
        return self.data.columns.get_level_values(ORDER_NAME).drop_duplicates()

    def get_clusters(self, orders: list):
        """ Index of cluster numbers for the given orders. """
        data = self.data.loc[:, orders]
        return data.columns.get_level_values(CLUST_NAME).drop_duplicates()


# Table by Index (position/read/frequency) #############################

class PosTable(RelTypeTable, ABC):
    """ Table indexed by position. """

    @classmethod
    def by_read(cls):
        return False

    @cached_property
    def seq(self):
        return index_to_seq(self.data.index)

    @property
    def positions(self):
        return index_to_pos(self.data.index)

    @property
    def end5(self):
        return int(self.positions[0])

    @property
    def end3(self):
        return int(self.positions[-1])


class ReadTable(RelTypeTable, ABC):
    """ Table indexed by read. """

    @classmethod
    def by_read(cls):
        return True

    @property
    def reads(self):
        return self.data.index.values


# Table by Source and Index ############################################

class RelPosTable(RelTable, PosTable, ABC):

    @classmethod
    def kind(cls):
        return path.RELATE_POS_TAB


class MaskPosTable(MaskTable, PosTable, ABC):

    @classmethod
    def kind(cls):
        return path.MASKED_POS_TAB


class ClustPosTable(ClustTable, PosTable, ABC):

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

    @property
    def clusters(self):
        return self.data.index.values

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
