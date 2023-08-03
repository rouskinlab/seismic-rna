from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..cluster.names import CLS_NAME, ORD_NAME
from ..core import path
from ..core.mu import winsorize
from ..core.sect import index_to_pos, index_to_seq

# General fields
READ_TITLE = "Read Name"
R_OBS_TITLE = "Reads Observed"
R_ADJ_TITLE = "Reads Adjusted"
REL_NAME = "Relationship"
CLUST_INDEX_NAMES = ORD_NAME, CLS_NAME, REL_NAME

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


# Table Base Classes ###################################################

class Table(ABC):
    """ Table base class. """

    @property
    @abstractmethod
    def out_dir(self):
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
        return (path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                path.MutTabSeg)

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
        return {path.TOP: self.out_dir,
                path.MOD: path.MOD_TABLE,
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
        return f"{self.__class__.__name__} at {self.path}"


class RelTypeTable(Table, ABC):
    """ Table with multiple types of relationships. """

    @property
    def _rel_level_index(self):
        """ Index of the column level indicating the relationship. """
        return self.data.columns.names.index(REL_NAME)

    def _switch_rel(self, column: tuple, new_rel: str):
        """ Switch the relationship in a column label. """
        return new_rel,

    @cache
    def _count_col(self, column: tuple):
        """ Count the bits for a column. """
        return self.data.loc[:, column].squeeze()

    @cache
    def _ratio_col(self, column: tuple,
                   quantile: float | None = None,
                   precision: int | None = None):
        """ Compute the ratio for a column. """
        # Determine the relationship to use as the numerator.
        numer_rel = column[self._rel_level_index]
        # Determine the relationship to use as the denominator.
        denom_rel = COVER_REL if numer_rel == INFOR_REL else INFOR_REL
        # Determine the column to use as the denominator.
        denom_col = self._switch_rel(column, denom_rel)
        # Compute the ratio of the numerator and the denominator.
        ratio = self._count_col(column) / self._count_col(denom_col)
        # If a quantile was given, then winsorize to it.
        if quantile is not None:
            ratio = winsorize(ratio, quantile)
        # Round the ratio to the desired precision.
        if precision is not None:
            ratio = ratio.round(precision)
        return ratio

    @staticmethod
    def _format_selection(**kwargs) -> dict[str, list]:
        """ Format keyword arguments into a valid column selection. """
        selection = dict()
        for key, level in dict(order=ORD_NAME,
                               cluster=CLS_NAME,
                               rels=REL_NAME).items():
            try:
                selection[level] = kwargs[key]
            except KeyError:
                pass
        return selection

    def _get_indexer(self, selection: dict[str, list]):
        """ Format a column selection into a column indexer. """
        return selection.get(REL_NAME, slice(None))

    def select(self, ratio: bool,
               quantile: float = -1.,
               precision: int | None = None,
               **kwargs: list):
        """ Output selected data from the table as a DataFrame. """
        # Instantiate an empty DataFrame with the index and columns.
        indexer = self._get_indexer(self._format_selection(**kwargs))
        columns = self.data.loc[:, indexer].columns
        data = pd.DataFrame(index=self.data.index, columns=columns, dtype=float)
        # Fill in the DataFrame one column at a time.
        for column in columns:
            data.loc[:, column] = (self._ratio_col(column, quantile, precision)
                                   if ratio else self._count_col(column))
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
        return self.data.columns.get_level_values(ORD_NAME).drop_duplicates()

    def get_clusters(self, orders: list):
        """ Index of cluster numbers for the given orders. """
        data = self.data.loc[:, orders]
        return data.columns.get_level_values(CLS_NAME).drop_duplicates()

    def _switch_rel(self, column: tuple, new_rel: str):
        return (column[: self._rel_level_index]
                + (new_rel,)
                + column[self._rel_level_index + 1:])

    def _get_indexer(self, selection: dict[str, list]):
        return tuple(selection.get(level, slice(None))
                     for level in self.data.columns.names)


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
