from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..cluster.names import CLS_NAME, ORD_NAME
from ..core import path
from ..core.sect import index_to_pos, index_to_seq

# General fields
READ_TITLE = "Read Name"
R_OBS_TITLE = "Reads Observed"
R_ADJ_TITLE = "Reads Adjusted"
REL_NAME = "Relationship"
POPAVG_TITLE = "pop-avg"
CLUST_INDEX_NAMES = ORD_NAME, CLS_NAME, REL_NAME

# Count relationships
TOTAL_REL = "Called"
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
    'b': TOTAL_REL,
    'n': INFOR_REL,
    'r': MATCH_REL,
    'm': MUTAT_REL,
    'd': DELET_REL,
    'i': INSRT_REL,
    's': SUBST_REL,
    'a': SUB_A_REL,
    'c': SUB_C_REL,
    'g': SUB_G_REL,
    't': SUB_T_REL,
}


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


class RelTypeTable(Table, ABC):
    """ Table with multiple types of relationships. """

    @abstractmethod
    def _count_col(self, col: str):
        """ Get the counts from a column in the data table. """
        return self.data[col]

    def _count_rel(self, rel: str):
        """ Count the bits for a relationship. """
        if rel == INFOR_REL:
            # Sum the matched and mutated bits.
            return self._count_col(MATCH_REL) + self._count_col(MUTAT_REL)
        # The relationship should appear in the table, so return it.
        return self._count_col(rel)

    @abstractmethod
    def _fract_rel(self, rel: str):
        """ Compute the fraction for a relationship. """
        # Determine the relationship to use as the denominator.
        denom = TOTAL_REL if rel == INFOR_REL else INFOR_REL
        # Compute the ratio of the numerator and the denominator.
        return self._count_rel(rel) / self._count_rel(denom)

    def count_rel(self, code: str):
        """ Count the bits for a relationship given its code. """
        return self._count_rel(REL_CODES[code])

    def fract_rel(self, code: str):
        """ Compute the fraction for a relationship given its code. """
        return self._fract_rel(REL_CODES[code])


# Table by Source (relate/mask/cluster) ################################

class RelTable(RelTypeTable, ABC):

    @cache
    def _count_col(self, col: str):
        return super()._count_col(col).rename(col)

    @cache
    def _fract_rel(self, rel: str):
        return super()._fract_rel(rel).rename(rel)


class MaskTable(RelTypeTable, ABC):

    @cache
    def _count_col(self, col: str):
        return super()._count_col(col).rename(col)

    @cache
    def _fract_rel(self, rel: str):
        return super()._fract_rel(rel).rename(rel)


class ClustTable(RelTypeTable, ABC):

    @cached_property
    def clusters(self):
        """ MultiIndex of all clusters. """
        return self.data.columns.droplevel(-1).drop_duplicates()

    @cache
    def _count_col(self, col: str):
        # Make a slicer that pulls the given column from all clusters.
        slicer = (slice(None),) * self.clusters.nlevels + (col,)
        # Copy the given column from all clusters into a new DataFrame.
        counts = self.data.loc[:, slicer].copy()
        # Drop the last level (the given column) from the column index.
        counts.columns = counts.columns.droplevel(-1)
        return counts

    @cache
    def _fract_rel(self, rel: str):
        return super()._fract_rel(rel)


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
