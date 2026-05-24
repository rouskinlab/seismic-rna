from abc import ABC
from functools import cache, cached_property
from typing import Iterable

from .dataset import load_idmut_dataset
from .io import IDmutFile
from ..core.header import RelHeader
from ..core.seq import DNA, Region
from ..core.table import (
    Tabulator,
    BatchTabulator,
    CountTabulator,
    DatasetTabulator,
    PositionTable,
    PositionTableLoader,
    PositionTableWriter,
    ReadTable,
    ReadTableLoader,
    ReadTableWriter,
    RelTypeTable,
)


class AverageTable(RelTypeTable, ABC):
    """Average over an ensemble of RNA structures."""

    @classmethod
    def get_header_type(cls):
        return RelHeader


class IDmutTable(AverageTable, IDmutFile, ABC):
    @classmethod
    def get_load_function(cls):
        return load_idmut_dataset


class IDmutPositionTable(IDmutTable, PositionTable, ABC):
    def _iter_profiles(
        self,
        *,
        regions: Iterable[Region] | None,
        quantile: float,
        rel: str,
        k: int | None,
        clust: int | None,
    ):
        """
        Yield RNA profiles from this table (always empty for IDmut).

        IDmut tables contain unfiltered reads and are unsuitable for
        generating RNA mutational profiles.

        Parameters
        ----------
        regions: Iterable[Region] | None
            Regions for which to generate profiles; unused.
        quantile: float
            Quantile for normalization; unused.
        rel: str
            Relationship type to profile; unused.
        k: int | None
            Number of clusters; unused.
        clust: int | None
            Cluster index; unused.
        """
        # IDmut tables have unfiltered reads and are thus unsuitable
        # for making RNA profiles: do not generate any.
        yield from ()


class IDmutReadTable(IDmutTable, ReadTable, ABC):
    pass


class IDmutPositionTableLoader(PositionTableLoader, IDmutPositionTable):
    """Load IDmut data indexed by position."""


class IDmutReadTableLoader(ReadTableLoader, IDmutReadTable):
    """Load IDmut data indexed by read."""


class IDmutPositionTableWriter(PositionTableWriter, IDmutPositionTable):
    pass


class IDmutReadTableWriter(ReadTableWriter, IDmutReadTable):
    pass


class FullTabulator(Tabulator, ABC):
    @classmethod
    def get_null_value(cls):
        return 0

    def __init__(self, *, ref: str, refseq: DNA, count_ends: bool = False, **kwargs):
        """
        Initialize a full-reference tabulator.

        Parameters
        ----------
        ref: str
            Name of the reference sequence.
        refseq: DNA
            Full reference sequence; used to construct the region
            spanning the entire reference.
        count_ends: bool
            Whether to count read end positions.
        **kwargs
            Additional keyword arguments forwarded to the parent class.
        """
        # For a full tabulator, the full reference sequence must be used
        # as the region.
        super().__init__(region=Region(ref, refseq), count_ends=count_ends, **kwargs)


class AverageTabulator(Tabulator, ABC):
    @cached_property
    def data_per_clust(self):
        # An ensemble average tabulator has no per-cluster data.
        return None


class IDmutTabulator(FullTabulator, AverageTabulator, ABC):
    @classmethod
    def table_types(cls):
        return [IDmutPositionTableWriter, IDmutReadTableWriter]


class IDmutCountTabulator(CountTabulator, IDmutTabulator):
    pass


class IDmutBatchTabulator(BatchTabulator, IDmutTabulator):
    pass


class IDmutDatasetTabulator(DatasetTabulator, IDmutTabulator):
    @classmethod
    @cache
    def init_kws(cls):
        return super().init_kws() + cls._list_args(FullTabulator.__init__)
