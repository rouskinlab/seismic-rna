from abc import ABC
from functools import cache, cached_property
from typing import Iterable

from .dataset import load_relate_dataset
from .io import RelateFile
from ..core.header import RelHeader
from ..core.seq import DNA, Region
from ..core.table import (Tabulator,
                          BatchTabulator,
                          CountTabulator,
                          DatasetTabulator,
                          PositionTable,
                          PositionTableLoader,
                          PositionTableWriter,
                          ReadTable,
                          ReadTableLoader,
                          ReadTableWriter,
                          RelTypeTable)


class AverageTable(RelTypeTable, ABC):
    """ Average over an ensemble of RNA structures. """

    @classmethod
    def get_header_type(cls):
        return RelHeader


class RelateTable(AverageTable, RelateFile, ABC):

    @classmethod
    def get_load_function(cls):
        return load_relate_dataset


class RelatePositionTable(RelateTable, PositionTable, ABC):

    def _iter_profiles(self, *,
                       regions: Iterable[Region] | None,
                       quantile: float,
                       rel: str,
                       k: int | None,
                       clust: int | None):
        # Relate tables have unfiltered reads and are thus unsuitable
        # for making RNA profiles: do not generate any.
        yield from ()


class RelateReadTable(RelateTable, ReadTable, ABC):
    pass


class RelatePositionTableLoader(PositionTableLoader, RelatePositionTable):
    """ Load relate data indexed by position. """


class RelateReadTableLoader(ReadTableLoader, RelateReadTable):
    """ Load relate data indexed by read. """


class RelatePositionTableWriter(PositionTableWriter, RelatePositionTable):
    pass


class RelateReadTableWriter(ReadTableWriter, RelateReadTable):
    pass


class FullTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return 0

    def __init__(self, *,
                 ref: str,
                 refseq: DNA,
                 count_ends: bool = False,
                 **kwargs):
        # For a full tabulator, the full reference sequence must be used
        # as the region.
        super().__init__(region=Region(ref, refseq),
                         count_ends=count_ends,
                         **kwargs)


class AverageTabulator(Tabulator, ABC):

    @cached_property
    def data_per_clust(self):
        # An ensemble average tabulator has no per-cluster data.
        return None


class RelateTabulator(FullTabulator, AverageTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [RelatePositionTableWriter, RelateReadTableWriter]


class RelateCountTabulator(CountTabulator, RelateTabulator):
    pass


class RelateBatchTabulator(BatchTabulator, RelateTabulator):
    pass


class RelateDatasetTabulator(DatasetTabulator, RelateTabulator):

    @classmethod
    @cache
    def init_kws(cls):
        return super().init_kws() + cls._list_args(FullTabulator.__init__)
