from abc import ABC
from functools import cache, cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .dataset import load_relate_dataset
from .io import RelateFile
from ..core import path
from ..core.header import RelHeader, parse_header
from ..core.logs import logger
from ..core.seq import DNA, Region
from ..core.table import (Table,
                          Tabulator,
                          BatchTabulator,
                          CountTabulator,
                          DatasetTabulator,
                          PositionTable,
                          PositionTableWriter,
                          ReadTable,
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


class TableLoader(Table, ABC):
    """ Load a table from a file. """

    @classmethod
    def find_tables(cls, paths: Iterable[str | Path]):
        """ Yield files of the tables within the given paths. """
        for file in path.find_files_chain(paths, cls.get_path_seg_types()):
            if file.name.startswith(cls.get_step()):
                yield file

    @classmethod
    def load_tables(cls, paths: Iterable[str | Path], **kwargs):
        """ Yield tables within the given paths. """
        for file in cls.find_tables(paths):
            try:
                yield cls(file, **kwargs)
            except Exception as error:
                logger.error(error)

    def __init__(self, table_file: str | Path, **kwargs):
        load_function = self.get_load_function()
        report_file = path.cast_path(table_file,
                                     self.get_path_seg_types(),
                                     load_function.report_path_seg_types,
                                     load_function.report_path_auto_fields)
        self._dataset = load_function(report_file, **kwargs)

    @property
    def _source(self):
        return self._dataset


class RelTypeTableLoader(TableLoader, RelTypeTable, ABC):
    """ Load a table of relationship types. """

    @cached_property
    def data(self) -> pd.DataFrame:
        data = pd.read_csv(self.path,
                           index_col=self.get_index_cols(),
                           header=self.get_header_rows())
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.columns)
        # The columns must be replaced with the header index for the
        # type casting to take effect.
        data.columns = header.index
        return data


class PositionTableLoader(RelTypeTableLoader, PositionTable, ABC):
    """ Load data indexed by position. """


class ReadTableLoader(RelTypeTableLoader, ReadTable, ABC):
    """ Load data indexed by read. """


class RelatePositionTableLoader(PositionTableLoader, RelatePositionTable):
    """ Load relate data indexed by position. """


class RelateReadTableLoader(ReadTableLoader, RelateReadTable):
    """ Load relate data indexed by read. """
