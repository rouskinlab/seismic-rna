from abc import ABC
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .data import load_relate_dataset
from .report import RelateReport
from ..core import path
from ..core.header import RelHeader, parse_header
from ..core.logs import logger
from ..core.seq import FULL_NAME, DNA, Section
from ..core.table import (Table,
                          Tabulator,
                          DatasetTabulator,
                          PrecountTabulator,
                          PosTable,
                          PosTableWriter,
                          ReadTable,
                          ReadTableWriter,
                          RelTypeTable)


class AvgTable(RelTypeTable, ABC):
    """ Average over an ensemble of RNA structures. """

    @classmethod
    def header_type(cls):
        return RelHeader


class FullTable(Table, ABC):

    @property
    def path_fields(self):
        return {path.TOP: self.top,
                path.CMD: self.kind(),
                path.SAMP: self.sample,
                path.REF: self.ref,
                path.TABLE: self.kind(),
                path.EXT: self.ext()}


class FullPosTable(FullTable, PosTable, ABC):

    @classmethod
    def path_segs(cls):
        return path.REF_DIR_SEGS + (path.PositionTableSeg,)


class FullReadTable(FullTable, ReadTable, ABC):

    @classmethod
    def path_segs(cls):
        return path.REF_DIR_SEGS + (path.ReadTableSeg,)


class RelTable(AvgTable, ABC):

    @classmethod
    def kind(cls):
        return path.CMD_REL_DIR


class RelPosTable(RelTable, FullPosTable, ABC):

    def _iter_profiles(self, *,
                       sections: Iterable[Section] | None,
                       quantile: float,
                       rel: str,
                       k: int | None,
                       clust: int | None):
        # Relation table loaders have unmasked, unfiltered reads and are
        # thus unsuitable for making RNA profiles. Yield no profiles.
        yield from ()


class RelReadTable(RelTable, FullReadTable, ABC):
    pass


class RelPosTableWriter(PosTableWriter, RelPosTable):
    pass


class RelReadTableWriter(ReadTableWriter, RelReadTable):
    pass


class FullTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return 0

    def __init__(self, *, ref: str, refseq: DNA, **kwargs):
        # For the full dataset, the section must be the full reference
        # sequence, and no pattern is used to filter the data.
        super().__init__(refseq=refseq,
                         section=Section(ref, refseq),
                         pattern=None,
                         **kwargs)


class AvgTabulator(Tabulator, ABC):

    def __init__(self, **kwargs):
        # An ensemble average tabulator has no clusters (ks=None).
        super().__init__(ks=None, **kwargs)


class RelateTabulator(FullTabulator, AvgTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [RelPosTableWriter, RelReadTableWriter]


class RelatePrecountTabulator(PrecountTabulator, RelateTabulator):
    pass


class RelateDatasetTabulator(DatasetTabulator, RelateTabulator):
    pass


class TableLoader(Table, ABC):
    """ Load a table from a file. """

    @classmethod
    def find_tables(cls, paths: Iterable[str | Path]):
        """ Yield files of the tables within the given paths. """
        for file in path.find_files_chain(paths, cls.path_segs()):
            if file.name.startswith(cls.kind()):
                yield file

    @classmethod
    def load_tables(cls, paths: Iterable[str | Path]):
        """ Yield tables within the given paths. """
        for file in cls.find_tables(paths):
            try:
                yield cls(file)
            except Exception as error:
                logger.error(error)

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        self._sect = fields.get(path.SECT, FULL_NAME)
        if not self.path.with_suffix(table_file.suffix).samefile(table_file):
            raise ValueError(f"{type(self).__name__} got path {table_file}, "
                             f"but expected {self.path}")

    @property
    def top(self) -> Path:
        return self._out_dir

    @property
    def sample(self) -> str:
        return self._sample

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def sect(self) -> str:
        return self._sect

    @cached_property
    def refseq(self):
        dataset = load_relate_dataset(RelateReport.build_path(
            top=self.top, sample=self.sample, ref=self.ref)
        )
        return dataset.refseq


class RelTypeTableLoader(TableLoader, RelTypeTable, ABC):
    """ Load a table of relationship types. """

    @cached_property
    def data(self) -> pd.DataFrame:
        data = pd.read_csv(self.path,
                           index_col=self.index_cols(),
                           header=self.header_rows())
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.columns)
        # The columns must be replaced with the header index for the
        # type casting to take effect.
        data.columns = header.index
        return data


class PosTableLoader(RelTypeTableLoader, PosTable, ABC):
    """ Load data indexed by position. """


class ReadTableLoader(RelTypeTableLoader, ReadTable, ABC):
    """ Load data indexed by read. """


class RelPosTableLoader(PosTableLoader, RelPosTable):
    """ Load relate data indexed by position. """


class RelReadTableLoader(ReadTableLoader, RelReadTable):
    """ Load relate data indexed by read. """
