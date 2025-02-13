from abc import ABC
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .. import path
from ..header import parse_header
from ..logs import logger
from ..table import (Table,
                     PositionTable,
                     ReadTable,
                     RelTypeTable)


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
    def _attrs(self):
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
