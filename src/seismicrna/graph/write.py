from abc import ABC, abstractmethod
from logging import getLogger
from functools import cache, cached_property
from pathlib import Path
from typing import Iterable

from .base import GraphBase
from ..table.load import load

logger = getLogger(__name__)


class GraphWriter(ABC):
    """ Write the proper types of graphs for a given table(s). """

    def __init__(self, table_files: Iterable[Path]):
        self.table_files = list(table_files)

    @abstractmethod
    def iter(self, *_, **__):
        """ Yield every graph for the table. """
        yield GraphBase()

    def write(self, *args, csv: bool, html: bool, pdf: bool, **kwargs):
        """ Generate and write every graph for the table. """
        # Get the paths for every graph.
        paths = list()
        for graph in self.iter(*args, **kwargs):
            paths.extend(graph.write(csv=csv, html=html, pdf=pdf))
        return paths


class OneTableGraphWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for a given table. """

    def __init__(self, table_file: Path):
        super().__init__([table_file])

    @property
    def table_file(self) -> Path:
        try:
            table_file, = self.table_files
        except ValueError:
            raise ValueError(f"{self.__class__.__name__} requires exactly 1 "
                             f"table file, but got {len(self.table_files)}")
        return table_file

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return load(self.table_file)


class TwoTableGraphWriter(GraphWriter, ABC):
    """ Write the proper types of graphs for two given tables. """

    def __init__(self, table1_file: Path, table2_file: Path):
        super().__init__([table1_file, table2_file])

    @cache
    def _get_table_files(self) -> tuple[Path, Path]:
        try:
            table1_file, table2_file = self.table_files
        except ValueError:
            raise ValueError(f"{self.__class__.__name__} requires exactly 2 "
                             f"table file, but got {len(self.table_files)}")
        return table1_file, table2_file

    @property
    def table1_file(self) -> Path:
        table1_file, _ = self._get_table_files()
        return table1_file

    @property
    def table2_file(self) -> Path:
        _, table2_file = self._get_table_files()
        return table2_file

    @cached_property
    def table1(self):
        """ The first table providing the data for the graph(s). """
        return load(self.table1_file)

    @cached_property
    def table2(self):
        """ The second table providing the data for the graph(s). """
        return load(self.table2_file)
