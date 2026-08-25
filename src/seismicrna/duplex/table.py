from __future__ import annotations
from abc import ABC
from functools import cached_property

from .dataset import load_duplex_dataset
from .io import DuplexFile
from ..core.header import parse_header
from ..core.table.load import PositionTableLoader
from ..filter.table import PartialPositionTable
from ..idmut.table import AverageTable


class DuplexTable(AverageTable, DuplexFile, ABC):
    @classmethod
    def get_load_function(cls):
        return load_duplex_dataset


class DuplexPositionTable(DuplexTable, PartialPositionTable, ABC):
    """Position table of a duplex of two sources into one duplex."""

    @property
    def cut(self) -> int:
        """Length of the 5' strand (position of the strand break)."""
        return self._attrs.cut

    @property
    def table1(self):
        """Position table of the 5' strand's source."""
        return self._attrs.table1

    @property
    def table2(self):
        """Position table of the 3' strand's source, or None if the 3'
        strand has no data (a bare sequence)."""
        return self._attrs.table2

    @property
    def data1(self):
        """Per-position data of the 5' strand's source."""
        return self.table1.data

    @property
    def data2(self):
        """Per-position data of the 3' strand's source, or None."""
        return self.table2.data if self.table2 is not None else None


class DuplexPositionTableLoader(PositionTableLoader, DuplexPositionTable):
    @cached_property
    def data(self):
        # A duplex CSV has either 1 header row (unclustered: relationship
        # only) or 3 (a cluster product: relationship, K, cluster). The
        # number is recorded in the duplex report, so read that many rows
        # rather than the fixed count implied by the class's header type.
        import pandas as pd

        data = pd.read_csv(
            self.path,
            index_col=self.get_index_cols(),
            header=list(range(self._dataset.header_depth)),
        )
        header = parse_header(data.columns)
        data.columns = header.index
        return data
