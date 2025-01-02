from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any, Generator, Iterable

from .base import GraphBase, GraphRunner, GraphWriter
from ..cluster.table import ClusterPositionTableLoader
from ..core.arg import opt_use_ratio, opt_quantile
from ..core.table import Table, PositionTable
from ..mask.table import MaskPositionTableLoader, MaskReadTableLoader
from ..relate.table import RelatePositionTableLoader, RelateReadTableLoader


class TableGraph(GraphBase, ABC):
    """ Graph based on one or more tables. """

    def __init__(self, *, use_ratio: bool, quantile: float, **kwargs):
        """
        Parameters
        ----------
        use_ratio: bool
            Use the ratio of the number of times the relationship occurs
            to the number of occurrances of another kind of relationship
            (which is Covered for Covered and Informative; Informative
            for all other relationships), rather than the raw count.
        quantile: float
            If `use_ratio` is True, then normalize the ratios to this
            quantile and then winsorize them to the interval [0, 1].
            Passing 0.0 disables normalization and winsorization.
        """
        super().__init__(**kwargs)
        self.use_ratio = use_ratio
        self.quantile = quantile

    @property
    def data_kind(self):
        """ Kind of data being used: either "ratio" or "count". """
        return "ratio" if self.use_ratio else "count"

    @cached_property
    def _fetch_kwargs(self) -> dict[str, Any]:
        """ Keyword arguments for self._fetch_data. """
        return dict(rel=self.rel_names)

    def _fetch_data(self, table: PositionTable, **kwargs):
        """ Fetch data from the table. """
        kwargs = self._fetch_kwargs | kwargs
        return (table.fetch_ratio(quantile=self.quantile, **kwargs)
                if self.use_ratio
                else table.fetch_count(**kwargs))

    @cached_property
    def _title_main(self):
        return [f"{self.what()} of {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over reference {repr(self.ref)} "
                f"region {repr(self.reg)}"]

    @cached_property
    def details(self):
        return ([f"quantile = {round(self.quantile, 3)}"] if self.use_ratio
                else list())

    @cached_property
    def predicate(self):
        fields = [self.codestring, self.data_kind]
        if self.use_ratio:
            fields.append(f"q{round(self.quantile * 100.)}")
        return "-".join(fields)


class TableGraphWriter(GraphWriter, ABC):

    def __init__(self, *tables: Table):
        self.tables = list(tables)

    @abstractmethod
    def iter_graphs(self,
                    *args,
                    **kwargs) -> Generator[TableGraph, None, None]:
        pass


def load_pos_tables(input_paths: Iterable[str | Path]):
    """ Load position tables. """
    paths = list(input_paths)
    for table_type in [RelatePositionTableLoader,
                       MaskPositionTableLoader,
                       ClusterPositionTableLoader]:
        yield from table_type.load_tables(paths)


def load_read_tables(input_paths: Iterable[str | Path]):
    """ Load read tables. """
    paths = list(input_paths)
    for table_type in [RelateReadTableLoader,
                       MaskReadTableLoader]:
        yield from table_type.load_tables(paths)


class TableGraphRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[TableGraphWriter]:
        pass

    @classmethod
    def var_params(cls):
        return [opt_use_ratio,
                opt_quantile]


class PosGraphRunner(TableGraphRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_pos_tables


class ReadGraphRunner(TableGraphRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_read_tables

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
