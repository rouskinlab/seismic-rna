from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any, Generator, Iterable

from .base import BaseGraph, BaseRunner, BaseWriter
from .rel import RelGraph, RelRunner
from ..cluster.data import (ClusterPositionTableLoader,
                            ClusterAbundanceTableLoader)
from ..core.arg import opt_use_ratio, opt_quantile, opt_verify_times
from ..core.table import Table, PositionTable
from ..mask.table import MaskPositionTableLoader, MaskReadTableLoader
from ..relate.table import RelatePositionTableLoader, RelateReadTableLoader


class TableGraph(BaseGraph, ABC):
    """ Graph based on one or more tables. """

    def __init__(self, *, use_ratio: bool, **kwargs):
        """
        Parameters
        ----------
        use_ratio: bool
            Use the ratio of the number of times the relationship occurs
            to the number of occurrances of another kind of relationship
            (which is Covered for Covered and Informative; Informative
            for all other relationships), rather than the raw count.
        """
        super().__init__(**kwargs)
        self.use_ratio = use_ratio

    @property
    def data_kind(self):
        """ Kind of data being used: either "ratio" or "count". """
        return "ratio" if self.use_ratio else "count"


class RelTableGraph(TableGraph, RelGraph, ABC):

    def __init__(self, *, quantile: float, **kwargs):
        """
        Parameters
        ----------
        quantile: float
            If `use_ratio` is True, then normalize the ratios to this
            quantile and then winsorize them to the interval [0, 1].
            Passing 0.0 disables normalization and winsorization.
        """
        super().__init__(**kwargs)
        self.quantile = quantile

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
        return super().details + ([f"quantile = {round(self.quantile, 3)}"]
                                  if self.use_ratio
                                  else list())

    @cached_property
    def predicate(self):
        fields = [self.codestring, self.data_kind]
        if self.use_ratio:
            fields.append(f"q{round(self.quantile * 100.)}")
        return super().predicate + ["-".join(fields)]


class TableWriter(BaseWriter, ABC):

    def __init__(self, *tables: Table, **kwargs):
        super().__init__(**kwargs)
        self.tables = list(tables)

    @abstractmethod
    def iter_graphs(self,
                    *args,
                    **kwargs) -> Generator[TableGraph, None, None]:
        pass


def load_pos_tables(input_paths: Iterable[str | Path], **kwargs):
    """ Load position tables. """
    paths = list(input_paths)
    for table_type in [RelatePositionTableLoader,
                       MaskPositionTableLoader,
                       ClusterPositionTableLoader]:
        yield from table_type.load_tables(paths, **kwargs)


def load_read_tables(input_paths: Iterable[str | Path], **kwargs):
    """ Load read tables. """
    paths = list(input_paths)
    for table_type in [RelateReadTableLoader,
                       MaskReadTableLoader]:
        yield from table_type.load_tables(paths, **kwargs)


def load_abundance_tables(input_paths: Iterable[str | Path], **kwargs):
    """ Load read tables. """
    paths = list(input_paths)
    for table_type in [ClusterAbundanceTableLoader]:
        yield from table_type.load_tables(paths, **kwargs)


class TableRunner(BaseRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[TableWriter]:
        pass

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_verify_times, opt_use_ratio]


class RelTableRunner(RelRunner, TableRunner, ABC):

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_quantile]


class PositionTableRunner(RelTableRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_pos_tables


class ReadTableRunner(RelTableRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_read_tables


class AbundanceTableRunner(TableRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_abundance_tables
