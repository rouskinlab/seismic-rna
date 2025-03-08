from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from inspect import Parameter, signature
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from .base import (DELET_REL,
                   INSRT_REL,
                   MATCH_REL,
                   MUTAT_REL,
                   SUBST_REL,
                   SUB_A_REL,
                   SUB_C_REL,
                   SUB_G_REL,
                   SUB_T_REL,
                   INFOR_REL,
                   TABLE_RELS,
                   Table,
                   PositionTable,
                   ReadTable,
                   AbundanceTable,
                   all_patterns)
from .. import path
from ..batch import accumulate_batches, accumulate_counts
from ..dataset import MutsDataset
from ..header import Header, make_header
from ..logs import logger
from ..seq import Region
from ..write import need_write

# These relationships are of all subtypes of mutations.
SUBMUTS = [SUBST_REL,
           SUB_A_REL,
           SUB_C_REL,
           SUB_G_REL,
           SUB_T_REL,
           DELET_REL,
           INSRT_REL]

PRECISION = 1


class Tabulator(ABC):
    """ Base class for tabulating data for one or more tables. """

    @classmethod
    @abstractmethod
    def table_types(cls) -> list[type[TableWriter]]:
        """ Types of tables that this tabulator can write. """

    @classmethod
    @cache
    def get_load_function(cls):
        """ LoadFunction for all Dataset types for this Tabulator. """
        load_functions = list(set(table_type.get_load_function()
                                  for table_type in cls.table_types()))
        assert len(load_functions) == 1
        return load_functions[0]

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    @classmethod
    def _format_table(cls, counts: pd.DataFrame, header: Header):
        """ Format the count data with the proper header. """
        # Initialize an empty table.
        table = pd.DataFrame(cls.get_null_value(), counts.index, header.index)
        # Count reads with each relationship at each position.
        for rel, rel_counts in counts.items():
            table.loc[:, rel] = rel_counts
        # Add a column for informative relationships.
        table.loc[:, INFOR_REL] = (table.loc[:, MATCH_REL].values
                                   +
                                   table.loc[:, MUTAT_REL].values)
        return table

    def __init__(self, *,
                 top: Path,
                 branches: dict[str, str],
                 sample: str,
                 region: Region,
                 count_ends: bool,
                 count_pos: bool,
                 count_read: bool,
                 validate: bool = True):
        self.top = top
        self.branches = branches
        self.sample = sample
        self.refseq = region.seq
        self.region = region
        self.pattern = None
        self.ks = None
        self.count_ends = count_ends
        self.count_pos = count_pos
        self.count_read = count_read
        self.validate = validate

    @property
    def ref(self):
        """ Name of the reference. """
        return self.region.ref

    @cached_property
    def pos_header(self):
        """ Header of the per-position data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def read_header(self):
        """ Header of the per-read data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def _accum_kwargs(self):
        """ Keyword arguments for the accumulate function. """
        return dict(refseq=self.refseq,
                    pos_nums=self.region.unmasked_int,
                    patterns=all_patterns(self.pattern),
                    ks=self.ks,
                    count_ends=self.count_ends,
                    count_pos=self.count_pos,
                    count_read=self.count_read,
                    validate=self.validate)

    @cached_property
    @abstractmethod
    def _counts(self):
        """ All counts for all table(s). """

    @property
    def num_reads(self):
        """ Raw number of reads. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return num_reads

    @property
    def counts_per_pos(self):
        """ Raw counts per position. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return per_pos

    @property
    def counts_per_read(self):
        """ Raw counts per read. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return per_read

    @property
    def end_counts(self):
        """ Raw counts for each pair of end coordinates. """
        num_reads, end_counts, per_pos, per_read = self._counts
        return end_counts

    @cached_property
    def data_per_pos(self):
        """ DataFrame of per-position data. """
        return self._format_table(self.counts_per_pos,
                                  self.pos_header).reindex(self.region.range)

    @cached_property
    def data_per_read(self):
        """ DataFrame of per-read data. """
        return self._format_table(self.counts_per_read,
                                  self.read_header)

    @cached_property
    @abstractmethod
    def data_per_clust(self) -> pd.Series | None:
        """ Series of per-cluster data (or None if no clusters). """

    def generate_tables(self, *,
                        pos: bool = True,
                        read: bool = True,
                        clust: bool = True):
        """ Generate the tables from this data. """
        for table_type in self.table_types():
            if issubclass(table_type, PositionTableWriter):
                if pos:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, ReadTableWriter):
                if read:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            elif issubclass(table_type, AbundanceTableWriter):
                if clust:
                    yield table_type(self)
                else:
                    logger.detail(f"Skipped {table_type} for {self}")
            else:
                # This should never happen; checking just in case.
                raise TypeError(table_type)

    def write_tables(self, *, force: bool = False, **kwargs):
        files = list()
        for table in self.generate_tables(**kwargs):
            files.append(table.write(force))
        return files

    @cached_property
    def _str_dict(self):
        return {path.TOP: self.top,
                path.SAMPLE: self.sample,
                path.BRANCHES: self.branches,
                path.REF: self.ref,
                path.REG: self.region.name}

    def __str__(self):
        return f"{type(self).__name__} of {self._str_dict}"


class CountTabulator(Tabulator, ABC):
    """ Tabulator that accepts pre-counted data from batches. """

    def __init__(self, *,
                 batch_counts: Iterable[tuple[Any, Any, Any, Any]],
                 **kwargs):
        super().__init__(**kwargs)
        self._batch_counts = iter(batch_counts)

    @cached_property
    def _counts(self):
        logger.routine(f"Began tabulating {self}")
        counts = accumulate_counts(self._batch_counts, **self._accum_kwargs)
        logger.routine(f"Ended tabulating {self}")
        return counts


class BatchTabulator(Tabulator, ABC):
    """ Tabulator that accepts batches as the input data. """

    def __init__(self, *,
                 get_batch_count_all: Callable,
                 num_batches: int,
                 num_cpus: int = 1,
                 **kwargs):
        super().__init__(**kwargs)
        self._get_batch_count_all = get_batch_count_all
        self.num_batches = num_batches
        self.num_cpus = num_cpus

    @cached_property
    def _counts(self):
        logger.routine(f"Began tabulating {self}")
        counts = accumulate_batches(self._get_batch_count_all,
                                    self.num_batches,
                                    num_cpus=self.num_cpus,
                                    **self._accum_kwargs)
        logger.routine(f"Ended tabulating {self}")
        return counts


class DatasetTabulator(BatchTabulator, ABC):
    """ Tabulator made from one dataset. """

    @classmethod
    def get_dataset_types(cls):
        """ Types of Dataset this Tabulator can process. """
        return cls.get_load_function().dataset_types

    @classmethod
    def _list_args(cls, func: Callable):
        """ List a function's keyword arguments with no defaults. """
        return [name for name, param in signature(func).parameters.items()
                if (param.kind is Parameter.KEYWORD_ONLY
                    and param.default is param.empty)]

    @classmethod
    @cache
    def init_kws(cls):
        """ Attributes of the dataset to use as keyword arguments in
        super().__init__(). """
        return ["top",
                "sample",
                "branches",
                "get_batch_count_all",
                "num_batches"]

    def __init__(self, *,
                 dataset: MutsDataset,
                 validate: bool = False,
                 **kwargs):
        if not isinstance(dataset, self.get_dataset_types()):
            raise TypeError(
                f"Dataset must be one of {self.get_dataset_types()}, "
                f"but got {type(dataset).__name__}"
            )
        # Since the batches come from a Dataset, they do not need to be
        # validated, so make validate False by default.
        super().__init__(**{attr: getattr(dataset, attr)
                            for attr in self.init_kws()},
                         validate=validate,
                         **kwargs)


class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: Tabulator):
        self._tabulator = tabulator

    @property
    def _attrs(self):
        return self._tabulator

    def write(self, force: bool):
        """ Write the table's rounded data to the table's CSV file. """
        if need_write(self.path, force):
            self.data.round(decimals=PRECISION).to_csv(self.path)
            logger.action(f"Wrote {self} to {self.path}")
        return self.path


class PositionTableWriter(TableWriter, PositionTable, ABC):

    @cached_property
    def data(self):
        return self._tabulator.data_per_pos


class ReadTableWriter(TableWriter, ReadTable, ABC):

    @cached_property
    def data(self):
        return self._tabulator.data_per_read


class AbundanceTableWriter(TableWriter, AbundanceTable, ABC):

    @cached_property
    def data(self):
        return self._tabulator.data_per_clust
