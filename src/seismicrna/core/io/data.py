from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from logging import getLogger
from pathlib import Path
from typing import Generator, Generic, Iterable, TypeVar

from .report import (Field,
                     SampleF,
                     RefF,
                     SectF,
                     End5F,
                     End3F,
                     NumBatchF,
                     ChecksumsF,
                     RefseqChecksumF)
from .seq import RefseqIO
from .. import path
from ..batch import list_batch_nums
from ..seq import DNA, Section

# Type variable for reports.
R = TypeVar('R')
# Type variable for data.
D = TypeVar('D')
# Type variables for datasets.
S1 = TypeVar('S1')
S2 = TypeVar('S2')

logger = getLogger(__name__)


class Dataset(Generic[D], ABC):
    """ Base class for a dataset.

    The purpose of the Dataset class is to load a dataset, which may or
    may not be split over multiple batches, each batch in one file.
    """

    @classmethod
    @abstractmethod
    def get_data_type(cls) -> type[D]:
        """ Type of the data for this dataset. """

    @property
    @abstractmethod
    def sample(self) -> str:
        """ Name of the sample. """

    @property
    @abstractmethod
    def ref(self) -> str:
        """ Name of the reference. """

    @property
    @abstractmethod
    def refseq(self) -> DNA:
        """ Sequence of the reference. """

    @property
    @abstractmethod
    def section(self) -> Section:
        """ Section of the reference. """


class UnifiedDataset(Dataset[D], ABC):
    """ Dataset in one whole piece (exactly one file). """

    @abstractmethod
    def _get_data(self) -> D:
        """ Get the data for this dataset. """

    @cached_property
    def data(self):
        data = self._get_data()
        # Verify the type of the data before returning it.
        if not isinstance(data, self.get_data_type()):
            raise TypeError(f"Expected {self.get_data_type()}, "
                            f"but got {type(data).__name__}")
        return data


class BatchedDataset(Dataset[D], ABC):
    """ Dataset split into batches (one file per batch). """

    @property
    @abstractmethod
    def num_batches(self) -> int:
        """ Number of batches. """

    @property
    def batch_nums(self):
        """ Numbers of the batches. """
        return list_batch_nums(self.num_batches)

    @abstractmethod
    def _iter_batches(self) -> Generator[D, None, None]:
        """ Yield each batch. """

    def iter_batches(self):
        """ Yield each batch. """
        for batch in self._iter_batches():
            # Verify the type of the batch before yielding it.
            if not isinstance(batch, self.get_data_type()):
                raise TypeError(f"Expected {self.get_data_type()}, "
                                f"but got {type(batch).__name__}")
            yield batch


class Loader(Generic[R], ABC):

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[R]:
        """ Type of the report for this loader. """

    def __init__(self, report: R, top: Path):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._report = report
        self.top = top


class LoadedDataset(Dataset[D], Loader[R], ABC):
    """ Dataset created by loading directly from a Report. """

    @classmethod
    def get_refseq_type(cls):
        """ Type of the reference sequence file for this loader. """
        return RefseqIO

    def _get_report_field(self, field: Field, missing_ok: bool = False):
        try:
            return self._report.get_field(field)
        except AttributeError:
            if missing_ok:
                return None
            raise

    @property
    def sample(self) -> str:
        return self._get_report_field(SampleF)

    @property
    def ref(self) -> str:
        return self._get_report_field(RefF)

    @cached_property
    def refseq(self):
        return self.get_refseq_type().load(
            self.get_refseq_type().build_path(top=self.top,
                                              sample=self.sample,
                                              ref=self.ref),
            self._get_report_field(RefseqChecksumF)
        ).refseq

    @property
    def _sect(self):
        """ Name of the section. """
        return self._get_report_field(SectF, missing_ok=True)

    @property
    def _end5(self):
        """ 5' end of the section. """
        return self._get_report_field(End5F, missing_ok=True)

    @property
    def _end3(self):
        """ 3' end of the section. """
        return self._get_report_field(End3F, missing_ok=True)

    @cached_property
    def section(self):
        """ Section of the reference. """
        return Section(ref=self.ref,
                       refseq=self.refseq,
                       end5=self._end5,
                       end3=self._end3,
                       name=self._sect)

    @classmethod
    def load(cls, report_file: Path):
        """ Create a new DatasetLoader from a report file. """
        report = cls.get_report_type().load(report_file)
        top, _ = cls.get_report_type().parse_file_path(report_file)
        return cls(report, top)

    def __str__(self):
        return (f"{type(self).__name__} for sample {repr(self.sample)} "
                f"over {self.section}")


class UnifiedLoadedDataset(LoadedDataset[D, R], UnifiedDataset[D], ABC):
    """ Dataset loaded as one piece from a report. """


class BatchedLoadedDataset(LoadedDataset[D, R], BatchedDataset[D], ABC):
    """ Dataset loaded in batches from a report. """

    @classmethod
    def get_btype_name(cls):
        """ Name of the type of batch for this Loader. """
        return cls.get_data_type().btype()

    @cached_property
    def num_batches(self):
        return self._get_report_field(NumBatchF)

    @cache
    def get_batch_path(self, batch: int):
        """ Get the path to a batch of a specific number. """
        fields = self._report.path_fields(self.top,
                                          self.get_data_type().auto_fields())
        return self.get_data_type().build_path(batch=batch, **fields)

    @cache
    def report_checksum(self, batch: int):
        """ Get the checksum of a specific batch from the report. """
        return self._report.get_field(ChecksumsF)[self.get_btype_name()][batch]

    def load_batch(self, batch: int):
        """ Load a specific batch of data. """
        return self.get_data_type().load(self.get_batch_path(batch),
                                         self.report_checksum(batch))

    def _iter_batches(self):
        for batch in self.batch_nums:
            yield self.load_batch(batch)


class Linker(Generic[S1, S2], ABC):

    @classmethod
    @abstractmethod
    def get_data1_type(cls) -> type[S1]:
        """ Type of Dataset 1. """

    @classmethod
    @abstractmethod
    def get_data2_type(cls) -> type[S2]:
        """ Type of Dataset 2. """

    @classmethod
    def verify_data_types(cls, data1: S1, data2: S2):
        if not isinstance(data1, cls.get_data1_type()):
            raise TypeError(f"Expected a {cls.get_data1_type().__name__} for "
                            f"data1, but got {type(data1).__name__}")
        if not isinstance(data2, cls.get_data2_type()):
            raise TypeError(f"Expected a {cls.get_data2_type().__name__} for "
                            f"data2, but got {type(data2).__name__}")

    def __init__(self, data1: S1, data2: S2):
        self.verify_data_types(data1, data2)
        self.data1 = data1
        self.data2 = data2


class LinkedDataset(Dataset[D], Linker[S1, S2], ABC):
    """ A Dataset created with a function that accepts two datasets and
    returns a third "linked" dataset. """

    @classmethod
    def _unify_attrs(cls, attr1, attr2):
        if attr1 != attr2:
            raise ValueError(
                f"Inconsistent attribute: {repr(attr1)} ≠ {repr(attr2)}")
        return attr1

    @property
    def sample(self):
        return self._unify_attrs(self.data1.sample, self.data2.sample)

    @property
    def ref(self):
        return self._unify_attrs(self.data1.ref, self.data2.ref)

    @property
    def refseq(self):
        return self._unify_attrs(self.data1.refseq, self.data2.refseq)

    @property
    def section(self):
        return self._unify_attrs(self.data1.section, self.data2.section)

    @abstractmethod
    def link(self):
        """ Link the datasets. """


class UnifiedLinkedDataset(LinkedDataset[D, S1, S2], UnifiedDataset[D], ABC):
    """ Linked unified dataset. """

    @classmethod
    @abstractmethod
    def _link(cls, data1, data2) -> D:
        """ Link the data in the two datasets. """

    def _get_data(self):
        return self._link(self.data1, self.data2)


class BatchedLinkedDataset(LinkedDataset[D, S1, S2], BatchedDataset[D], ABC):
    """ Linked batched dataset. """

    @classmethod
    @abstractmethod
    def _link(cls, batch1, batch2) -> D:
        """ Link corresponding batches of data. """

    @property
    def num_batches(self):
        return self._unify_attrs(self.data1.num_batches, self.data2.num_batches)

    def _iter_batches(self):
        for batch1, batch2 in zip(self.data1.iter_batches(),
                                  self.data2.iter_batches(),
                                  strict=True):
            yield self._link(batch1, batch2)


def load_data(report_files: Iterable[Path], loader_type: type[LoadedDataset]):
    """ Load the data for each report file. """
    for report_file in path.deduplicated(report_files):
        try:
            yield loader_type.load(report_file)
        except Exception as error:
            logger.error(f"Failed to open {report_file}: {error}")

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
