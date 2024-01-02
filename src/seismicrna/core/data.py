from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Generic, Iterable, TypeVar

from . import path
from .batch import list_batch_nums
from .io import RefseqIO, convert_path
from .rel import RelPattern
from .report import (SampleF,
                     RefF,
                     SectF,
                     End5F,
                     End3F,
                     NumBatchF,
                     ChecksumsF,
                     RefseqChecksumF)
from .seq import FULL_NAME, DNA, Section, hyphenate_ends

# Type variable for reports.
R = TypeVar("R")
# Type variable for data.
D = TypeVar("D")
# Type variables for datasets.
S1 = TypeVar("S1")
S2 = TypeVar("S2")

logger = getLogger(__name__)


class Dataset(Generic[D], ABC):
    """ Handle to a dataset comprising batches of data. """

    @classmethod
    @abstractmethod
    def get_batch_type(cls) -> type[D]:
        """ Type of the batch for this dataset. """

    @property
    @abstractmethod
    def top(self) -> Path:
        """ Top-level directory of the dataset. """

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
    def end5(self) -> int:
        """ 5' end of the section. """

    @property
    @abstractmethod
    def end3(self) -> int:
        """ 3' end of the section. """

    @property
    @abstractmethod
    def sect(self) -> str:
        """ Name of the section. """

    @property
    @abstractmethod
    def pattern(self) -> RelPattern | None:
        """ Pattern of mutations to count. """

    @property
    @abstractmethod
    def num_batches(self) -> int:
        """ Number of batches. """

    @property
    def batch_nums(self):
        """ Numbers of the batches. """
        return list_batch_nums(self.num_batches)

    @abstractmethod
    def get_batch(self, batch_num: int) -> D:
        """ Get a specific batch of data. """

    def iter_batches(self):
        """ Yield each batch. """
        for batch_num in self.batch_nums:
            yield self.get_batch(batch_num)

    def __str__(self):
        return f"{type(self).__name__} for sample {repr(self.sample)}"


class MutsDataset(Dataset[D], ABC):
    """ Dataset with explicit mutational data. """

    @property
    @abstractmethod
    def refseq(self) -> DNA:
        """ Sequence of the reference. """

    @property
    def reflen(self):
        """ Length of the reference sequence. """
        return len(self.refseq)

    @cached_property
    def section(self):
        """ Section of the dataset. """
        return Section(ref=self.ref,
                       seq=self.refseq,
                       end5=self.end5,
                       end3=self.end3,
                       name=self.sect)


class Loader(Generic[R], ABC):

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[R]:
        """ Type of the report for this loader. """


class LoadedDataset(Dataset[D], Loader[R], ABC):
    """ Dataset created by loading directly from a Report. """

    @classmethod
    def get_btype_name(cls):
        """ Name of the type of batch for this Loader. """
        return cls.get_batch_type().btype()

    @classmethod
    def load(cls, report_file: Path):
        """ Create a new DatasetLoader from a report file. """
        report = cls.get_report_type().load(report_file)
        top, _ = cls.get_report_type().parse_path(report_file)
        return cls(report, top)

    def __init__(self, report: R, top: Path):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self.report = report
        self._top = top

    @property
    def top(self):
        return self._top

    @property
    def sample(self):
        return self.report.get_field(SampleF)

    @property
    def ref(self):
        return self.report.get_field(RefF)

    @property
    def end5(self):
        return self.report.get_field(End5F)

    @property
    def end3(self):
        return self.report.get_field(End3F)

    @property
    def sect(self):
        return self.report.get_field(SectF)

    @cached_property
    def num_batches(self):
        return self.report.get_field(NumBatchF)

    def get_batch_path(self, batch: int):
        """ Get the path to a batch of a specific number. """
        fields = self.report.path_fields(self.top,
                                         self.get_batch_type().auto_fields())
        return self.get_batch_type().build_path(batch=batch, **fields)

    def report_checksum(self, batch: int):
        """ Get the checksum of a specific batch from the report. """
        return self.report.get_field(ChecksumsF)[self.get_btype_name()][batch]

    def get_batch(self, batch_num: int):
        batch = self.get_batch_type().load(self.get_batch_path(batch_num),
                                           self.report_checksum(batch_num))
        if batch.batch != batch_num:
            raise ValueError(f"Expected batch to have number {batch_num}, "
                             f"but got {batch.batch}")
        return batch


class LoadedMutsDataset(LoadedDataset[D, R], MutsDataset, ABC):

    @classmethod
    def get_refseq_type(cls):
        """ Type of the reference sequence file for this loader. """
        return RefseqIO

    @cached_property
    def refseq(self):
        return self.get_refseq_type().load(
            self.get_refseq_type().build_path(top=self.top,
                                              sample=self.sample,
                                              ref=self.ref),
            self.report.get_field(RefseqChecksumF)
        ).refseq

    @property
    def end5(self):
        try:
            return super().end5
        except AttributeError:
            return 1

    @property
    def end3(self):
        try:
            return super().end3
        except AttributeError:
            return self.reflen

    @property
    def sect(self):
        try:
            return super().sect
        except AttributeError:
            return (FULL_NAME if self.end5 == 1 and self.end3 == self.reflen
                    else hyphenate_ends(self.end5, self.end3))


class Linker(Generic[S1, S2], ABC):

    @classmethod
    @abstractmethod
    def get_dataset1_type(cls) -> type[S1]:
        """ Type of Dataset 1. """

    @classmethod
    @abstractmethod
    def get_dataset2_type(cls) -> type[S2]:
        """ Type of Dataset 2. """

    @classmethod
    def verify_data_types(cls, data1: S1, data2: S2):
        if not isinstance(data1, cls.get_dataset1_type()):
            raise TypeError(f"Expected a {cls.get_dataset1_type().__name__} for "
                            f"data1, but got {type(data1).__name__}")
        if not isinstance(data2, cls.get_dataset2_type()):
            raise TypeError(f"Expected a {cls.get_dataset2_type().__name__} for "
                            f"data2, but got {type(data2).__name__}")


class LinkedDataset(Dataset[D], Linker[S1, S2], ABC):
    """ A Dataset created with a function that accepts two datasets and
    returns a third "linked" dataset. """

    @classmethod
    def get_report_type(cls):
        return cls.get_dataset2_type().get_report_type()

    @classmethod
    def load(cls, report_file: Path):
        """ Create a new LinkedDataset from a report file. """
        data1_type = cls.get_dataset1_type()
        data1 = data1_type.load(convert_path(report_file,
                                             cls.get_report_type(),
                                             data1_type.get_report_type()))
        data2 = cls.get_dataset2_type().load(report_file)
        return cls(data1, data2)

    def __init__(self, data1: S1, data2: S2):
        self.verify_data_types(data1, data2)
        self.data1 = data1
        self.data2 = data2

    def _get_data_attr(self, name: str):
        val1 = getattr(self.data1, name)
        val2 = getattr(self.data2, name)
        if val1 != val2:
            raise ValueError(f"Inconsistent values for {repr(name)} in "
                             f"{self.data1} ({val1}) and {self.data2} ({val2})")
        return val1

    @property
    def top(self):
        return self._get_data_attr("top")

    @property
    def sample(self):
        return self._get_data_attr("sample")

    @property
    def ref(self):
        return self._get_data_attr("ref")

    @property
    def refseq(self):
        return self.data1.refseq

    @abstractmethod
    def _link(self, batch1, batch2) -> D:
        """ Link corresponding batches of data. """

    @property
    def num_batches(self):
        return self._get_data_attr("num_batches")

    def get_batch(self, batch_num: int):
        return self._link(self.data1.get_batch(batch_num),
                          self.data2.get_batch(batch_num))


class LinkedMutsDataset(LinkedDataset[D, S1, S2], MutsDataset, ABC):

    @property
    def end5(self):
        return self.data2.end5

    @property
    def end3(self):
        return self.data2.end3

    @property
    def sect(self):
        return self.data2.sect


class PooledDataset(Dataset[D], ABC):
    """ Dataset made of one or more other datasets (all the same type)
    with the same top-level directory, reference, and section (but not
    necessarily the same sample) pooled together. """

    @classmethod
    @abstractmethod
    def get_dataset_type(cls) -> type[Dataset[D]]:
        """ Type of the datasets that are pooled together. """

    @classmethod
    def get_batch_type(cls):
        return cls.get_dataset_type().get_batch_type()

    def __init__(self, sample: str, datasets: Iterable[Dataset[D]]):
        self._sample = sample
        self._datasets = list()
        for dataset in datasets:
            if not isinstance(dataset, self.get_dataset_type()):
                raise TypeError(f"{type(self).__name__} expected each dataset "
                                f"to be {self.get_dataset_type().__name__}, "
                                f"but got {type(dataset).__name__}")
            self._datasets.append(dataset)
        if not self._datasets:
            raise ValueError(f"{type(self).__name__} got no datasets")

    def _list_dataset_attr(self, name: str):
        """ Get a list of an attribute for each dataset. """
        return [getattr(dataset, name) for dataset in self._datasets]

    def _get_common_attr(self, name: str):
        """ Get a common attribute among datasets. """
        attrs = list(set(self._list_dataset_attr(name)))
        if len(attrs) != 1:
            raise ValueError(f"{type(self).__name__} got multiple values for "
                             f"attribute {repr(name)}: {attrs}")
        return attrs[0]

    @property
    def top(self):
        return self._get_common_attr("top")

    @property
    def sample(self):
        return self._sample

    @cached_property
    def samples(self) -> list[str]:
        """ Name of the sample for each dataset in the pool. """
        return self._list_dataset_attr("sample")

    @property
    def ref(self):
        return self._get_common_attr("ref")

    @property
    def end5(self):
        return self._get_common_attr("end5")

    @property
    def end3(self):
        return self._get_common_attr("end3")

    @property
    def sect(self):
        return self._get_common_attr("sect")

    @property
    def pattern(self):
        return self._get_common_attr("pattern")

    @cached_property
    def nums_batches(self) -> list[int]:
        """ Number of batches in each dataset in the pool. """
        return self._list_dataset_attr("num_batches")

    @cached_property
    def num_batches(self):
        return sum(self.nums_batches)

    def _translate_batch_num(self, batch_num: int):
        """ Translate a batch number into the numbers of the dataset and
        of the batch within the dataset. """
        dataset_batch_num = batch_num
        for dataset_num, num_batches in enumerate(self.nums_batches):
            if dataset_batch_num < 0:
                raise ValueError(f"Invalid batch number: {dataset_batch_num}")
            if dataset_batch_num < num_batches:
                return dataset_num, dataset_batch_num
            dataset_batch_num -= num_batches
        raise ValueError(f"Batch {batch_num} is invalid for {self} with "
                         f"{self.num_batches} batches")

    def get_batch(self, batch_num: int):
        dataset_num, dataset_batch_num = self._translate_batch_num(batch_num)
        return self._datasets[dataset_num].get_batch(dataset_batch_num)


def load_data(report_files: Iterable[Path], loader_type: type[LoadedDataset]):
    """ Load the data for each report file. """
    for report_file in path.deduplicated(report_files):
        try:
            yield loader_type.load(report_file)
        except Exception as error:
            logger.error(f"Failed to open {report_file}: {error}")

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
