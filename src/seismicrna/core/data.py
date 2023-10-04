from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Generator, TypeVar

from .batch import Batch, list_batch_nums
from .refseq import RefseqFile
from .report import (SampleF,
                     RefF,
                     SectF,
                     End5F,
                     End3F,
                     NumBatchF,
                     ChecksumsF,
                     RefseqChecksumF,
                     RefseqReport,
                     BatchReport)
from .sect import Section
from .seq import DNA


AnyData = TypeVar("AnyData")


class Data(ABC):
    """ Base class for a dataset. """

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


class Batches(Data, ABC):
    """ Dataset split into batches. """

    @property
    @abstractmethod
    def num_batches(self) -> int:
        """ Number of batches. """

    @property
    def batch_nums(self):
        """ Numbers of the batches. """
        return list_batch_nums(self.num_batches)

    @abstractmethod
    def iter_batches(self) -> Generator[Batch, None, None]:
        """ Yield each batch. """


class Loader(Data, ABC):
    """ Data loaded from a report. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[RefseqReport]:
        """ Type of the report for this loader. """

    def __init__(self, report: RefseqReport, top: Path):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._report: RefseqReport | BatchReport = report
        self._top = top

    @property
    def sample(self) -> str:
        return self._report.get_field(SampleF)

    @property
    def ref(self) -> str:
        return self._report.get_field(RefF)

    @cached_property
    def refseq(self):
        return RefseqFile.load(RefseqFile.build_path(top=self._top,
                                                     sample=self.sample,
                                                     ref=self.ref),
                               self._report.get_field(RefseqChecksumF)).refseq

    @property
    def _sect(self):
        """ Name of the section. """
        try:
            return self._report.get_field(SectF)
        except AttributeError:
            return None

    @property
    def _end5(self):
        """ 5' end of the section. """
        try:
            return self._report.get_field(End5F)
        except AttributeError:
            return None

    @property
    def _end3(self):
        """ 3' end of the section. """
        try:
            return self._report.get_field(End3F)
        except AttributeError:
            return None

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
        """ Create a new Loader from a report file. """
        report = cls.get_report_type().load(report_file)
        top, _ = cls.get_report_type().parse_file_path(report_file)
        return cls(report, top)

    def __str__(self):
        return (f"{type(self).__name__} for sample {repr(self.sample)} "
                f"over {self.section}")


class BatchLoader(Loader, Batches, ABC):
    """ Data loaded in batches from a report. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[BatchReport]:
        """ Type of the report for this Loader. """

    @classmethod
    def get_batch_type(cls) -> type[Batch]:
        """ Type of the batch for this Loader. """

    @classmethod
    def get_btype_name(cls):
        """ Name of the type of batch for this Loader. """
        return cls.get_batch_type().btype()

    @cached_property
    def num_batches(self):
        return self._report.get_field(NumBatchF)

    @cache
    def get_batch_path(self, batch: int):
        """ Get the path to a batch of a specific number. """
        fields = self._report.path_fields(self._top,
                                          self.get_batch_type().auto_fields())
        return self.get_batch_type().build_path(batch=batch, **fields)

    @cache
    def report_checksum(self, batch: int):
        """ Get the checksum of a specific batch from the report. """
        return self._report.get_field(ChecksumsF)[self.get_btype_name()][batch]

    def load_batch(self, batch: int):
        """ Load a specific batch of data. """
        return self.get_batch_type().load(self.get_batch_path(batch),
                                          self.report_checksum(batch))

    def iter_batches(self):
        for batch in self.batch_nums:
            yield self.load_batch(batch)


class Linker(Data, ABC):
    """ Data created by linking two other datasets. """

    @classmethod
    @abstractmethod
    def data1_type(cls) -> type[AnyData]:
        """ Type of dataset 1. """

    @classmethod
    @abstractmethod
    def data2_type(cls) -> type[AnyData]:
        """ Type of dataset 2. """

    @classmethod
    def verify_data_types(cls, data1: AnyData, data2: AnyData):
        if not isinstance(data1, cls.data1_type()):
            raise TypeError(f"Expected a {cls.data1_type().__name__} for "
                            f"data1, but got {type(data1).__name__}")
        if not isinstance(data2, cls.data2_type()):
            raise TypeError(f"Expected a {cls.data2_type().__name__} for "
                            f"data2, but got {type(data2).__name__}")

    def __init__(self, data1: AnyData, data2: AnyData):
        self.verify_data_types(data1, data2)
        self._data1 = data1
        self._data2 = data2

    @abstractmethod
    def _link(self, *args, **kwargs):
        """ Link the datasets. """


class BatchLinker(Linker, Batches, ABC):

    @classmethod
    @abstractmethod
    def data1_type(cls) -> type[Batches]:
        pass

    @classmethod
    @abstractmethod
    def data2_type(cls) -> type[Batches]:
        pass

    @abstractmethod
    def _link(self, batch1: Batch, batch2: Batch, *args, **kwargs):
        """ Link corresponding batches of datasets. """

    def iter_batches(self, *args, **kwargs):
        for batch1, batch2 in zip(self._data1.iter_batches(),
                                  self._data2.iter_batches(),
                                  strict=True):
            yield self._link(batch1, batch2, *args, **kwargs)


'''
class ChainLoader(DataLoader, ABC):
    """ Load data via a DataLoader from a previous step. """

    @classmethod
    @abstractmethod
    def get_import_type(cls) -> type[ChainLoader]:
        """ Type of the data loader that is immediately before this data
        loader in the chain and from which data are thus imported. """

    @property
    def import_path_fields(self):
        """ Fields for creating the imported data loader. """
        return {path.TOP: self._top, path.SAMP: self.sample, path.REF: self.ref}

    @cached_property
    def import_loader(self):
        """ Data loader that is immediately before this data loader in
        the chain and from which data are thus imported. """
        itype = self.get_import_type()
        if itype is None:
            return None
        rtype = itype.get_report_type()
        return itype.load(rtype.build_path(**self.import_path_fields))


class BatchChainLoader(ChainLoader, BatchLoader, ABC):
    """ Load data via a BatchLoader from a previous step. """

    @abstractmethod
    def process_batch(self, imported_batch, personal_batch, **kwargs):
        """ Return a batch of processed data from one of data imported
        from another DataLoader and one batch of personal data. """

    @abstractmethod
    def iter_batches_processed(self, **kwargs):
        # Keyword arguments of self.import_loader.iter_batches_processed
        imp_kwonly = _get_kwonly(self.import_loader.iter_batches_processed)
        imp_kwargs = {name: kwargs.pop(name) for name in list(kwargs.keys())
                      if name in imp_kwonly}
        imp_batches = self.import_loader.iter_batches_processed(**imp_kwargs)
        # Keyword arguments of self.iter_batches_personal
        pers_kwonly = _get_kwonly(self.iter_batches_personal)
        pers_kwargs = {name: kwargs.pop(name) for name in list(kwargs.keys())
                       if name in pers_kwonly}
        pers_batches = self.iter_batches_personal(**pers_kwargs)
        # Keyword arguments of self.process_batch
        proc_kwonly = _get_kwonly(self.process_batch)
        proc_kwargs = {name: kwargs.pop(name) for name in list(kwargs.keys())
                       if name in proc_kwonly}
        # Check for extraneous keyword arguments.
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {kwargs}")
        # Yield every batch of processed data.
        for imported, personal in zip(imp_batches, pers_batches, strict=True):
            yield self.process_batch(imported, personal, **proc_kwargs)
'''

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
