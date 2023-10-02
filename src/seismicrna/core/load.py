from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path

from . import path
from .batch import ReadBatch
from .files import load_pkl_br
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


class DataLoader(ABC):
    """ Base class for loading data from a report. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[RefseqReport]:
        """ Type of the report for this Loader. """

    def __init__(self, report: RefseqReport, top: Path):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._report: RefseqReport | BatchReport = report
        self._top = top

    @property
    def sample(self) -> str:
        """ Name of sample. """
        return self._report.get_field(SampleF)

    @property
    def ref(self) -> str:
        """ Name of reference. """
        return self._report.get_field(RefF)

    @cached_property
    def refseq_path(self):
        return path.build(*self._report.refseq_seg_types(),
                          **self._report.refseq_auto_fields(),
                          top=self._top,
                          sample=self.sample,
                          ref=self.ref)

    @cached_property
    def refseq(self):
        """ Sequence of the reference. """
        return RefseqFile.load(self.refseq_path,
                               self._report.get_field(RefseqChecksumF)).refseq

    @cached_property
    def _sect(self):
        """ Name of the section. """
        try:
            return self._report.get_field(SectF)
        except AttributeError:
            return None

    @cached_property
    def _end5(self):
        """ 5' end of the section. """
        try:
            return self._report.get_field(End5F)
        except AttributeError:
            return None

    @cached_property
    def _end3(self):
        """ 3' end of the section. """
        try:
            return self._report.get_field(End3F)
        except AttributeError:
            return None

    @property
    def section(self):
        """ Section of the reference. """
        return Section(ref=self.ref,
                       refseq=self.refseq,
                       end5=self._end5,
                       end3=self._end3,
                       name=self._sect)

    @property
    def sect(self) -> str:
        """ Name of the section. """
        return self.section.name if self._sect is None else self._sect

    @property
    def end5(self) -> int:
        """ 5' end of the section. """
        return self.section.end5 if self._end5 is None else self._end5

    @property
    def end3(self) -> int:
        """ 3' end of the section. """
        return self.section.end3 if self._end3 is None else self._end3

    @property
    def seq(self):
        """ Sequence of the section. """
        return self.section.seq

    @classmethod
    def load(cls, report_file: Path):
        """ Create a new DataLoader from a report file. """
        report = cls.get_report_type().load(report_file)
        top, _ = cls.get_report_type().parse_file_path(report_file)
        return cls(report, top)

    def __str__(self):
        return (f"{type(self).__name__} for sample {repr(self.sample)} "
                f"over reference {repr(self.ref)} section {repr(self.sect)}")

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self._report == other._report


class BatchLoader(DataLoader, ABC):
    """ Load a dataset that is split into batches. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[BatchReport]:
        """ Type of the report for this Loader. """

    @classmethod
    def get_batch_type(cls) -> type[ReadBatch]:
        """ Type of the batch for this Loader. """

    @classmethod
    def get_btype_name(cls):
        """ Name of the type of batch for this Loader. """
        return cls.get_batch_type().btype()

    @property
    def batch_nums(self):
        """ Batch numbers. """
        return range(self._report.get_field(NumBatchF))

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
        return load_pkl_br(self.get_batch_path(batch),
                           check_type=self.get_batch_type(),
                           checksum=self.report_checksum(batch))

    def iter_batches(self):
        """ Yield every batch of personal data. """
        for batch in self.batch_nums:
            yield self.load_batch(batch)


class Chain(ABC):

    @classmethod
    @abstractmethod
    def loader1_type(cls) -> type[BatchLoader]:
        """ Type of the first loader. """

    @classmethod
    @abstractmethod
    def loader2_type(cls) -> type[BatchLoader]:
        """ Type of the second loader. """

    def __init__(self, loader1: BatchLoader, loader2: BatchLoader):
        if not isinstance(loader1, self.loader1_type()):
            raise TypeError(f"Expected a {self.loader1_type().__name__} for "
                            f"loader1, but got {type(loader1).__name__}")
        self._loader1 = loader1
        if not isinstance(loader2, self.loader2_type()):
            raise TypeError(f"Expected a {self.loader2_type().__name__} for "
                            f"loader2, but got {type(loader2).__name__}")
        self._loader2 = loader2

    @abstractmethod
    def _process(self, batch1: ReadBatch, batch2: ReadBatch, *args, **kwargs):
        """ Process corresponding batches from both batch loaders. """

    def process_batches(self, *args, **kwargs):
        """ For each pair of batches from loader1 and loader2, yield a
        processed batch. """
        for batch1, batch2 in zip(self._loader1.iter_batches(),
                                  self._loader2.iter_batches(),
                                  strict=True):
            yield self._process(batch1, batch2, *args, **kwargs)


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
