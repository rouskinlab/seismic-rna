from __future__ import annotations

from abc import ABC, abstractmethod
from inspect import signature
from functools import cached_property, wraps
from pathlib import Path
from typing import Callable

from . import path
from .report import (OutDirF, SampleF, RefF, SectF, End5F, End3F,
                     Report, BatchReport)
from .sect import Section
from .seq import DNA


def no_kwargs(func: Callable):
    """ Prevent a function/method from accepting **kwargs. """

    @wraps(func)
    def wrapper(*args, **kwargs):
        if kwargs:
            raise TypeError(f"{func.__name__} accepts no keyword arguments, "
                            f"but got {kwargs}")
        return func(*args)

    return wrapper


def _get_kwonly(func: Callable):
    """ Get the names of the keyword-only parameters of a function. """
    return tuple(name for name, param in signature(func).parameters.items()
                 if param.kind == param.KEYWORD_ONLY)


class DataLoader(ABC):
    """ Base class for loading data from a step. """

    def __init__(self, report: Report | BatchReport):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._report: Report | BatchReport = report

    @classmethod
    @abstractmethod
    def get_report_type(cls):
        """ Type of the report for this Loader. """
        return Report

    @classmethod
    def _open_report(cls, report_file: Path):
        """ Open a report of the correct type for this class. """
        return cls.get_report_type().open(report_file)

    @property
    def out_dir(self) -> Path:
        """ Output directory. """
        return self._report.get_field(OutDirF)

    @property
    def sample(self) -> str:
        """ Name of sample. """
        return self._report.get_field(SampleF)

    @property
    def ref(self) -> str:
        """ Name of reference. """
        return self._report.get_field(RefF)

    @abstractmethod
    def get_refseq(self):
        """ Sequence of the reference. """
        return DNA(b"")

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
        return Section(ref=self.ref, refseq=self.get_refseq(),
                       end5=self._end5, end3=self._end3, name=self._sect)

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

    @abstractmethod
    def load_data_personal(self, data_file: Path, **kwargs):
        """ Load a dataset of the DataLoader itself. """

    @classmethod
    def open(cls, report_file: Path):
        """ Create a new DataLoader from a report file. """
        return cls(cls._open_report(report_file))

    def __str__(self):
        return (f"{self.__class__.__name__} for sample '{self.sample}' "
                f"over reference '{self.ref}' section '{self.sect}'")

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self._report == other._report


class BatchLoader(DataLoader, ABC):
    """ Load a dataset that is split into batches. """

    @property
    def num_batches(self):
        return self._report.num_batches

    def build_batch_path(self, batch: int):
        """ Path to the batch with the given number. """
        return self._report.build_batch_path(self.out_dir, batch,
                                             **self._report.path_fields())

    @classmethod
    @abstractmethod
    def get_report_type(cls):
        """ Type of the report for this Loader. """
        return BatchReport

    @abstractmethod
    def iter_batches_personal(self, **kwargs):
        """ Yield every batch of personal data. """
        for batch_file in self._report.iter_batch_paths():
            yield self.load_data_personal(batch_file, **kwargs)

    @abstractmethod
    def iter_batches_processed(self, **kwargs):
        """ Yield every batch of processed data. """
        # For the base BatchLoader class, no processing is performed, so
        # the "processed" data are the same as the "personal" data. The
        # distinction matters for the BatchChainLoader subclass.
        yield from self.iter_batches_personal(**kwargs)


class ChainLoader(DataLoader, ABC):
    """ Load data via a DataLoader from a previous step. """

    @classmethod
    @abstractmethod
    def get_import_type(cls):
        """ Type of the data loader that is immediately before this data
        loader in the chain and from which data are thus imported. """

    @property
    def import_path_fields(self):
        """ Fields for creating the imported data loader. """
        return {path.SAMP: self.sample, path.REF: self.ref}

    @cached_property
    def import_loader(self):
        """ Data loader that is immediately before this data loader in
        the chain and from which data are thus imported. """
        imp_type = self.get_import_type()
        if imp_type is None:
            return None
        rep_type = imp_type.get_report_type()
        return imp_type.open(rep_type.build_path(self.out_dir,
                                                 **self.import_path_fields))


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
