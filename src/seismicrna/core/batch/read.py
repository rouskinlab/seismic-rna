from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
import pandas as pd

from .index import RB_INDEX_NAMES, NO_READ, get_inverse, get_length
from ..types import fit_uint_type


class ReadBatch(ABC):
    """ Batch of reads. """

    def __init__(self, *, batch: int, **kwargs):
        super().__init__(**kwargs)
        self.batch = batch

    @cached_property
    @abstractmethod
    def read_nums(self) -> np.ndarray:
        """ Read numbers in use. """

    @cached_property
    @abstractmethod
    def max_read(self) -> int:
        """ Maximum possible value for a read index. """

    @cached_property
    @abstractmethod
    def num_reads(self) -> int:
        """ Number of reads that are actually in use. """

    @cached_property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(max(self.max_read, 0))

    @cached_property
    @abstractmethod
    def read_indexes(self) -> np.ndarray:
        """ Map each read number to its index in self.read_nums. """

    @cached_property
    def batch_read_index(self):
        """ MultiIndex of the batch number and read numbers. """
        return pd.MultiIndex.from_arrays(
            [np.broadcast_to(self.batch,
                             get_length(self.read_nums, "read_nums")),
             self.read_nums],
            names=RB_INDEX_NAMES)

    def __str__(self):
        return (f"{type(self).__name__} {self.batch} with "
                f"{get_length(self.read_nums, 'read_nums')} reads")


class AllReadBatch(ReadBatch, ABC):

    @cached_property
    def read_nums(self):
        return np.arange(self.num_reads, dtype=self.read_dtype)

    @cached_property
    def max_read(self):
        return self.num_reads - 1

    @cached_property
    def read_indexes(self):
        return self.read_nums


class PartialReadBatch(ReadBatch, ABC):

    @cached_property
    def num_reads(self):
        return get_length(self.read_nums, "read_nums")

    @cached_property
    def max_read(self):
        return self.read_nums.max(initial=NO_READ)

    @cached_property
    def read_indexes(self):
        return get_inverse(self.read_nums, "read_nums")