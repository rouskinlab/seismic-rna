from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
import pandas as pd

from .util import INDEX_NAMES, get_length
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

    @property
    @abstractmethod
    def num_reads(self) -> int:
        """ Number of reads that are actually in use. """

    @property
    @abstractmethod
    def max_read(self) -> int:
        """ Maximum possible value for a read index. """

    @property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(self.max_read)

    @cached_property
    @abstractmethod
    def read_idx(self) -> np.ndarray:
        """ Map each read number to its index in self.read_nums. """

    @cached_property
    def multiindex(self):
        """ MultiIndex of the batch number and read numbers. """
        return pd.MultiIndex.from_arrays([np.broadcast_to(self.batch,
                                                          self.num_reads),
                                          self.read_nums],
                                         names=INDEX_NAMES)

    def __str__(self):
        return f"{type(self).__name__} {self.batch} with {self.num_reads} reads"


class AllReadsBatch(ReadBatch, ABC):

    @property
    def max_read(self):
        return self.num_reads - 1

    @cached_property
    def read_nums(self):
        return np.arange(self.num_reads)

    @cached_property
    def read_idx(self):
        return self.read_nums


class MaskReadBatch(ReadBatch):

    def __init__(self, *, read_nums: np.ndarray, max_read: int, **kwargs):
        self._max_read = max_read
        super().__init__(**kwargs)
        if read_nums.size > 0 and np.max(read_nums) > max_read:
            raise ValueError(f"All read numbers must be ≤ {max_read}, but got "
                             f"{read_nums[read_nums > max_read]}")
        self._read_nums = np.asarray(read_nums, fit_uint_type(max_read))

    @property
    def read_nums(self):
        return self._read_nums

    @property
    def num_reads(self):
        return get_length(self.read_nums, "read_nums")

    @property
    def max_read(self):
        return self._max_read

    @cached_property
    def read_idx(self):
        read_map = np.zeros(self.max_read + 1, dtype=self.read_dtype)
        read_map[self.read_nums] = np.arange(self.num_reads)
        return read_map
