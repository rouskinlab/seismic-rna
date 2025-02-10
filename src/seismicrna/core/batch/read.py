from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
import pandas as pd

from .index import RB_INDEX_NAMES
from ..array import get_length
from ..types import fit_uint_type


class ReadBatch(ABC):
    """ Batch of reads. """

    def __init__(self, *, batch: int, **kwargs):
        super().__init__(**kwargs)
        self.batch = batch
        self.masked_read_nums = None

    @cached_property
    @abstractmethod
    def max_read(self) -> int:
        """ Maximum possible value for a read index. """

    @cached_property
    @abstractmethod
    def num_reads(self) -> int | pd.Series:
        """ Number of reads. """

    @cached_property
    @abstractmethod
    def read_nums(self) -> np.ndarray:
        """ Read numbers. """

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

    @property
    def masked_reads_bool(self):
        masked_reads_bool = np.zeros_like(self.read_nums, dtype=bool)
        if self.masked_read_nums is not None:
            masked_reads_bool[self.read_indexes[self.masked_read_nums]] = 1
        return masked_reads_bool

    def __str__(self):
        return (f"{type(self).__name__} {self.batch} with "
                f"{get_length(self.read_nums, 'read_nums')} reads")
