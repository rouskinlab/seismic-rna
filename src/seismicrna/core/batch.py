import re
from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .output import PickleOutput
from .rel import REL_TYPE
from .sect import seq_pos_to_index
from .seq import DNA
from .types import fit_uint_type

logger = getLogger(__name__)

BATCH_INDEX = 0
POS_INDEX = 1
READ_INDEX = 0

READ_NUM = "Read Number"


def get_length(array: np.ndarray, what: str) -> int:
    if array.ndim != 1:
        raise ValueError(f"{what} must have 1 dimension, but got {array.ndim}")
    length, = array.shape
    return length


def ensure_same_length(arr1: np.ndarray,
                       arr2: np.ndarray,
                       what1: str,
                       what2: str):
    if (len1 := get_length(arr1, what1)) != (len2 := get_length(arr2, what2)):
        raise ValueError(
            f"Lengths differ between {what1} ({len1}) and {what2} ({len2})")
    return len1


def ensure_order(read_nums: np.ndarray,
                 array1: np.ndarray,
                 array2: np.ndarray,
                 what1: str,
                 what2: str,
                 gt_eq: bool = False):
    ensure_same_length(read_nums, array1, "read numbers", what1)
    ensure_same_length(array1, array2, what1, what2)
    ineq_func, ineq_sign = (np.less, '<') if gt_eq else (np.greater, '>')
    if np.any(is_err := ineq_func(array1, array2)):
        index = pd.Index(read_nums[is_err], name=READ_NUM)
        errors = pd.DataFrame.from_dict({what1: pd.Series(array1[is_err],
                                                          index=index),
                                         what2: pd.Series(array2[is_err],
                                                          index=index)})
        raise ValueError(f"Got {what1} {ineq_sign} {what2}:\n{errors}")


def sanitize_ends(read_nums: np.ndarray,
                  end5s: list[int] | np.ndarray,
                  mid5s: list[int] | np.ndarray,
                  mid3s: list[int] | np.ndarray,
                  end3s: list[int] | np.ndarray,
                  max_pos: int):
    pos_dtype = fit_uint_type(max_pos)
    end5s = np.asarray(end5s, pos_dtype)
    mid5s = np.asarray(mid5s, pos_dtype)
    mid3s = np.asarray(mid3s, pos_dtype)
    end3s = np.asarray(end3s, pos_dtype)
    # Verify 5' end ≥ min position
    ensure_order(read_nums,
                 end5s,
                 np.broadcast_to(POS_INDEX, end5s.shape),
                 "5' end positions",
                 f"minimum position ({POS_INDEX})",
                 gt_eq=True)
    # Verify 5' end ≤ 5' mid
    ensure_order(read_nums,
                 end5s,
                 mid5s,
                 "5' end positions",
                 "5' middle positions")
    # Verify 5' end ≤ 3' mid
    ensure_order(read_nums,
                 end5s,
                 mid3s,
                 "5' end positions",
                 "3' middle positions")
    # Verify 5' mid ≤ 3' end
    ensure_order(read_nums,
                 mid5s,
                 end3s,
                 "5' middle positions",
                 "3' end positions")
    # Verify 3' mid ≤ 3' end
    ensure_order(read_nums,
                 mid3s,
                 end3s,
                 "3' middle positions",
                 "3' end positions")
    # Verify 3' end ≤ max position
    ensure_order(read_nums,
                 end3s,
                 np.broadcast_to(max_pos, end3s.shape),
                 "3' end positions",
                 f"maximum position ({max_pos})")
    return end5s, mid5s, mid3s, end3s


class ReadBatch(PickleOutput, ABC):
    """ """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z]*)batch$", cls.__name__.lower()).groups()
        return btype

    def __init__(self, batch: int, sample: str, ref: str):
        self.batch = batch
        self.sample = sample
        self.ref = ref

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
    def read_nums(self) -> np.ndarray:
        """ Read numbers in use. """


class ArrayBatch(ReadBatch, ABC):

    def __init__(self,
                 *args,
                 max_pos: int,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.max_pos = max_pos
        (self.end5s,
         self.mid5s,
         self.mid3s,
         self.end3s) = sanitize_ends(self.read_nums,
                                     end5s,
                                     mid5s,
                                     mid3s,
                                     end3s,
                                     self.max_pos)

    @cached_property
    @abstractmethod
    def pos(self) -> np.ndarray:
        """ Array of 1-indexed positions. """

    @property
    @abstractmethod
    def num_pos(self) -> int:
        """ Number of positions that are actually in use. """

    @property
    def pos_dtype(self):
        """ Data type for position numbers. """
        dtypes = list(set(array.dtype for array in (self.end5s,
                                                    self.mid5s,
                                                    self.mid3s,
                                                    self.end3s)))
        if len(dtypes) != 1:
            raise

    @cached_property
    @abstractmethod
    def array(self) -> np.ndarray:
        """ Array of data in use. """

    def get_index(self, refseq: DNA):
        """ Positions and bases in use. """
        if len(refseq) != self.max_pos:
            raise ValueError(f"Expected reference sequence of {self.max_pos} "
                             f"nt, but got {len(refseq)} nt")
        return seq_pos_to_index(refseq, self.pos, POS_INDEX)

    def get_table(self, refseq: DNA):
        """ Table of data in use. """
        return pd.DataFrame(self.array,
                            self.read_nums,
                            self.get_index(refseq),
                            copy=False)


class MutsBatch(ArrayBatch, ABC):

    def __init__(self,
                 *args,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.muts = {pos: {REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                           for rel, reads in muts.get(pos, dict()).items()}
                     for pos in self.pos}

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
