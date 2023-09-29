import pickle
import re
from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Any

import brotli
import numpy as np
import pandas as pd

from .cli import opt_brotli_level
from .output import Output
from .rel import REL_TYPE
from .sect import seq_pos_to_index
from .seq import DNA
from .types import fit_uint_type

logger = getLogger(__name__)

PICKLE_PROTOCOL = 5

BATCH_INDEX = 0
POS_INDEX = 1
READ_INDEX = 0

READ_NUM = "Read Number"


class BaseBatch(Output, ABC):
    """ """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z]*)batch$", cls.__name__.lower()).groups()
        return btype

    @classmethod
    def load(cls, file: Path):
        """ Load from a compressed pickle file. """
        with open(file, 'rb') as f:
            obj = pickle.loads(brotli.decompress(f.read()))
        if not isinstance(obj, cls):
            raise TypeError(
                f"Expected {cls.__name__}, but got {obj.__class__.__name__}")
        return obj

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

    def save(self,
             top: Path,
             brotli_level: int = opt_brotli_level.default,
             overwrite: bool = False):
        """ Save to a pickle file compressed with Brotli. """
        data = brotli.compress(pickle.dumps(self, PICKLE_PROTOCOL),
                               quality=brotli_level)
        save_path = self.get_path(top)
        with open(save_path, 'wb' if overwrite else 'xb') as f:
            f.write(data)
        logger.debug(f"Wrote {self} to {save_path}")
        return save_path

    def __getstate__(self):
        # Copy the __dict__ to avoid modifying this object's state.
        state = self.__dict__.copy()
        # Do not pickle cached properties.
        for name, value in list(state.items()):
            if isinstance(value, cached_property):
                state.pop(name)
        return state

    def __setstate__(self, state: dict[str, Any]):
        self.__dict__.update(state)


def ensure_1d(array: np.ndarray, what: str):
    if array.ndim != 1:
        raise ValueError(f"{what} must have 1 dimension, but got {array.ndim}")


def ensure_same_dims(array1: np.ndarray,
                     array2: np.ndarray,
                     what1: str,
                     what2: str):
    if array1.shape != array2.shape:
        raise ValueError(f"Shapes differ between {what1} {array1.shape} "
                         f"and {what2} {array2.shape}")


def ensure_order(read_nums: np.ndarray,
                 array1: np.ndarray,
                 array2: np.ndarray,
                 what1: str,
                 what2: str,
                 gt_eq: bool = False):
    ensure_same_dims(array1, array2, what1, what2)
    ensure_same_dims(array1, read_nums, what1, "read numbers")
    ineq_func, ineq_sign = (np.less, '<') if gt_eq else (np.greater, '>')
    if np.any(is_err := ineq_func(array1, array2)):
        index = pd.Index(read_nums[is_err], name=READ_NUM)
        errors = pd.DataFrame.from_dict({what1: pd.Series(array1[is_err],
                                                          index=index),
                                         what2: pd.Series(array2[is_err],
                                                          index=index)})
        raise ValueError(f"Got {what1} {ineq_sign} {what2}:\n{errors}")


class ArrayBatch(BaseBatch, ABC):

    def __init__(self,
                 *args,
                 refseq: DNA,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.refseq = refseq
        # 5' end positions
        self.end5s = np.asarray(end5s, self.pos_dtype)
        ensure_1d(self.end5s, "5' end positions")
        # 5' middle positions
        self.mid5s = np.asarray(mid5s, self.pos_dtype)
        ensure_1d(self.mid5s, "5' middle positions")
        # 3' middle positions
        self.mid3s = np.asarray(mid3s, self.pos_dtype)
        ensure_1d(self.mid3s, "3' middle positions")
        # 3' end positions
        self.end3s = np.asarray(end3s, self.pos_dtype)
        ensure_1d(self.end3s, "3' end positions")
        # Verify 5' end ≥ min position
        ensure_order(self.read_nums,
                     self.end5s,
                     np.broadcast_to(POS_INDEX, self.end5s.shape),
                     "5' end positions",
                     f"minimum position ({POS_INDEX})",
                     gt_eq=True)
        # Verify 5' end ≤ 5' mid
        ensure_order(self.read_nums,
                     self.end5s,
                     self.mid5s,
                     "5' end positions",
                     "5' middle positions")
        # Verify 5' end ≤ 3' mid
        ensure_order(self.read_nums,
                     self.end5s,
                     self.mid3s,
                     "5' end positions",
                     "3' middle positions")
        # Verify 5' mid ≤ 3' end
        ensure_order(self.read_nums,
                     self.mid5s,
                     self.end3s,
                     "5' middle positions",
                     "3' end positions")
        # Verify 3' mid ≤ 3' end
        ensure_order(self.read_nums,
                     self.mid3s,
                     self.end3s,
                     "3' middle positions",
                     "3' end positions")
        # Verify 3' end ≤ max position
        ensure_order(self.read_nums,
                     self.end3s,
                     np.broadcast_to(self.max_pos, self.end3s.shape),
                     "3' end positions",
                     f"maximum position ({self.max_pos})")
        # Mutations
        self.muts = {pos: {REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                           for rel, reads in muts.get(pos, dict()).items()}
                     for pos in self.pos}

    @cached_property
    @abstractmethod
    def pos(self) -> np.ndarray:
        """ Array of 1-indexed positions. """

    @property
    @abstractmethod
    def num_pos(self) -> int:
        """ Number of positions that are actually in use. """

    @property
    def max_pos(self):
        """ Maximum possible value for a position. """
        return len(self.refseq) + POS_INDEX - 1

    @property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(self.max_pos)

    @cached_property
    def index(self):
        """ Positions and bases in use. """
        return seq_pos_to_index(self.refseq, self.pos, POS_INDEX)

    @cached_property
    @abstractmethod
    def array(self) -> np.ndarray:
        """ Array of data in use. """

    @cached_property
    def table(self) -> pd.DataFrame | pd.Series:
        """ Table of data in use. """
        return pd.DataFrame(self.array,
                            self.read_nums,
                            self.index,
                            copy=False)

    def __getstate__(self):
        state = super().__getstate__()
        # Compress the reference sequence before pickling.
        state["refseq"] = self.refseq.compress()
        return state

    def __setstate__(self, state: dict[str, Any]):
        super().__setstate__(state)
        # Decompress the reference sequence when unpickling.
        self.refseq = state["refseq"].decompress()

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
