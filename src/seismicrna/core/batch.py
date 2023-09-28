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
from .sect import seq_pos_to_index
from .seq import DNA
from .types import fit_uint_type

logger = getLogger(__name__)

PICKLE_PROTOCOL = 5

BATCH_INDEX = 0
POS_INDEX = 1
READ_INDEX = 0


class Batch(Output, ABC):
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

    def __init__(self, sample: str, ref: str, refseq: DNA, batch: int):
        self.sample = sample
        self.ref = ref
        self.refseq = refseq
        self.batch = batch

    @property
    @abstractmethod
    def num_pos(self) -> int:
        """ Number of positions that are actually in use. """

    @property
    @abstractmethod
    def num_reads(self) -> int:
        """ Number of reads that are actually in use. """

    @property
    def max_pos(self):
        """ Maximum possible value for a position. """
        return len(self.refseq) + POS_INDEX - 1

    @property
    @abstractmethod
    def max_read(self) -> int:
        """ Maximum possible value for a read index. """

    @property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(self.max_pos)

    @property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(self.max_read)

    @cached_property
    @abstractmethod
    def pos(self) -> np.ndarray:
        """ Array of 1-indexed positions. """

    @cached_property
    def index(self):
        """ Positions and bases in use. """
        return seq_pos_to_index(self.refseq, self.pos, POS_INDEX)

    @cached_property
    @abstractmethod
    def read_nums(self) -> np.ndarray:
        """ Read numbers in use. """

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
        state = self.__dict__.copy()
        # Do not pickle cached properties.
        for name, value in list(state.items()):
            if isinstance(value, cached_property):
                state.pop(name)
        # Compress the reference sequence before pickling.
        state["refseq"] = self.refseq.compress()
        return state

    def __setstate__(self, state: dict[str, Any]):
        self.__dict__.update(state)
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
