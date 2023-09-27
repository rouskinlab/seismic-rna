import pickle
from collections import defaultdict
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Iterable

import brotli
import numpy as np
import pandas as pd

from .rel import MATCH, NOCOV, REL_TYPE
from .sect import seq_pos_to_index
from .seq import DNA
from .types import fit_uint_type

logger = getLogger(__name__)

PICKLE_PROTOCOL = 5

BATCH_INDEX = 0
POS_INDEX = 1
READ_INDEX = 0


class ArrayBatch(object):
    """ """

    NO_PICKLE = "pos", "index", "names", "matrix", "df"

    def __init__(self,
                 batch: int,
                 refseq: DNA,
                 names: list[str] | np.ndarray,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 muts: dict[int, dict[int, list[int] | np.ndarray]]):
        self.batch = batch
        self.refseq = refseq
        # Encode read names as bytes, which occupy 1/4 of the memory in
        # NumPy arrays compared to unicode characters.
        self._names = np.char.encode(names)
        # 5' ends
        self.end5s = np.asarray(end5s, dtype=self.pos_dtype)
        if self.end5s.shape != (self.num_reads,):
            raise ValueError(
                f"Expected {self.num_reads} end5s, but got {self.end5s.shape}")
        if np.any(self.end5s < POS_INDEX):
            errors = self.names[self.end5s < POS_INDEX]
            raise ValueError(f"Got end5 < {POS_INDEX} for reads {errors}")
        # 3' ends
        self.end3s = np.asarray(end3s, dtype=self.pos_dtype)
        if self.end3s.shape != (self.num_reads,):
            raise ValueError(
                f"Expected {self.num_reads} 3' ends, but got {self.end3s.shape}")
        if np.any(self.end3s < self.end5s):
            errors = self.names[self.end3s < self.end5s]
            raise ValueError(f"Got end3 < end5 for reads {errors}")
        # 5' middles
        self.mid5s = np.asarray(mid5s, dtype=self.pos_dtype)
        if self.mid5s.shape != (self.num_reads,):
            raise ValueError(
                f"Expected {self.num_reads} mid5s, but got {self.mid5s.shape}")
        if np.any(self.mid5s < self.end5s):
            errors = self.names[self.mid5s < self.end5s]
            raise ValueError(f"Got mid5 < end5 for reads {errors}")
        # 3' middles
        self.mid3s = np.asarray(mid3s, dtype=self.pos_dtype)
        if self.mid3s.shape != (self.num_reads,):
            raise ValueError(
                f"Expected {self.num_reads} mid3s, but got {self.mid3s.shape}")
        if np.any(self.mid3s > self.end3s):
            errors = self.names[self.mid3s > self.end3s]
            raise ValueError(f"Got mid3 > end3 for reads {errors}")
        # Mutations
        self.muts = dict(((self.pos_dtype(pos),
                           {REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                            for rel, reads in rels.items()})
                          for pos, rels in muts.items()))

    @property
    def num_pos(self):
        """ Number of positions that are actually in use. """
        return len(self.muts)

    @property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(len(self.refseq))

    @property
    def num_reads(self):
        """ Number of reads that are actually in use. """
        return len(self._names)

    @property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(self.num_reads)

    @cached_property
    def pos(self):
        """ Array of 1-indexed positions. """
        return np.array(list(self.muts))

    @cached_property
    def index(self):
        """ Positions and bases in use. """
        return seq_pos_to_index(self.refseq, self.pos, POS_INDEX)

    @cached_property
    def names(self):
        """ Read names in use. """
        return np.char.decode(self._names)

    @cached_property
    def matrix(self):
        """ Matrix of vectors in use. """
        matrix = np.full((self.num_reads, self.num_pos), NOCOV)
        # Mark the covered positions as matches.
        for read, (end5i, mid5i, mid3p, end3p) in enumerate(zip(self.end5s - 1,
                                                                self.mid5s - 1,
                                                                self.mid3s,
                                                                self.end3s,
                                                                strict=True)):
            matrix[read, end5i: mid3p] = MATCH
            matrix[read, mid5i: end3p] = MATCH
        # Mark the mutated positions.
        for pos, rels in self.muts.items():
            for rel, reads in rels.items():
                matrix[reads, pos - 1] = rel
        return matrix

    @cached_property
    def df(self):
        """ DataFrame of vectors in use. """
        return pd.DataFrame(self.matrix, self.names, self.index, copy=False)

    def save(self, file: Path, brotli_level: int, overwrite: bool = False):
        """ Save to a compressed pickle file. """
        pickled = pickle.dumps(self, PICKLE_PROTOCOL)
        with open(file, 'wb' if overwrite else 'xb') as f:
            f.write(brotli.compress(pickled, quality=brotli_level))

    @classmethod
    def load(cls, file: Path):
        """ Load from a compressed pickle file. """
        with open(file, 'rb') as f:
            obj = pickle.loads(brotli.decompress(f.read()))
        if not isinstance(obj, cls):
            raise TypeError(
                f"Expected {cls.__name__}, but got {obj.__class__.__name__}")
        return obj

    def __getstate__(self):
        state = self.__dict__.copy()
        for attr in self.NO_PICKLE:
            state.pop(attr, None)
        return state


def from_reads(reads: Iterable[tuple[str, int, int, int, int, dict[int, int]]],
               batch: int,
               refseq: DNA):
    """ Accumulate reads into relation vectors. """
    names = list()
    end5s = list()
    mid5s = list()
    mid3s = list()
    end3s = list()
    muts = {pos: defaultdict(list) for pos in range(POS_INDEX,
                                                    len(refseq) + POS_INDEX)}
    read = READ_INDEX
    for name, end5, mid5, mid3, end3, poss in reads:
        names.append(name)
        end5s.append(end5)
        mid5s.append(mid5)
        mid3s.append(mid3)
        end3s.append(end3)
        for pos, mut in poss.items():
            muts[pos][mut].append(read)
        read += 1
    return ArrayBatch(batch, refseq, names, end5s, mid5s, mid3s, end3s, muts)

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
