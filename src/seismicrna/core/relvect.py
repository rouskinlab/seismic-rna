import gzip
from collections import defaultdict
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .nptype import (ENDIAN_STR, get_uint_type, get_uint_max, get_uint_size,
                     get_byte_dtype)
from .sect import seq_pos_to_index
from .seq import DNA

logger = getLogger(__name__)

# Data types.
NAME_TYPE = str
MAGIC_TYPE = get_uint_type(4)
POS_TYPE = get_uint_type(4)
READ_TYPE = get_uint_type(4)
REL_TYPE = get_uint_type(1)
SIZE_TYPE = get_uint_type(4)

# Minimum and maximum values of each data type.
MIN_MAGIC = 0
MAX_MAGIC = get_uint_max(MAGIC_TYPE)
MIN_POS = 1
MAX_POS = get_uint_max(POS_TYPE)
MIN_READ = 0
MAX_READ = get_uint_max(READ_TYPE)
MIN_REL = 0
MAX_REL = get_uint_max(REL_TYPE)
MIN_SIZE = 0
MAX_SIZE = get_uint_max(SIZE_TYPE)

# Magic numbers.
MAGIC_A = MAGIC_TYPE(int.from_bytes(b"relA", ENDIAN_STR))
MAGIC_B = MAGIC_TYPE(int.from_bytes(b"relB", ENDIAN_STR))
MAGIC_C = MAGIC_TYPE(int.from_bytes(b"relC", ENDIAN_STR))

# Integer encodings for relation vectors
IRREC = REL_TYPE(int("00000000", 2))  # 000: irreconcilable mates
MATCH = REL_TYPE(int("00000001", 2))  # 001: match
DELET = REL_TYPE(int("00000010", 2))  # 002: deletion
INS_5 = REL_TYPE(int("00000100", 2))  # 004: 5' of insertion
INS_3 = REL_TYPE(int("00001000", 2))  # 008: 3' of insertion
SUB_A = REL_TYPE(int("00010000", 2))  # 016: substitution to A
SUB_C = REL_TYPE(int("00100000", 2))  # 032: substitution to C
SUB_G = REL_TYPE(int("01000000", 2))  # 064: substitution to G
SUB_T = REL_TYPE(int("10000000", 2))  # 128: substitution to T
NOCOV = REL_TYPE(int("11111111", 2))  # 255: not covered by read
SUB_N = REL_TYPE(SUB_A | SUB_C | SUB_G | SUB_T)
SUB_B = REL_TYPE(SUB_N ^ SUB_A)
SUB_D = REL_TYPE(SUB_N ^ SUB_C)
SUB_H = REL_TYPE(SUB_N ^ SUB_G)
SUB_V = REL_TYPE(SUB_N ^ SUB_T)
ANY_N = REL_TYPE(SUB_N | MATCH)
ANY_B = REL_TYPE(SUB_B | MATCH)
ANY_D = REL_TYPE(SUB_D | MATCH)
ANY_H = REL_TYPE(SUB_H | MATCH)
ANY_V = REL_TYPE(SUB_V | MATCH)
INS_8 = REL_TYPE(INS_5 | INS_3)
INDEL = REL_TYPE(DELET | INS_8)
MINS5 = REL_TYPE(INS_5 | MATCH)
MINS3 = REL_TYPE(INS_3 | MATCH)
ANY_8 = REL_TYPE(INS_8 | MATCH)


class RelationFileFormatError(Exception):
    """ Error during parsing of a relation format file. """


class RelationArray(object):
    """ """

    def __init__(self,
                 refseq: DNA,
                 names: list[str] | np.ndarray,
                 end5s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 muts: list[dict[int, list[int] | np.ndarray]]):
        self.refseq = refseq
        self.names = np.asarray(names, dtype=NAME_TYPE)
        self._end5s = np.asarray(end5s, dtype=POS_TYPE)
        self._end3s = np.asarray(end3s, dtype=POS_TYPE)
        self._muts = [{REL_TYPE(rel): np.asarray(reads, dtype=READ_TYPE)
                       for rel, reads in rels.items()}
                      for rels in muts]
        if np.any(self.end5s > self.end3s):
            errors = self.names[self.end5s > self.end3s]
            raise ValueError(f"Got end5 > end3 for reads {errors}")

    @cached_property
    def npos(self):
        npos = len(self.refseq)
        if npos > MAX_POS:
            raise ValueError(f"Reference sequence must be ≤ {MAX_POS} nt, "
                             f"but got {npos} nt")
        return npos

    @cached_property
    def pos(self):
        """ Array of 1-indexed positions. """
        return np.arange(MIN_POS, self.npos + MIN_POS)

    @cached_property
    def nreads(self):
        nreads = len(self.names)
        # The maximum integer label for a read is MAX_READ. Because the
        # first read is labeled 0, the maximum possible number of reads
        # such that the last one's label is ≤ MAX_READ is MAX_READ + 1.
        if nreads > MAX_READ + 1:
            raise ValueError(f"Number of reads must be ≤ {MAX_READ + 1}, "
                             f"but got {nreads}")
        return nreads

    @cached_property
    def names_bytes(self):
        """ Read names as bytes. """
        return np.char.encode(self.names)

    @cached_property
    def end5s(self):
        if self._end5s.shape != (self.nreads,):
            raise ValueError(
                f"Expected {self.nreads} 5' ends, but got {self._end5s.shape}")
        if self.nreads and np.min(self._end5s) < MIN_POS:
            errors = self.names[self._end5s < MIN_POS]
            raise ValueError(f"Got end5 < {MIN_POS} for reads {errors}")
        return self._end5s

    @cached_property
    def end3s(self):
        if self._end3s.shape != (self.nreads,):
            raise ValueError(
                f"Expected {self.nreads} 3' ends, but got {self._end3s.shape}")
        if self.nreads and np.max(self._end3s) > self.npos:
            errors = self.names[self._end3s > self.npos]
            raise ValueError(f"Got end3 > {self.npos} for reads {errors}")
        return self._end3s

    @cached_property
    def muts(self):
        if len(self._muts) != self.npos:
            raise ValueError(f"Got {len(self._muts)} mutation position but "
                             f"{self.npos} positions in the reference sequence")
        return self._muts

    @cached_property
    def columns(self):
        return seq_pos_to_index(self.refseq, self.pos, MIN_POS)

    @cached_property
    def df(self):
        """ DataFrame of relation vectors. """
        # Initialize an empty DataFrame.
        df = pd.DataFrame(NOCOV,
                          index=self.names,
                          columns=self.columns,
                          dtype=REL_TYPE)
        # Mark the covered positions as matches.
        for name, end5, end3 in zip(self.names,
                                    self.end5s,
                                    self.end3s,
                                    strict=True):
            df.loc[name, end5: end3] = MATCH
        # Mark the mutated positions.
        for pos, rels in enumerate(self.muts, start=MIN_POS):
            for rel, reads in rels.items():
                df.loc[self.names[reads], pos] = rel
        return df

    def save(self, file: Path, overwrite: bool = False):
        """  """
        with gzip.open(file, 'wb' if overwrite else 'xb') as f:
            # Magic number A (beginning of file).
            f.write(MAGIC_A.tobytes())
            # Number of reads.
            f.write(READ_TYPE(self.nreads).tobytes())
            # Size of each read name.
            f.write(SIZE_TYPE(self.names_bytes.itemsize).tobytes())
            # Read names.
            f.write(self.names_bytes.tobytes())
            # 5' positions.
            f.write(self.end5s.tobytes())
            # 3' positions.
            f.write(self.end3s.tobytes())
            # Magic number B (separating read and position data).
            f.write(MAGIC_B.tobytes())
            # Number of positions.
            f.write(POS_TYPE(self.npos).tobytes())
            # Iterate over the positions.
            for rels in self.muts:
                # Number of relationships for the position.
                f.write(REL_TYPE(len(rels)).tobytes())
                # Iterate over the relationships.
                for rel, reads in rels.items():
                    # Type of relationship.
                    f.write(rel.tobytes())
                    # Number of reads with the relationship.
                    f.write(READ_TYPE(reads.size).tobytes())
                    # Numbers of the reads with the relationship.
                    f.write(reads.tobytes())
            # Magic number C (end of file).
            f.write(MAGIC_C.tobytes())


def from_file(file: Path, refseq: DNA):
    """ """
    with gzip.open(file, 'rb') as f:

        def read_uint(uint_type: type):
            return uint_type(int.from_bytes(f.read(get_uint_size(uint_type)),
                                            ENDIAN_STR))

        def verify_magic_num(expect: int):
            if (magic := read_uint(MAGIC_TYPE)) != expect:
                raise RelationFileFormatError(
                    f"Bad magic number: expected {expect}, but got {magic}")

        def read_data(n_items: int,
                      item_size: int,
                      dtype_func: Callable = get_uint_type):
            items = np.frombuffer(f.read(n_items * item_size),
                                  dtype_func(item_size))
            if items.size != n_items:
                raise ValueError(
                    f"Expected {n_items} items, but got {items.size}")
            return items

        verify_magic_num(MAGIC_A)
        nreads = read_uint(READ_TYPE)
        names = np.char.decode(read_data(nreads,
                                         read_uint(SIZE_TYPE),
                                         get_byte_dtype))
        end5s = read_data(nreads, get_uint_size(POS_TYPE))
        end3s = read_data(nreads, get_uint_size(POS_TYPE))
        verify_magic_num(MAGIC_B)
        muts = [{read_uint(REL_TYPE): read_data(read_uint(READ_TYPE),
                                                get_uint_size(READ_TYPE))
                 for __ in range(read_uint(REL_TYPE))}
                for _ in range(read_uint(POS_TYPE))]
        verify_magic_num(MAGIC_C)
    return RelationArray(refseq, names, end5s, end3s, muts)


def from_reads(refseq: DNA,
               reads: Iterable[tuple[str, int, int, dict[int, int]]]):
    """ Accumulate reads into relation vectors. """
    npos = len(refseq)
    if not MIN_POS <= npos <= MAX_POS:
        raise ValueError(f"Reference sequence must be {MIN_POS}-{MAX_POS} nt, "
                         f"but got {len(refseq)} nt")
    names = list()
    end5s = list()
    end3s = list()
    muts = [defaultdict(list) for _ in refseq]
    read = MIN_READ
    for name, end5, end3, poss in reads:
        try:
            if not MIN_POS <= end5 <= npos:
                raise ValueError(
                    f"Invalid 5' position for sequence of {npos} nt: {end5}")
            if not MIN_POS <= end3 <= npos:
                raise ValueError(
                    f"Invalid 3' position for sequence of {npos} nt: {end3}")
            if end5 > end3:
                raise ValueError(f"5' position ({end5}) > 3' position ({end3})")
            if poss:
                if min(poss) < end5:
                    errors = {pos for pos in poss if pos < end5}
                    raise ValueError(f"Position(s) < 5' end ({end5}): {errors}")
                if max(poss) > end3:
                    errors = {pos for pos in poss if pos > end3}
                    raise ValueError(f"Position(s) > 5' end ({end3}): {errors}")
                if min(poss.values()) < MIN_REL:
                    errors = {rel for rel in poss.values() if rel < MIN_REL}
                    raise ValueError(f"Relationships(s) < {MIN_REL}: {errors}")
                if max(poss.values()) > MAX_REL:
                    errors = {rel for rel in poss.values() if rel > MAX_REL}
                    raise ValueError(f"Relationships(s) > {MAX_REL}: {errors}")
        except ValueError as error:
            logger.error(f"Failed to process read {repr(name)}: {error}")
        else:
            names.append(name)
            end5s.append(end5)
            end3s.append(end3)
            for pos, mut in poss.items():
                muts[pos - 1][mut].append(read)
            read += 1
    return RelationArray(refseq, names, end5s, end3s, muts)

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
