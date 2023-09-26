import gzip
from collections import defaultdict
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .sect import seq_pos_to_index
from .seq import DNA
from .types import (fit_uint_type, get_uint_type, get_max_uint, get_uint_size,
                    get_byte_dtype, ENDIAN_STR, UINT_NBYTES)

logger = getLogger(__name__)

# Data types.
MAGIC_SIZE = 4
MAGIC_TYPE = get_uint_type(MAGIC_SIZE)
REL_SIZE = 1
REL_TYPE = get_uint_type(REL_SIZE)
SIZE_TYPE = fit_uint_type(max(UINT_NBYTES))

# Minimum and maximum values of each data type.
MIN_POS = 1
MAX_POS = get_max_uint(get_uint_type(4))
MIN_READ = 0
MAX_READ = get_max_uint(get_uint_type(4))
MIN_REL = 0
MAX_REL = get_max_uint(REL_TYPE)

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
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 muts: list[dict[int, list[int] | np.ndarray]]):
        self.refseq = refseq
        self.names = np.asarray(names, dtype=str)
        self._end5s = np.asarray(end5s, dtype=self.pos_dtype)
        self._mid5s = np.asarray(mid5s, dtype=self.pos_dtype)
        self._mid3s = np.asarray(mid3s, dtype=self.pos_dtype)
        self._end3s = np.asarray(end3s, dtype=self.pos_dtype)
        self._muts = [{REL_TYPE(rel): np.asarray(reads, dtype=self.read_dtype)
                       for rel, reads in rels.items()}
                      for rels in muts]
        if np.any(self.end5s > self.end3s):
            errors = self.names[self.end5s > self.end3s]
            raise ValueError(f"Got end5 > end3 for reads {errors}")

    @property
    def npos(self):
        return len(self.refseq)

    @cached_property
    def pos(self):
        """ Array of 1-indexed positions. """
        return np.arange(MIN_POS, self.npos + MIN_POS)

    @cached_property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(self.npos)

    @property
    def nreads(self):
        return len(self.names)

    @cached_property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(self.nreads)

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
    def mid5s(self):
        if self._mid5s.shape != (self.nreads,):
            raise ValueError(
                f"Expected {self.nreads} 5' mids, but got {self._mid5s.shape}")
        if np.any(self._mid5s < self.end5s):
            errors = self.names[self._mid5s < self.end5s]
            raise ValueError(f"Got mid5 < end5 for reads {errors}")
        return self._mid5s

    @cached_property
    def mid3s(self):
        if self._mid3s.shape != (self.nreads,):
            raise ValueError(
                f"Expected {self.nreads} 3' mids, but got {self._mid3s.shape}")
        if np.any(self._mid3s > self.end3s):
            errors = self.names[self._mid3s > self.end3s]
            raise ValueError(f"Got mid3 > end3 for reads {errors}")
        return self._mid3s

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
        values = df.values
        # Mark the covered positions as matches.
        for read, (end5, mid5, mid3, end3) in enumerate(zip(self.end5s - 1,
                                                            self.mid5s - 1,
                                                            self.mid3s,
                                                            self.end3s,
                                                            strict=True)):
            values[read, end5: mid3] = MATCH
            values[read, mid5: end3] = MATCH
        # Mark the mutated positions.
        for idx, rels in enumerate(self.muts):
            for rel, reads in rels.items():
                values[reads, idx] = rel
        return df

    def save(self, file: Path, overwrite: bool = False):
        """  """
        with gzip.open(file, 'wb' if overwrite else 'xb') as f:
            # Magic number A (beginning of file).
            f.write(MAGIC_A.tobytes())
            # Size of each read name.
            names_bytes = np.char.encode(self.names)
            f.write(SIZE_TYPE(names_bytes.itemsize).tobytes())
            # Size of each read number.
            f.write(SIZE_TYPE(get_uint_size(self.read_dtype)).tobytes())
            # Size of each position number.
            f.write(SIZE_TYPE(get_uint_size(self.pos_dtype)).tobytes())
            # Number of reads.
            f.write(self.read_dtype(self.nreads).tobytes())
            # Read names.
            f.write(names_bytes.tobytes())
            # 5' end positions.
            f.write(self.end5s.tobytes())
            # 5' mid positions.
            f.write(self.mid5s.tobytes())
            # 3' mid positions.
            f.write(self.mid3s.tobytes())
            # 3' end positions.
            f.write(self.end3s.tobytes())
            # Magic number B (dividing per-read and per-position data).
            f.write(MAGIC_B.tobytes())
            # Number of positions.
            f.write(self.pos_dtype(self.npos).tobytes())
            # Iterate over the positions.
            for rels in self.muts:
                # Number of relationships for the position.
                f.write(REL_TYPE(len(rels)).tobytes())
                # Iterate over the relationships.
                for rel, reads in rels.items():
                    # Type of relationship.
                    f.write(rel.tobytes())
                    # Number of reads with the relationship.
                    f.write(self.read_dtype(reads.size).tobytes())
                    # Numbers of the reads with the relationship.
                    f.write(reads.tobytes())
            # Magic number C (end of file).
            f.write(MAGIC_C.tobytes())


def from_file(file: Path, refseq: DNA):
    """ """
    with gzip.open(file, 'rb') as f:

        def read_uint(uint_type: type):
            return int.from_bytes(f.read(get_uint_size(uint_type)), ENDIAN_STR)

        def verify_magic_num(expect: int):
            if (magic := read_uint(MAGIC_TYPE)) != expect:
                raise RelationFileFormatError(
                    f"Bad magic number: expected {expect}, but got {magic}")

        def _read_data(dtype: type, item_size: int, n_items: int):
            items = np.frombuffer(f.read(n_items * item_size), dtype)
            if items.size != n_items:
                raise ValueError(
                    f"Expected {n_items} items, but got {items.size}")
            return items

        def read_uint_data(uint_type: type, n_items: int):
            return _read_data(uint_type, get_uint_size(uint_type), n_items)

        def read_byte_data(size: int, n_items: int):
            return _read_data(get_byte_dtype(size), size, n_items)

        # Magic number A (beginning of file).
        verify_magic_num(MAGIC_A)
        # Size of each read name.
        name_size = read_uint(SIZE_TYPE)
        # Type of each read number.
        read_type = get_uint_type(read_uint(SIZE_TYPE))
        # Type of each position number.
        pos_type = get_uint_type(read_uint(SIZE_TYPE))
        # Number of reads.
        nreads = read_uint(read_type)
        # Read names.
        names = np.char.decode(read_byte_data(name_size, nreads))
        # 5' end positions.
        end5s = read_uint_data(pos_type, nreads)
        # 5' mid positions.
        mid5s = read_uint_data(pos_type, nreads)
        # 3' mid positions.
        mid3s = read_uint_data(pos_type, nreads)
        # 3' end positions.
        end3s = read_uint_data(pos_type, nreads)
        # Magic number B (separating read and position data).
        verify_magic_num(MAGIC_B)
        # Reads with each type of relationship at each position.
        muts = [{read_uint(REL_TYPE): read_uint_data(read_type,
                                                     read_uint(read_type))
                 for __ in range(read_uint(REL_TYPE))}
                for _ in range(read_uint(pos_type))]
        # Magic number C (end of file).
        verify_magic_num(MAGIC_C)
    return RelationArray(refseq, names, end5s, mid5s, mid3s, end3s, muts)


def from_reads(refseq: DNA,
               reads: Iterable[tuple[str, int, int, int, int, dict[int, int]]]):
    """ Accumulate reads into relation vectors. """
    npos = len(refseq)
    if not MIN_POS <= npos <= MAX_POS:
        raise ValueError(f"Reference sequence must be {MIN_POS}-{MAX_POS} nt, "
                         f"but got {len(refseq)} nt")
    names = list()
    end5s = list()
    mid5s = list()
    mid3s = list()
    end3s = list()
    muts = [defaultdict(list) for _ in refseq]
    read = MIN_READ
    for name, end5, mid5, mid3, end3, poss in reads:
        try:
            if not MIN_POS <= end5 <= mid5 <= npos:
                raise ValueError(f"Invalid 5' end/mid for {npos} nt sequence: "
                                 f"{end5}/{mid5}")
            if not MIN_POS <= mid3 <= end3 <= npos:
                raise ValueError(f"Invalid 3' mid/end for {npos} nt sequence: "
                                 f"{mid3}/{end3}")
            if end5 > end3:
                raise ValueError(f"5' position ({end5}) > 3' position ({end3})")
            if poss:
                if min(poss) < end5:
                    errors = {pos for pos in poss if pos < end5}
                    raise ValueError(f"Position(s) < 5' end ({end5}): {errors}")
                if max(poss) > end3:
                    errors = {pos for pos in poss if pos > end3}
                    raise ValueError(f"Position(s) > 3' end ({end3}): {errors}")
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
            mid5s.append(mid5)
            mid3s.append(mid3)
            end3s.append(end3)
            for pos, mut in poss.items():
                muts[pos - 1][mut].append(read)
            read += 1
    return RelationArray(refseq, names, end5s, mid5s, mid3s, end3s, muts)

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
