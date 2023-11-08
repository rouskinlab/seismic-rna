from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .count import (get_count_per_pos,
                    get_count_per_read,
                    get_coverage_matrix,
                    get_cover_per_pos,
                    get_cover_per_read,
                    get_reads_per_pos,
                    get_rels_per_pos,
                    get_rels_per_read)
from .index import POS_INDEX, iter_windows, sanitize_ends, sanitize_pos
from .read import ReadBatch, PartialReadBatch
from ..rel import REL_TYPE, RelPattern
from ..seq import BASE_NAME, DNA, seq_pos_to_index
from ..types import fit_uint_type

logger = getLogger(__name__)


class MutsBatch(ReadBatch, ABC):
    """ Batch of mutational data. """

    def __init__(self, *,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 sanitize: bool = True,
                 **kwargs):
        super().__init__(**kwargs)
        (self.end5s,
         self._mid5s,
         self._mid3s,
         self.end3s) = (sanitize_ends(self.max_pos, end5s, mid5s, mid3s, end3s)
                        if sanitize
                        else (end5s, mid5s, mid3s, end3s))
        self.muts = {pos: {REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                           for rel, reads in muts[pos].items()}
                     for pos in (sanitize_pos(muts, self.max_pos)
                                 if sanitize
                                 else muts)}

    @property
    def num_pos(self):
        """ Number of positions in use. """
        return self.pos_nums.size

    @property
    @abstractmethod
    def max_pos(self) -> int:
        """ Maximum allowed position. """

    @cached_property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(self.max_pos)

    @cached_property
    @abstractmethod
    def pos_nums(self) -> np.ndarray:
        """ Positions in use. """

    @property
    def mid5s(self):
        return self._mid5s if self._mid5s is not None else self.end5s

    @property
    def mid3s(self):
        return self._mid3s if self._mid3s is not None else self.end3s

    @property
    def read_weights(self) -> pd.DataFrame | None:
        """ Weights for each read when computing counts. """
        return None


class ReflenMutsBatch(MutsBatch, ABC):
    """ Batch of mutational data with only a known reference length. """

    def __init__(self, *, reflen: int, **kwargs):
        self.seqlen = reflen
        super().__init__(**kwargs)

    @property
    def max_pos(self) -> int:
        return self.seqlen


class RefseqMutsBatch(MutsBatch, ABC):
    """ Batch of mutational data with a known reference sequence. """
    
    def __init__(self, *, refseq: DNA, **kwargs):
        super().__init__(**kwargs)
        self.refseq = refseq

    @property
    def max_pos(self):
        return len(self.refseq)

    @cached_property
    def pos_index(self):
        return seq_pos_to_index(self.refseq, self.pos_nums, POS_INDEX)

    @cached_property
    def count_base_types(self):
        bases = self.pos_index.get_level_values(BASE_NAME)
        base_types, counts = np.unique(bases, return_counts=True)
        return pd.Series(counts, base_types)

    def iter_windows(self, size: int):
        """ Yield the positions in each window of size positions of the
        section. """
        yield from iter_windows(self.pos_nums, size)

    @cached_property
    def coverage_matrix(self):
        return get_coverage_matrix(self.pos_index,
                                   self.end5s,
                                   self.mid5s,
                                   self.mid3s,
                                   self.end3s,
                                   self.read_nums)

    @cached_property
    def cover_per_pos(self):
        """ Number of reads covering each position. """
        return get_cover_per_pos(self.coverage_matrix, self.read_weights)

    @cached_property
    def cover_per_read(self):
        """ Number of positions covered by each read. """
        return get_cover_per_read(self.coverage_matrix)

    @cached_property
    def rels_per_pos(self):
        """ For each relationship, the number of reads at each position
        with that relationship. """
        return get_rels_per_pos(self.muts,
                                self.num_reads,
                                self.cover_per_pos,
                                self.read_indexes,
                                self.read_weights)

    @cached_property
    def rels_per_read(self):
        """ For each relationship, the number of positions in each read
        with that relationship. """
        return get_rels_per_read(self.muts,
                                 self.pos_index,
                                 self.cover_per_read,
                                 self.read_indexes)

    def reads_per_pos(self, pattern: RelPattern):
        """ For each position, find all reads matching a relationship
        pattern. """
        return get_reads_per_pos(pattern, self.muts, self.pos_index)

    def count_per_pos(self, pattern: RelPattern):
        """ Count the reads that fit a relationship pattern at each
        position in a section. """
        return get_count_per_pos(pattern,
                                 self.cover_per_pos,
                                 self.rels_per_pos)

    def count_per_read(self, pattern: RelPattern):
        """ Count the positions in a section that fit a relationship
        pattern in each read. """
        return get_count_per_read(pattern,
                                  self.cover_per_read,
                                  self.rels_per_read,
                                  self.read_weights)

    def nonprox_muts(self, pattern: RelPattern, min_gap: int):
        """ List the reads with non-proximal mutations. """
        if min_gap < 0:
            raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
        if min_gap == 0:
            # No reads can have proximal mutations.
            return self.read_nums
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(pattern)
        # Track which reads have no proximal mutations.
        nonprox = np.ones(self.num_reads, dtype=bool)
        # Count the number of times each read occurs in a window.
        read_counts = np.zeros_like(nonprox, dtype=self.pos_dtype)
        # Track which positions are being counted.
        pos_counted = set()
        # Iterate over all windows in the section.
        for _, win_pos in self.iter_windows(min_gap + 1):
            win_pos_set = set(win_pos)
            for pos in win_pos_set - pos_counted:
                # These positions have entered the window: count them.
                read_counts[self.read_indexes[reads_per_pos[pos]]] += 1
                pos_counted.add(pos)
            for pos in pos_counted - win_pos_set:
                # These positions have exited the window: drop them.
                read_counts[self.read_indexes[reads_per_pos[pos]]] -= 1
                pos_counted.remove(pos)
            if len(pos_counted) > 1:
                # Check for reads counted more than once in the window,
                # which means that they have proximal mutations.
                nonprox &= read_counts <= 1
        return self.read_nums[nonprox]

    def iter_reads(self, pattern: RelPattern):
        """ Yield the 5'/3' end/middle positions and the positions that
        are mutated in each read. """
        reads_per_pos = self.reads_per_pos(pattern)
        # Find the maximum number of mutations in any read.
        max_mut_count = max(sum(map(Counter, reads_per_pos.values()),
                                start=Counter()).values(),
                            default=0)
        # Initialize a matrix of the positions mutated in each read.
        mut_pos = np.zeros((self.num_reads, max_mut_count),
                           dtype=self.pos_dtype)
        # Fill in the matrix one position at a time.
        for pos, reads in reads_per_pos.items():
            # Append the position to the end of each row corresponding
            # to a read that has a mutation at the position.
            row_idxs = self.read_indexes[reads]
            mut_pos[row_idxs, np.count_nonzero(mut_pos[row_idxs], axis=1)] = pos
        # Count the mutations in each read.
        mut_counts = np.count_nonzero(mut_pos, axis=1)
        # For each read, yield the 5'/3' end/middle positions and the
        # mutated positions as a tuple.
        for read_num, end5, mid5, mid3, end3, muts, count in zip(self.read_nums,
                                                                 self.end5s,
                                                                 self.mid5s,
                                                                 self.mid3s,
                                                                 self.end3s,
                                                                 mut_pos,
                                                                 mut_counts,
                                                                 strict=True):
            yield read_num, ((end5, mid5, mid3, end3), tuple(muts[:count]))


class PartialMutsBatch(PartialReadBatch, MutsBatch, ABC):

    @cached_property
    def pos_nums(self):
        return np.array(list(self.muts))

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
