from abc import ABC
from collections import Counter
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .count import (count_end_coords,
                    calc_count_per_pos,
                    calc_count_per_read,
                    calc_coverage,
                    calc_reads_per_pos,
                    calc_rels_per_pos,
                    calc_rels_per_read)
from .ends import (find_contiguous_reads,
                   find_read_end5s,
                   find_read_end3s,
                   mask_segment_ends,
                   match_reads_segments,
                   sanitize_segment_ends)
from .index import iter_windows
from .read import ReadBatch, PartialReadBatch
from ..rel import REL_TYPE, RelPattern
from ..seq import Section

logger = getLogger(__name__)


class MutsBatch(ReadBatch, ABC):
    """ Batch of mutational data. """

    def __init__(self, *,
                 section: Section,
                 seg_end5s: np.ndarray,
                 seg_end3s: np.ndarray,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 sanitize: bool = True,
                 **kwargs):
        super().__init__(**kwargs)
        # Validate and store the read end coordinates.
        self._end5s, self._end3s = sanitize_segment_ends(seg_end5s,
                                                         seg_end3s,
                                                         section.end5,
                                                         section.end3,
                                                         sanitize)
        # Validate and store the mutations.
        self._muts = {pos: ({REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                             for rel, reads in muts[pos].items()}
                            if sanitize
                            else muts[pos])
                      for pos in section.unmasked_int}

    @property
    def muts(self):
        """ Reads with each type of mutation at each position. """
        return self._muts

    @cached_property
    def num_segments(self):
        """ Number of segments in each read. """
        _, num_segs = match_reads_segments(self._end5s, self._end3s)
        return num_segs

    @cached_property
    def _seg_ends(self):
        """ 5' and 3' ends of each segment in each read. """
        return mask_segment_ends(self._end5s, self._end3s)

    @property
    def seg_end5s(self):
        """ 5' end of each segment in each read. """
        seg_end5s, _ = self._seg_ends
        return seg_end5s

    @property
    def seg_end3s(self):
        """ 3' end of each segment in each read. """
        _, seg_end3s = self._seg_ends
        return seg_end3s

    @cached_property
    def read_end5s(self):
        """ 5' end of each read. """
        return find_read_end5s(self.seg_end5s)

    @cached_property
    def read_end3s(self):
        """ 3' end of each read. """
        return find_read_end3s(self.seg_end3s)

    @property
    def pos_dtype(self):
        """ Data type for positions. """
        if self._end5s.dtype is not self._end3s.dtype:
            raise ValueError("Data types differ for 5' and 3' ends")
        return self._end5s.dtype

    @cached_property
    def pos_nums(self):
        """ Positions in use. """
        return np.fromiter(self.muts, self.pos_dtype)

    @property
    def read_weights(self) -> pd.DataFrame | None:
        """ Weights for each read when computing counts. """
        return None

    @cached_property
    def read_end_counts(self):
        """ Counts of read end coordinates. """
        return count_end_coords(self.read_end5s,
                                self.read_end3s,
                                self.read_weights)

    @cached_property
    def contiguous(self):
        """ Whether the segments of each read are contiguous. """
        return find_contiguous_reads(self.seg_end5s, self.seg_end3s)

    @cached_property
    def num_contiguous(self):
        """ Number of contiguous reads. """
        return np.count_nonzero(self.contiguous)

    @cached_property
    def num_discontiguous(self):
        """ Number of discontiguous reads. """
        return self.num_reads - self.num_contiguous


class SectionMutsBatch(MutsBatch, ABC):
    """ Batch of mutational data that knows its section. """

    def __init__(self, *, section: Section, **kwargs):
        self.section = section
        super().__init__(section=section, **kwargs)

    def iter_windows(self, size: int):
        """ Yield the positions in each window of size positions of the
        section. """
        yield from iter_windows(self.pos_nums, size)

    @cached_property
    def pos_index(self):
        """ Index of unmasked positions and bases. """
        return self.section.unmasked

    @cached_property
    def _coverage(self):
        """ Coverage per position and per read. """
        return calc_coverage(self.pos_index,
                             self.read_nums,
                             self.ends,
                             self.read_weights)

    @property
    def cover_per_pos(self):
        """ Number of reads covering each position. """
        per_pos, per_read = self._coverage
        return per_pos

    @property
    def cover_per_read(self):
        """ Number of positions covered by each read. """
        per_pos, per_read = self._coverage
        return per_read

    @cached_property
    def rels_per_pos(self):
        """ For each relationship, the number of reads at each position
        with that relationship. """
        return calc_rels_per_pos(self.muts,
                                 self.num_reads,
                                 self.cover_per_pos,
                                 self.read_indexes,
                                 self.read_weights)

    @cached_property
    def rels_per_read(self):
        """ For each relationship, the number of positions in each read
        with that relationship. """
        return calc_rels_per_read(self.muts,
                                  self.pos_index,
                                  self.cover_per_read,
                                  self.read_indexes)

    def reads_per_pos(self, pattern: RelPattern):
        """ For each position, find all reads matching a relationship
        pattern. """
        return calc_reads_per_pos(pattern, self.muts, self.pos_index)

    def count_per_pos(self, pattern: RelPattern):
        """ Count the reads that fit a relationship pattern at each
        position in a section. """
        return calc_count_per_pos(pattern,
                                  self.cover_per_pos,
                                  self.rels_per_pos)

    def count_per_read(self, pattern: RelPattern):
        """ Count the positions in a section that fit a relationship
        pattern in each read. """
        return calc_count_per_read(pattern,
                                   self.cover_per_read,
                                   self.rels_per_read,
                                   self.read_weights)

    def reads_noclose_muts(self, pattern: RelPattern, min_gap: int):
        """ List the reads with no two mutations too close. """
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

    def iter_reads(self,
                   pattern: RelPattern,
                   only_read_ends: bool = False,
                   require_contiguous: bool = False):
        """ End coordinates and mutated positions in each read. """
        if require_contiguous and self.num_discontiguous:
            raise ValueError("This function requires contiguous reads, but got "
                             f"{self.num_discontiguous} discontiguous read(s)")
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
        # For each read, yield the end coordinates and mutated positions
        # as a tuple.
        iter_ends = (zip(self.end5s, self.end3s, strict=True)
                     if only_read_ends or require_contiguous
                     else self.ends)
        for read_num, ends, muts, count in zip(self.read_nums,
                                               iter_ends,
                                               mut_pos,
                                               mut_counts,
                                               strict=True):
            yield read_num, (tuple(ends), tuple(muts[:count]))


class PartialMutsBatch(PartialReadBatch, MutsBatch, ABC):

    @cached_property
    def pos_nums(self):
        return np.array(list(self.muts))

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
