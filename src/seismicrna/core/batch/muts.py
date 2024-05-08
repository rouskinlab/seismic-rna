from abc import ABC
from collections import Counter
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .count import (count_end_coords,
                    get_count_per_pos,
                    get_count_per_read,
                    calc_coverage,
                    get_reads_per_pos,
                    get_rels_per_pos,
                    get_rels_per_read)
from .index import (iter_windows,
                    sanitize_ends,
                    split_ends)
from .read import ReadBatch, PartialReadBatch
from ..rel import REL_TYPE, RelPattern
from ..seq import Section

logger = getLogger(__name__)


class MutsBatch(ReadBatch, ABC):
    """ Batch of mutational data. """

    def __init__(self, *,
                 section: Section,
                 ends: np.ndarray,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 sanitize: bool = True,
                 **kwargs):
        super().__init__(**kwargs)
        # Validate and store the read end coordinates.
        self.ends = sanitize_ends(ends, section.end5, section.end3, sanitize)
        # Validate and store the mutations.
        self.muts = {pos: ({REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                            for rel, reads in muts[pos].items()}
                           if sanitize
                           else muts[pos])
                     for pos in section.unmasked_int}

    @cached_property
    def end5s(self):
        """ 5' end coordinates of reads. """
        end5s, _ = split_ends(self.ends)
        return end5s.min(axis=1)

    @cached_property
    def end3s(self):
        """ 3' end coordinates of reads. """
        _, end3s = split_ends(self.ends)
        return end3s.max(axis=1)

    @property
    def pos_dtype(self):
        """ Data type for positions. """
        return self.ends.dtype

    @cached_property
    def pos_nums(self):
        """ Positions in use. """
        return np.fromiter(self.muts, self.pos_dtype)

    @property
    def read_weights(self) -> pd.DataFrame | None:
        """ Weights for each read when computing counts. """
        return None

    @cached_property
    def end_counts(self):
        """ Counts of end coordinates. """
        return count_end_coords(self.end5s, self.end3s, self.read_weights)


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
        for read_num, end5, end3, muts, count in zip(self.read_nums,
                                                     self.end5s,
                                                     self.end3s,
                                                     mut_pos,
                                                     mut_counts,
                                                     strict=True):
            yield read_num, ((end5, end3), tuple(muts[:count]))


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
