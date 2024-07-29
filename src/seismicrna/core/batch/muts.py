from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from numba import jit

from .count import (calc_count_per_pos,
                    calc_count_per_read,
                    calc_coverage,
                    calc_reads_per_pos,
                    calc_rels_per_pos,
                    calc_rels_per_read,
                    count_end_coords)
from .ends import EndCoords, match_reads_segments
from .index import iter_windows
from .read import ReadBatch, PartialReadBatch
from ..array import calc_inverse, find_dims
from ..header import REL_NAME
from ..rel import MATCH, NOCOV, REL_TYPE, RelPattern
from ..seq import Section, index_to_pos
from ..types import fit_uint_type

rng = np.random.default_rng()

logger = getLogger(__name__)

NUM_READS = "reads"
NUM_SEGMENTS = "segments"


def sanitize_muts(muts: dict[int, dict[int, list[int] | np.ndarray]],
                  section: Section,
                  data_type: type,
                  sanitize: bool = True):
    return {pos: ({REL_TYPE(rel): np.asarray(reads, data_type)
                   for rel, reads in muts[pos].items()}
                  if sanitize
                  else muts[pos])
            for pos in section.unmasked_int}


def simulate_muts(pmut: pd.DataFrame,
                  seg_end5s: np.ndarray,
                  seg_end3s: np.ndarray):
    """ Simulate mutation data.

    Parameters
    ----------
    pmut: pd.DataFrame
        Rate of each type of mutation at each position.
    seg_end5s:
        5' end coordinate of each segment.
    seg_end3s:
        3' end coordinate of each segment.
    """
    num_reads, _ = match_reads_segments(seg_end5s, seg_end3s)
    read_nums = np.arange(num_reads, dtype=fit_uint_type(num_reads))
    rels = np.asarray(pmut.columns.get_level_values(REL_NAME), dtype=REL_TYPE)
    if MATCH not in rels:
        raise ValueError(f"Relationships omit matches ({MATCH}): {rels}")
    if NOCOV in rels:
        raise ValueError(f"Relationships include no coverage ({NOCOV}): {rels}")
    muts = dict()
    for pos in index_to_pos(pmut.index):
        muts[pos] = dict()
        # Find the reads that cover this position.
        usable_reads = read_nums[np.any(np.logical_and(seg_end5s <= pos,
                                                       pos <= seg_end3s),
                                        axis=1)]
        if usable_reads.size > 0:
            # Choose a number of reads for each type of relationship.
            num_reads_pos_rels = pd.Series(rng.multinomial(usable_reads.size,
                                                           pmut.loc[pos])[0],
                                           index=rels).drop(MATCH)
            for rel, num_reads_pos_rel in num_reads_pos_rels.items():
                if num_reads_pos_rel > 0:
                    # Randomly select reads with this relationship.
                    reads_pos_rel = rng.choice(usable_reads,
                                               num_reads_pos_rel,
                                               replace=False,
                                               shuffle=False)
                    muts[pos][rel] = reads_pos_rel
                    # Prevent those reads from being chosen for another
                    # relationship.
                    usable_reads = np.setdiff1d(usable_reads,
                                                reads_pos_rel,
                                                assume_unique=True)
    return muts


@jit()
def _fill_matches(matrix: np.ndarray,
                  index5s: np.ndarray,
                  index3s: np.ndarray,
                  unmasked_read_indexes: np.ndarray):
    """ Fill all covered positions with matches. """
    for i, read_index in enumerate(unmasked_read_indexes):
        matrix[read_index, index5s[i]: index3s[i]] = MATCH


def calc_muts_matrix(section: Section,
                     read_nums: np.ndarray,
                     seg_end5s: np.ndarray,
                     seg_end3s: np.ndarray,
                     muts: dict[int, dict[int, np.ndarray]]):
    """ Matrix of relationships at each position in each read. """
    dims = find_dims([(NUM_READS,),
                      (NUM_READS, NUM_SEGMENTS),
                      (NUM_READS, NUM_SEGMENTS)],
                     [read_nums, seg_end5s, seg_end3s],
                     ["read_nums", "seg_end5s", "seg_end3s"])
    num_reads = dims[NUM_READS]
    num_segments = dims[NUM_SEGMENTS]
    section_unmasked = section.unmasked_int
    matrix = np.full((num_reads, section_unmasked.size), NOCOV)
    if matrix.size > 0:
        # Map each 5' and 3' end coordinate to its index in the unmasked
        # positions of the section.
        pos5_indexes = calc_inverse(section_unmasked,
                                    require=(section.end3 + 1),
                                    fill=True,
                                    fill_rev=True)
        pos3_indexes = calc_inverse(section_unmasked,
                                    require=section.end3,
                                    fill=True,
                                    fill_rev=False) + 1
        # Fill all covered positions with matches.
        read_indexes = np.arange(num_reads)
        for s in range(num_segments):
            end5s = seg_end5s[:, s]
            end3s = seg_end3s[:, s]
            if np.ma.is_masked(end5s) or np.ma.is_masked(end3s):
                unmasked_read_indexes = read_indexes[~(end5s.mask | end3s.mask)]
                end5s = end5s.data[unmasked_read_indexes]
                end3s = end3s.data[unmasked_read_indexes]
            else:
                unmasked_read_indexes = read_indexes
            if unmasked_read_indexes.size > 0:
                _fill_matches(matrix,
                              pos5_indexes[end5s],
                              pos3_indexes[end3s],
                              unmasked_read_indexes)
        # Overlay the mutation data.
        read_indexes = calc_inverse(read_nums)
        for pos in section_unmasked:
            if rels := muts.get(pos):
                column = matrix[:, pos5_indexes[pos]]
                for rel, reads in rels.items():
                    column[read_indexes[reads]] = rel
    return pd.DataFrame(matrix, read_nums, section.unmasked)


class MutsBatch(EndCoords, ReadBatch, ABC):
    """ Batch of mutational data. """

    def __init__(self, *,
                 section: Section,
                 sanitize: bool = True,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 **kwargs):
        super().__init__(section=section, sanitize=sanitize, **kwargs)
        # Validate and store the mutations.
        self._muts = sanitize_muts(muts, section, self.read_dtype, sanitize)

    @property
    def muts(self):
        """ Reads with each type of mutation at each position. """
        return self._muts

    @cached_property
    def pos_nums(self):
        """ Positions in use. """
        return np.fromiter(self.muts, self.pos_dtype)

    @property
    @abstractmethod
    def read_weights(self) -> pd.DataFrame | None:
        """ Weights for each read when computing counts. """

    @cached_property
    def read_end_counts(self):
        """ Counts of read end coordinates. """
        return count_end_coords(self.read_end5s,
                                self.read_end3s,
                                self.read_weights)


class SectionMutsBatch(MutsBatch, ABC):
    """ Batch of mutational data that knows its section. """

    def __init__(self, *, section: Section, **kwargs):
        self.section = section
        super().__init__(section=section, **kwargs)

    def iter_windows(self, size: int):
        """ Yield the positions in each window of the section. """
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
                             self.seg_end5s,
                             self.seg_end3s,
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

    @cached_property
    def matrix(self):
        """ Matrix of relationships at each position in each read. """
        return calc_muts_matrix(self.section,
                                self.read_nums,
                                self.seg_end5s,
                                self.seg_end3s,
                                self.muts)

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
                                   self.rels_per_read)

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
        if only_read_ends:
            end5s = self.read_end5s[:, np.newaxis]
            end3s = self.read_end3s[:, np.newaxis]
        else:
            end5s = self.seg_end5s
            end3s = self.seg_end3s
        for read_num, (end5, end3), muts, count in zip(self.read_nums,
                                                       zip(end5s,
                                                           end3s,
                                                           strict=True),
                                                       mut_pos,
                                                       mut_counts,
                                                       strict=True):
            yield read_num, ((tuple(end5), tuple(end3)), tuple(muts[:count]))


class PartialMutsBatch(MutsBatch, PartialReadBatch, ABC):
    pass

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
