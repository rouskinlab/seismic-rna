from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd
from numba import jit

from .confusion import calc_confusion_matrix
from .count import (calc_count_per_pos,
                    calc_count_per_read,
                    calc_coverage,
                    calc_covered_reads_per_pos,
                    calc_reads_per_pos,
                    calc_rels_per_pos,
                    calc_rels_per_read,
                    count_end_coords)
from .ends import EndCoords, match_reads_segments
from .read import ReadBatch
from ..array import calc_inverse, find_dims
from ..header import REL_NAME, make_header
from ..rel import MATCH, NOCOV, REL_TYPE, RelPattern
from ..seq import Region, index_to_pos
from ..types import fit_uint_type

rng = np.random.default_rng()

NUM_READS = "reads"
NUM_SEGMENTS = "segments"


def sanitize_muts(muts: dict[int, dict[int, list[int] | np.ndarray]],
                  region: Region,
                  data_type: type,
                  sanitize: bool = True):
    return {int(pos): ({int(rel): np.asarray(reads, data_type)
                        for rel, reads in muts[pos].items()}
                       if sanitize
                       else muts[pos])
            for pos in region.unmasked_int}


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
    num_reads, _ = match_reads_segments(seg_end5s, seg_end3s, None)
    read_nums = np.arange(num_reads, dtype=fit_uint_type(num_reads))
    rels = np.asarray(pmut.columns.get_level_values(REL_NAME), dtype=REL_TYPE)
    if MATCH not in rels:
        raise ValueError(f"Relationships omit matches ({MATCH}): {rels}")
    if NOCOV in rels:
        raise ValueError(f"Relationships include no coverage ({NOCOV}): {rels}")
    muts = dict()
    for pos in index_to_pos(pmut.index):
        muts[int(pos)] = dict()
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
                    muts[pos][int(rel)] = reads_pos_rel
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


def calc_muts_matrix(region: Region,
                     read_nums: np.ndarray,
                     seg_end5s: np.ndarray,
                     seg_end3s: np.ndarray,
                     seg_ends_mask: np.ndarray,
                     muts: dict[int, dict[int, np.ndarray]]):
    """ Matrix of relationships at each position in each read. """
    dims = find_dims([(NUM_READS,),
                      (NUM_READS, NUM_SEGMENTS),
                      (NUM_READS, NUM_SEGMENTS)],
                     [read_nums, seg_end5s, seg_end3s],
                     ["read_nums", "seg_end5s", "seg_end3s"])
    num_reads = dims[NUM_READS]
    num_segments = dims[NUM_SEGMENTS]
    region_unmasked = region.unmasked_int
    matrix = np.full((num_reads, region_unmasked.size), NOCOV)
    if matrix.size > 0:
        # Map each 5' and 3' end coordinate to its index in the unmasked
        # positions of the region.
        pos5_indexes = calc_inverse(region_unmasked,
                                    require=(region.end3 + 1),
                                    fill=True,
                                    fill_rev=True)
        pos3_indexes = calc_inverse(region_unmasked,
                                    require=region.end3,
                                    fill=True,
                                    fill_rev=False) + 1
        # Fill all covered positions with matches.
        read_indexes = np.arange(num_reads)
        for s in range(num_segments):
            if seg_ends_mask is not None:
                mask = seg_ends_mask[:, s]
                unmasked_read_indexes = read_indexes[~mask]
            else:
                unmasked_read_indexes = read_indexes
            end5s = seg_end5s[unmasked_read_indexes, s]
            end3s = seg_end3s[unmasked_read_indexes, s]
            if unmasked_read_indexes.size > 0:
                _fill_matches(matrix,
                              pos5_indexes[end5s],
                              pos3_indexes[end3s],
                              unmasked_read_indexes)
        # Overlay the mutation data.
        read_indexes = calc_inverse(read_nums)
        for pos in region_unmasked:
            if rels := muts.get(pos):
                column = matrix[:, pos5_indexes[pos]]
                for rel, reads in rels.items():
                    column[read_indexes[reads]] = rel
    return pd.DataFrame(matrix, read_nums, region.unmasked)


class MutsBatch(EndCoords, ReadBatch, ABC):
    """ Batch of mutational data. """

    def __init__(self, *,
                 region: Region,
                 sanitize: bool = True,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 masked_read_nums: np.ndarray | list[int] | None = None,
                 **kwargs):
        super().__init__(region=region, sanitize=sanitize, **kwargs)
        # Validate and store the mutations.
        self.muts = sanitize_muts(muts, region, self.read_dtype, sanitize)
        if masked_read_nums is not None:
            self.masked_read_nums = np.asarray(masked_read_nums, dtype=int)

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


def _add_to_column(added: pd.Series | pd.DataFrame,
                   frame: pd.DataFrame,
                   column: str):
    """ Add the values in `added` to the column `column` of `frame`. """
    frame_col = frame[column]
    if not frame_col.index.equals(added.index):
        raise ValueError(f"Got different indexes for frame {frame_col.index} "
                         f"and added values {added.index}")
    if (isinstance(added, pd.DataFrame)
            and not frame_col.columns.equals(added.columns)):
        raise ValueError(f"Got different columns for frame {frame_col.columns} "
                         f"and added values {added.columns}")
    frame[column] = (frame_col + added).values


class RegionMutsBatch(MutsBatch, ABC):
    """ Batch of mutational data that knows its region. """

    def __init__(self, *, region: Region, **kwargs):
        self.region = region
        super().__init__(region=region, **kwargs)

    @cached_property
    def pos_index(self):
        """ Index of unmasked positions and bases. """
        return self.region.unmasked

    @cached_property
    def _coverage(self):
        """ Coverage per position and per read. """
        return calc_coverage(self.pos_index,
                             self.read_nums,
                             self.seg_end5s,
                             self.seg_end3s,
                             self.seg_ends_mask,
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
    def covered_reads_per_pos(self):
        """ Reads covering each position. """
        return calc_covered_reads_per_pos(self.pos_index,
                                          self.read_nums,
                                          self.seg_end5s,
                                          self.seg_end3s,
                                          self.seg_ends_mask)

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
        return calc_muts_matrix(self.region,
                                self.read_nums,
                                self.seg_end5s,
                                self.seg_end3s,
                                self.seg_ends_mask,
                                self.muts)

    def reads_per_pos(self, pattern: RelPattern):
        """ For each position, find all reads matching a relationship
        pattern. """
        return calc_reads_per_pos(pattern, self.muts, self.pos_index)

    def count_per_pos(self, pattern: RelPattern):
        """ Count the reads that fit a relationship pattern at each
        position in a region. """
        return calc_count_per_pos(pattern,
                                  self.cover_per_pos,
                                  self.rels_per_pos)

    def count_per_read(self, pattern: RelPattern):
        """ Count the positions in a region that fit a relationship
        pattern in each read. """
        return calc_count_per_read(pattern,
                                   self.cover_per_read,
                                   self.rels_per_read)

    def count_all(self,
                  patterns: dict[str, RelPattern],
                  ks: Iterable[int] | None = None, *,
                  count_ends: bool = True,
                  count_pos: bool = True,
                  count_read: bool = True):
        """ Calculate all counts. """
        # Determine whether the data are clustered.
        header = make_header(rels=list(patterns), ks=ks)
        if header.get_is_clustered():
            zero = 0.
            rel_header = header.get_rel_header()
        else:
            zero = 0
            rel_header = header
        # Initialize the counts to 0.
        count_per_pos = (
            pd.DataFrame(zero, self.region.unmasked, header.index)
            if count_pos else None
        )
        count_per_read = (
            pd.DataFrame(0, self.batch_read_index, rel_header.index)
            if count_read else None
        )
        for column, pattern in patterns.items():
            if count_per_pos is not None:
                # Count the matching reads per position.
                _, count_per_pos_pattern = self.count_per_pos(pattern)
                _add_to_column(count_per_pos_pattern, count_per_pos, column)
            if count_per_read is not None:
                # Count the matching positions per read.
                _, count_per_read_pattern = self.count_per_read(pattern)
                count_per_read.loc[:, column] = count_per_read_pattern.values
        return (self.num_reads,
                (self.read_end_counts if count_ends else None),
                count_per_pos,
                count_per_read)

    def calc_min_mut_dist(self, pattern: RelPattern):
        """ For each read, calculate the smallest distance (i.e. the gap
        plus 1) between any two mutations. """
        # For each read, initialize the smallest distance between two
        # mutations to the length of the region, which is 1 more than
        # the maximum possible distance between two mutations.
        min_mut_dist = np.full(self.read_nums.size, self.region.length)
        # Keep track of the last mutated position in each read, with 0
        # meaning that 0 mutations have yet occurred.
        last_mut_pos = np.zeros_like(min_mut_dist)
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(pattern)
        for pos, reads in reads_per_pos.items():
            if pos < 1:
                # This algorithm relies on valid positions being ≥ 1,
                # otherwise last_mut_pos will break.
                raise ValueError(f"Position must be ≥ 1, but got {pos}")
            # Indexes of all reads with a mutation at this position.
            pos_indexes = self.read_indexes[reads]
            # Indexes of reads with a mutation at this position and at
            # any previous position.
            pos_mut_indexes = pos_indexes[last_mut_pos[pos_indexes] > 0]
            # For reads with a mutation at this position and a previous
            # position, update the minimum distance between mutations.
            min_mut_dist[pos_mut_indexes] = np.minimum(
                pos - last_mut_pos[pos_mut_indexes],
                min_mut_dist[pos_mut_indexes]
            )
            # Update the last mutated position.
            last_mut_pos[pos_indexes] = pos
        # Finally, to make it easy to distinguish reads with fewer than
        # two mutations, set min_mut_dist for all such reads to 0.
        min_mut_dist[min_mut_dist == self.region.length] = 0
        return min_mut_dist

    def reads_noclose_muts(self, pattern: RelPattern, min_gap: int):
        """ List the reads with no two mutations too close. """
        if min_gap < 0:
            raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
        if min_gap == 0:
            # No reads can have two mutations too close.
            return self.read_nums
        min_mut_dist = self.calc_min_mut_dist(pattern)
        return self.read_nums[np.logical_or(min_mut_dist == 0,
                                            min_mut_dist > min_gap)]

    def calc_confusion_matrix(self, pattern: RelPattern, min_gap: int = 0):
        """ Calculate the confusion matrix of mutations. """
        return calc_confusion_matrix(self.pos_index,
                                     self.covered_reads_per_pos,
                                     self.reads_per_pos(pattern),
                                     self.read_weights,
                                     min_gap=min_gap)

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
