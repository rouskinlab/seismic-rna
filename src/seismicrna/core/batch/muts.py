from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd
from numba import jit

from .confusion import calc_confusion_matrix
from .count import (
    calc_count_per_pos,
    calc_count_per_read,
    calc_coverage,
    calc_covered_reads_per_pos,
    calc_reads_per_pos,
    calc_rels_per_pos,
    calc_rels_per_read,
    count_end_coords,
)
from .ends import EndCoords, match_reads_segments
from .read import ReadBatch
from ..array import calc_inverse, find_dims
from ..header import REL_NAME, make_header
from ..rel import MATCH, NOCOV, REL_TYPE, RelPattern
from ..seq import DNA, Region, index_to_pos
from ..types import fit_uint_type

NUM_READS = "reads"
NUM_SEGMENTS = "segments"


def sanitize_muts(
    muts: dict[int, dict[int, list[int] | np.ndarray]],
    region: Region,
    data_type: type,
    sanitize: bool = True,
):
    """Keep only unmasked positions in the muts dictionary and convert
    the read lists to arrays.

    Parameters
    ----------
    muts: dict[int, dict[int, list[int] | np.ndarray]]
        Mapping from position to relationship code to read numbers.
    region: Region
        Region whose unmasked positions define which positions to keep.
    data_type: type
        NumPy dtype to use for the read number arrays.
    sanitize: bool = True
        If True, filter to unmasked positions and convert to arrays;
        if False, return the muts values as-is.

    Returns
    -------
    dict[int, dict[int, np.ndarray]]
        Sanitized mutation data.
    """
    return {
        int(pos): (
            {int(rel): np.asarray(reads, data_type) for rel, reads in muts[pos].items()}
            if sanitize
            else muts[pos]
        )
        for pos in region.unmasked_int
    }


def simulate_muts(
    pmut: pd.DataFrame, seg_end5s: np.ndarray, seg_end3s: np.ndarray, seed: int | None
):
    """Simulate mutation data.

    Parameters
    ----------
    pmut: pd.DataFrame
        Rate of each type of mutation at each position.
    seg_end5s: np.ndarray
        5' end coordinate of each segment.
    seg_end3s: np.ndarray
        3' end coordinate of each segment.
    seed: int | None
        Random seed for reproducibility; None for a random seed.

    Returns
    -------
    dict[int, dict[int, np.ndarray]]
        Mutation data: mapping from position to relationship code to
        read numbers.
    """
    rng = np.random.default_rng(seed)
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
        usable_reads = read_nums[
            np.any(np.logical_and(seg_end5s <= pos, pos <= seg_end3s), axis=1)
        ]
        if usable_reads.size > 0:
            # Choose a number of reads for each type of relationship.
            num_reads_pos_rels = pd.Series(
                rng.multinomial(usable_reads.size, pmut.loc[pos])[0], index=rels
            ).drop(MATCH)
            for rel, num_reads_pos_rel in num_reads_pos_rels.items():
                if num_reads_pos_rel > 0:
                    # Randomly select reads with this relationship.
                    reads_pos_rel = rng.choice(
                        usable_reads, num_reads_pos_rel, replace=False, shuffle=False
                    )
                    muts[pos][int(rel)] = reads_pos_rel
                    # Prevent those reads from being chosen for another
                    # relationship.
                    usable_reads = np.setdiff1d(
                        usable_reads, reads_pos_rel, assume_unique=True
                    )
    return muts


@jit()
def _fill_matches(
    matrix: np.ndarray,
    index5s: np.ndarray,
    index3s: np.ndarray,
    unmasked_read_indexes: np.ndarray,
):
    """Fill all covered positions with matches."""
    for i, read_index in enumerate(unmasked_read_indexes):
        matrix[read_index, index5s[i] : index3s[i]] = MATCH


def calc_muts_matrix(
    region: Region,
    read_nums: np.ndarray,
    seg_end5s: np.ndarray,
    seg_end3s: np.ndarray,
    seg_ends_mask: np.ndarray,
    muts: dict[int, dict[int, np.ndarray]],
):
    """Build a matrix of relationships at each position in each read.

    Parameters
    ----------
    region: Region
        Region providing unmasked positions.
    read_nums: np.ndarray
        Read numbers: 1D array (reads).
    seg_end5s: np.ndarray
        5' end coordinates of each segment: 2D array (reads x segments).
    seg_end3s: np.ndarray
        3' end coordinates of each segment: 2D array (reads x segments).
    seg_ends_mask: np.ndarray
        Boolean mask of segments to exclude: 2D array (reads x segments).
    muts: dict[int, dict[int, np.ndarray]]
        Mapping from position to relationship code to read numbers.

    Returns
    -------
    pd.DataFrame
        DataFrame of relationship codes, indexed by read number and
        columned by position index.
    """
    dims = find_dims(
        [(NUM_READS,), (NUM_READS, NUM_SEGMENTS), (NUM_READS, NUM_SEGMENTS)],
        [read_nums, seg_end5s, seg_end3s],
        ["read_nums", "seg_end5s", "seg_end3s"],
    )
    num_reads = dims[NUM_READS]
    num_segments = dims[NUM_SEGMENTS]
    region_unmasked = region.unmasked_int
    matrix = np.full((num_reads, region_unmasked.size), NOCOV)
    if matrix.size > 0:
        # Map each 5' and 3' end coordinate to its index in the unmasked
        # positions of the region.
        pos5_indexes = calc_inverse(
            region_unmasked, require=(region.end3 + 1), fill=True, fill_rev=True
        )
        pos3_indexes = (
            calc_inverse(
                region_unmasked, require=region.end3, fill=True, fill_rev=False
            )
            + 1
        )
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
                _fill_matches(
                    matrix,
                    pos5_indexes[end5s],
                    pos3_indexes[end3s],
                    unmasked_read_indexes,
                )
        # Overlay the mutation data.
        read_indexes = calc_inverse(read_nums)
        for pos in region_unmasked:
            if rels := muts.get(pos):
                column = matrix[:, pos5_indexes[pos]]
                for rel, reads in rels.items():
                    column[read_indexes[reads]] = rel
    return pd.DataFrame(matrix, read_nums, region.unmasked)


class MutsBatch(EndCoords, ReadBatch, ABC):
    """Batch of mutational data."""

    def __init__(
        self,
        *,
        region: Region,
        sanitize: bool = True,
        muts: dict[int, dict[int, list[int] | np.ndarray]],
        masked_read_nums: np.ndarray | list[int] | None = None,
        **kwargs,
    ):
        super().__init__(region=region, sanitize=sanitize, **kwargs)
        # Validate and store the mutations.
        self.muts = sanitize_muts(muts, region, self.read_dtype, sanitize)
        if masked_read_nums is not None:
            self.masked_read_nums = np.asarray(masked_read_nums, dtype=int)

    @cached_property
    def pos_nums(self):
        """Positions in use."""
        return np.fromiter(self.muts, self.pos_dtype)

    @property
    @abstractmethod
    def read_weights(self) -> pd.DataFrame | None:
        """Weights for each read when computing counts."""

    @cached_property
    def read_end_counts(self):
        """Counts of read end coordinates."""
        return count_end_coords(self.read_end5s, self.read_end3s, self.read_weights)


def _add_to_column(added: pd.Series | pd.DataFrame, frame: pd.DataFrame, column: str):
    """Add the values in `added` to the column `column` of `frame`."""
    frame_col = frame[column]
    if not frame_col.index.equals(added.index):
        raise ValueError(
            f"Got different indexes for frame {frame_col.index} "
            f"and added values {added.index}"
        )
    if isinstance(added, pd.DataFrame) and not frame_col.columns.equals(added.columns):
        raise ValueError(
            f"Got different columns for frame {frame_col.columns} "
            f"and added values {added.columns}"
        )
    frame[column] = (frame_col + added).values


class RegionMutsBatch(MutsBatch, ABC):
    """Batch of mutational data that knows its region."""

    def __init__(self, *, region: Region, **kwargs):
        self.region = region
        super().__init__(region=region, **kwargs)

    @cached_property
    def pos_index(self):
        """Index of unmasked positions and bases."""
        return self.region.unmasked

    @cached_property
    def _coverage(self):
        """Coverage per position and per read."""
        return calc_coverage(
            self.pos_index,
            self.read_nums,
            self.seg_end5s,
            self.seg_end3s,
            self.seg_ends_mask,
            self.read_weights,
        )

    @property
    def cover_per_pos(self):
        """Number of reads covering each position."""
        per_pos, _ = self._coverage
        return per_pos

    @property
    def cover_per_read(self):
        """Number of positions covered by each read."""
        _, per_read = self._coverage
        return per_read

    @cached_property
    def covered_reads_per_pos(self):
        """Reads covering each position."""
        return calc_covered_reads_per_pos(
            self.pos_index,
            self.read_nums,
            self.seg_end5s,
            self.seg_end3s,
            self.seg_ends_mask,
        )

    @cached_property
    def rels_per_pos(self):
        """For each relationship, the number of reads at each position
        with that relationship."""
        return calc_rels_per_pos(
            self.muts,
            self.num_reads,
            self.cover_per_pos,
            self.read_indexes,
            self.read_weights,
        )

    @cached_property
    def rels_per_read(self):
        """For each relationship, the number of positions in each read
        with that relationship."""
        return calc_rels_per_read(
            self.muts, self.pos_index, self.cover_per_read, self.read_indexes
        )

    @cached_property
    def matrix(self):
        """Matrix of relationships at each position in each read."""
        return calc_muts_matrix(
            self.region,
            self.read_nums,
            self.seg_end5s,
            self.seg_end3s,
            self.seg_ends_mask,
            self.muts,
        )

    def reads_per_pos(self, pattern: RelPattern):
        """For each position, find all reads matching a relationship
        pattern."""
        return calc_reads_per_pos(pattern, self.muts, self.pos_index)

    def count_per_pos(self, pattern: RelPattern):
        """Count the reads that fit a relationship pattern at each
        position in a region."""
        return calc_count_per_pos(pattern, self.cover_per_pos, self.rels_per_pos)

    def count_per_read(self, pattern: RelPattern):
        """Count the positions in a region that fit a relationship
        pattern in each read."""
        return calc_count_per_read(pattern, self.cover_per_read, self.rels_per_read)

    def count_all(
        self,
        patterns: dict[str, RelPattern],
        ks: Iterable[int] | None = None,
        *,
        count_ends: bool = True,
        count_pos: bool = True,
        count_read: bool = True,
    ):
        """Calculate all counts."""
        # Determine whether the data are clustered.
        header = make_header(rels=list(patterns), ks=ks)
        if header.get_is_clustered():
            zero = 0.0
            rel_header = header.get_rel_header()
        else:
            zero = 0
            rel_header = header
        # Initialize the counts to 0.
        count_per_pos = (
            pd.DataFrame(zero, self.region.unmasked, header.index)
            if count_pos
            else None
        )
        count_per_read = (
            pd.DataFrame(0, self.batch_read_index, rel_header.index)
            if count_read
            else None
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
        return (
            self.num_reads,
            (self.read_end_counts if count_ends else None),
            count_per_pos,
            count_per_read,
        )

    def calc_min_mut_dist(self, pattern: RelPattern):
        """For each read, calculate the smallest distance (i.e. the gap
        plus 1) between any two mutations."""
        # For each read, initialize the smallest distance between two
        # mutations to the length of the region, which is 1 more than
        # the maximum possible distance between two mutations.
        min_mut_dist = np.full(self.read_nums.size, self.region.length)
        # Keep track of the last mutated position in each read, with 0
        # meaning that no mutations have yet occurred.
        last_mut_pos = np.zeros(self.read_nums.size, dtype=int)
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(pattern)
        prev_pos = 0
        for pos, mut_reads in reads_per_pos.items():
            # This algorithm relies on valid positions being ≥ 1,
            # otherwise last_mut_pos will break.
            assert 1 <= pos <= self.region.end3, "Position out of bounds"
            # Positions must be in increasing order for this algorithm
            # to work correctly, but this should be guaranteed.
            assert pos > prev_pos, "Positions are not in increasing order"
            prev_pos = pos
            # Indexes of all reads with a mutation at this position.
            pos_mut_indexes = self.read_indexes[mut_reads]
            # Indexes of reads with a mutation at this position and at
            # any previous position.
            pos_multi_muts_indexes = pos_mut_indexes[last_mut_pos[pos_mut_indexes] > 0]
            # For reads with a mutation at this position and a previous
            # position, update the minimum distance between mutations.
            min_mut_dist[pos_multi_muts_indexes] = np.minimum(
                pos - last_mut_pos[pos_multi_muts_indexes],
                min_mut_dist[pos_multi_muts_indexes],
            )
            # Update the last mutated position.
            last_mut_pos[pos_mut_indexes] = pos
        # Finally, to make it easy to distinguish reads with fewer than
        # two mutations, set min_mut_dist for all such reads to 0.
        min_mut_dist[min_mut_dist == self.region.length] = 0
        return min_mut_dist

    def reads_noclose_muts(self, pattern: RelPattern, min_gap: int):
        """List the reads with no two mutations too close."""
        if min_gap < 0:
            raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
        if min_gap == 0:
            # No reads can have two mutations too close.
            return self.read_nums
        min_mut_dist = self.calc_min_mut_dist(pattern)
        return self.read_nums[np.logical_or(min_mut_dist == 0, min_mut_dist > min_gap)]

    def merge_close_muts(self, pattern: RelPattern, min_gap: int):
        """Return a new muts dictionary in which mutations closer than
        min_gap are merged into a single mutation, keeping only the
        3'-most mutation. This algorithm corrects for extra mutations
        occuring at non-modified positions shortly 5' of modifications,
        and therefore distinguishes probe-modified positions from all
        mutated positions."""
        mods = dict()
        # Keep track of the most recently encountered mutated position
        # that was identified as a probe-modified position, with 0
        # meaning that no modifications have yet occurred.
        last_mod_pos = np.zeros(self.read_nums.size, dtype=int)
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(pattern)
        # Iterate over positions from 3' to 5'.
        prev_pos = self.region.end3 + 1
        for pos in reversed(reads_per_pos.keys()):
            # This algorithm relies on valid positions being ≥ 1 and
            # ≤ self.region.end3.
            assert 1 <= pos <= self.region.end3, "Position out of bounds"
            # Positions must be in decreasing order for this algorithm
            # to work correctly, but this should be guaranteed.
            assert pos < prev_pos, "Positions are not in decreasing order"
            prev_pos = pos
            # Indexes of all reads with a mutation at this position.
            pos_mut_indexes = self.read_indexes[reads_per_pos[pos]]
            # Indexes of reads with a mutation at this position and no
            # true modification within min_gap positions on the 3' side.
            # These reads are considered to have a probe modification at
            # this position.
            last_mod_pos_of_pos_mut_indexes = last_mod_pos[pos_mut_indexes]
            pos_mod_indexes = pos_mut_indexes[
                (last_mod_pos_of_pos_mut_indexes == 0)
                | (last_mod_pos_of_pos_mut_indexes > (pos + min_gap))
            ]
            # Map the read indexes back to their read numbers.
            pos_mod_read_nums = self.read_nums[pos_mod_indexes]
            # Add the reads with probe modifications at this position to
            # the new mods dictionary. Artifact reads are those with
            # pattern-matching mutations at pos that are too close to a
            # true modification on the 3' side; all other reads (including
            # those with non-pattern relationships) are kept unchanged.
            artifact_reads = np.setdiff1d(
                reads_per_pos[pos], pos_mod_read_nums, assume_unique=True
            )
            mods[pos] = {
                rel: remaining
                for rel, reads in self.muts[pos].items()
                if (
                    remaining := np.setdiff1d(reads, artifact_reads, assume_unique=True)
                ).size
            }
            # Update the most recent probe-modified position.
            last_mod_pos[pos_mod_indexes] = pos
        # Restore the original order of positions.
        return {pos: mods[pos] for pos in self.muts.keys()}

    def inject_close_muts(
        self,
        pattern: RelPattern,
        injected_mut_probs: dict[int, float],
        seed: int | None,
    ):
        """Return a new muts dictionary with extra mutations injected
        at each offset 5' of existing mutations."""
        rng = np.random.default_rng(seed)
        # Make a deep copy of muts; for speed, make only shallow copies
        # of the NumPy arrays.
        new_muts = {pos: rels.copy() for pos, rels in self.muts.items()}
        # Validate the mutation probabilities.
        for offset, prob in injected_mut_probs.items():
            if offset < 1:
                raise ValueError(
                    f"All injected_mut_probs offsets must be ≥ 1, but got {offset}"
                )
            if not 0 <= prob <= 1:
                raise ValueError(
                    f"All injected_mut_probs values must be in [0, 1], but got {prob}"
                )
        if sum(injected_mut_probs.values()) == 0.0:
            # No mutations to inject.
            return new_muts
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(pattern)
        # Map each unmasked position to its reference base.
        pos_to_base = dict(self.pos_index)
        # Determine the valid mutations for each base.
        valid_muts = {
            base: np.array(
                [rel for rel in pattern.yes.mut_bits if pattern.yes.fits(base, rel)]
            )
            for base in DNA.alph()
        }
        # Iterate over positions from 3' to 5'.
        prev_pos = self.region.end3 + 1
        for pos in reversed(reads_per_pos.keys()):
            # This algorithm relies on valid positions being ≥ 1 and
            # ≤ self.region.end3.
            assert 1 <= pos <= self.region.end3, "Position out of bounds"
            # Positions must be in decreasing order for this algorithm
            # to work correctly, but this should be guaranteed.
            assert pos < prev_pos, "Positions are not in decreasing order"
            prev_pos = pos
            # Indexes of all reads with a mutation at this position.
            pos_mut_indexes = self.read_indexes[reads_per_pos[pos]]
            # 5'/3' ends of all reads with a mutation at this position.
            pos_mut_end5s = self.seg_end5s[pos_mut_indexes]
            pos_mut_end3s = self.seg_end3s[pos_mut_indexes]
            for offset, mut_prob in injected_mut_probs.items():
                new_pos = pos - offset
                if new_pos < 1 or new_pos not in new_muts:
                    # Skip out-of-bounds and masked positions.
                    continue
                # Whether each read has a segment that covers new_pos.
                # Note: self.seg_ends_mask does not need to be applied
                # because any masked segment has end5 > end3, and thus
                # the following inequality must be False automatically.
                covers_new_pos = np.any(
                    (pos_mut_end5s <= new_pos) & (new_pos <= pos_mut_end3s), axis=1
                )
                # Choose which reads to mutate at new_pos.
                mut_new_pos = np.logical_and(
                    covers_new_pos,
                    rng.binomial(n=1, p=mut_prob, size=covers_new_pos.size),
                )
                new_pos_indexes = pos_mut_indexes[mut_new_pos]
                if new_pos_indexes.size == 0:
                    continue
                # Map the indexes back to the read numbers.
                new_pos_reads = self.read_nums[new_pos_indexes]
                # Exclude reads that already have a mutation at new_pos
                # under any rel; each read must appear in at most one rel
                # per position.
                if new_muts[new_pos]:
                    already_mutated = np.unique(
                        np.concatenate(list(new_muts[new_pos].values()))
                    )
                    new_pos_reads = np.setdiff1d(
                        new_pos_reads, already_mutated, assume_unique=True
                    )
                if new_pos_reads.size == 0:
                    continue
                # Choose a valid mutation for each read in new_pos_reads.
                new_base = pos_to_base[new_pos]
                new_pos_valid_muts = valid_muts[new_base]
                if new_pos_valid_muts.size == 0:
                    raise ValueError(
                        "There are no relationships for "
                        f"base {repr(new_base)} at position "
                        f"{new_pos} that fit {pattern}"
                    )
                new_pos_muts = rng.choice(
                    new_pos_valid_muts, new_pos_reads.size, replace=True
                )
                unique_rels, inverse = np.unique(new_pos_muts, return_inverse=True)
                for i, rel in enumerate(unique_rels):
                    reads_for_rel = new_pos_reads[inverse == i]
                    try:
                        new_muts[new_pos][rel] = np.union1d(
                            new_muts[new_pos][rel], reads_for_rel
                        )
                    except KeyError:
                        new_muts[new_pos][rel] = reads_for_rel
        # Restore the original order of positions.
        return new_muts

    def calc_confusion_matrix(self, pattern: RelPattern, min_gap: int = 0):
        """Calculate the confusion matrix of mutations."""
        return calc_confusion_matrix(
            self.pos_index,
            self.covered_reads_per_pos,
            self.reads_per_pos(pattern),
            self.read_weights,
            min_gap=min_gap,
        )

    def iter_reads(
        self,
        pattern: RelPattern,
        only_read_ends: bool = False,
        require_contiguous: bool = False,
    ):
        """End coordinates and mutated positions in each read."""
        if require_contiguous and self.num_discontiguous:
            raise ValueError(
                "This function requires contiguous reads, but got "
                f"{self.num_discontiguous} discontiguous read(s)"
            )
        reads_per_pos = self.reads_per_pos(pattern)
        # Find the maximum number of mutations in any read.
        max_mut_count = max(
            sum(map(Counter, reads_per_pos.values()), start=Counter()).values(),
            default=0,
        )
        # Initialize a matrix of the positions mutated in each read.
        mut_pos = np.zeros((self.num_reads, max_mut_count), dtype=self.pos_dtype)
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
        for read_num, (end5, end3), muts, count in zip(
            self.read_nums,
            zip(end5s, end3s, strict=True),
            mut_pos,
            mut_counts,
            strict=True,
        ):
            yield read_num, ((tuple(end5), tuple(end3)), tuple(muts[:count]))
