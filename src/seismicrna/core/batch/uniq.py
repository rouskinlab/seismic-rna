from functools import cached_property
from typing import Iterable

import numpy as np

from .muts import MutsBatch
from .util import NO_READ, get_length, get_read_inverse
from ..rel import RelPattern
from ..seq import DNA, Section


class UniqReads(object):
    """ Collection of unique reads with only mutations. """

    def __init__(self,
                 section: Section,
                 uniq_reads_per_pos: dict[int, np.ndarray],
                 uniq_read_inverse: np.ndarray,
                 uniq_read_counts: np.ndarray):
        self.section = section
        self.muts = tuple(uniq_reads_per_pos.get(pos, np.array([], dtype=int))
                          for pos in section.unmasked_int)
        self.inverse = uniq_read_inverse
        self.counts = uniq_read_counts

    @cached_property
    def n_uniq(self):
        """ Number of unique reads. """
        return get_length(self.counts, "counts")

    @cached_property
    def n_pos(self):
        """ Number of positions in each bit vector. """
        return len(self.indexes)

    def get_full(self):
        """ Full boolean matrix of the unique bit vectors. """
        # Initialize an all-False matrix with one row for each unique
        # bit vector and one column for each position.
        full = np.zeros((self.n_uniq, self.n_pos), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in enumerate(self.indexes):
            full[indexes, j] = True
        return full

    def get_uniq_names(self):
        """ Return the unique bit vectors as byte strings. """
        # Get the full boolean matrix of the unique bit vectors and cast
        # the data from boolean to unsigned 8-bit integer type.
        chars = self.get_full().astype(REL_TYPE, copy=False)
        if chars.size > 0:
            # Add ord('0') to transform every 0 into b'0' and every 1
            # into # b'1', and convert each row (bit vector) into a
            # bytes object of b'0' and b'1' characters.
            names = np.apply_along_axis(np.ndarray.tobytes, 1, chars + ord('0'))
        else:
            # If there are no unique bit vectors, then apply_along_axis
            # will fail, so set names to an empty list.
            names = list()
        return pd.Index(names, name=BIT_VECTOR_NAME)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self.n_pos == other.n_pos
                and np.array_equal(self.counts, other.counts)
                and np.array_equal(self.inverse, other.inverse)
                and all(np.array_equal(self_idxs, other_idxs)
                        for self_idxs, other_idxs in zip(self.indexes,
                                                         other.indexes,
                                                         strict=True)))


def uniq(reads_per_pos: dict[int, np.ndarray],
         read_nums: np.ndarray,
         end_positions: np.ndarray):
    uniq_read_nums = list()
    uniq_read_counts = list()
    # Group the reads by their 5' and 3' ends.
    (uniq_ends,
     ends_indexes,
     ends_groups,
     ends_counts) = np.unique(end_positions,
                              return_index=True,
                              return_inverse=True,
                              return_counts=True,
                              axis=0)
    # All reads with unique 5' and 3' ends must be unique.
    has_uniq_ends = np.equal(ends_counts, 1)
    uniq_read_nums.append(read_nums[ends_indexes[has_uniq_ends]])
    uniq_read_counts.append(ends_counts[has_uniq_ends])
    # Check each group of reads that share 5' and 3' ends.
    for ends_group in np.arange(ends_indexes.size)[~has_uniq_ends]:
        # List the numbers of the reads in the group.
        group_read_nums = read_nums[ends_groups == ends_group]
        # Count the mutations in each read in the group.
        (mut_read_nums,
         mut_counts) = np.unique(np.hstack([np.intersect1d(pos_read_nums,
                                                           group_read_nums,
                                                           assume_unique=True)
                                            for pos_read_nums
                                            in reads_per_pos.values()]),
                                 return_counts=True)
        # Further group reads by their mutation counts.
        for mut_count in np.unique(mut_counts):
            # List the numbers of the reads with the mutation count.
            mut_count_read_nums = mut_read_nums[mut_counts == mut_count]
            mut_count_read_idxs = get_read_inverse(mut_count_read_nums)
            # Initialize a matrix of the mutated positions in each read.
            muts = np.zeros((mut_count_read_nums.size, mut_count), dtype=int)
            # Fill in the matrix one position at a time.
            for pos, pos_reads in reads_per_pos.items():
                # Reads at the position with this number of mutations.
                mut_count_pos_read_nums = np.intersect1d(pos_reads,
                                                         mut_count_read_nums,
                                                         assume_unique=True)
                # Find the rows corresponding to those reads.
                muts_rows = mut_count_read_idxs[mut_count_pos_read_nums]
                # Find the column of the first zero in each row.
                muts_cols = np.count_nonzero(muts[muts_rows], axis=1)
                # Fill in the position.
                muts[muts_rows, muts_cols] = pos
            # Find the reads with unique patterns of mutations.
            _, indexes, counts = np.unique(muts,
                                           return_index=True,
                                           return_counts=True,
                                           axis=0)
            # Record the original numbers of the unique reads.
            uniq_read_nums.append(mut_count_read_nums[indexes])
            # Record the number of times each unique read occurs.
            uniq_read_counts.append(counts)
        # Merge the unique reads over all numbers of mutations per read.
        uniq_read_inverse = get_read_inverse(np.hstack(uniq_read_nums))
        uniq_read_counts = np.hstack(uniq_read_counts)
        # Find the unique reads at each position.
        uniq_reads_per_pos = dict()
        for pos, pos_reads in reads_per_pos.items():
            uniq_pos_reads = uniq_read_inverse[pos_reads]
            uniq_reads_per_pos[pos] = uniq_pos_reads[uniq_pos_reads > NO_READ]
    return uniq_reads_per_pos, uniq_read_inverse, uniq_read_counts


def accum_uniq_reads_per_pos(refseq: DNA,
                             pattern: RelPattern,
                             batches: Iterable[MutsBatch]):
    for batch in batches:
        get_uniq_reads_per_pos(batch.reads_per_pos(refseq, pattern))
