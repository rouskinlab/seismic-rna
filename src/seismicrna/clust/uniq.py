from collections import defaultdict
from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from ..core.batch import RefseqMutsBatch, get_length
from ..core.rel import REL_TYPE, RelPattern
from ..core.seq import Section
from ..mask.data import MaskMerger

BIT_VECTOR_NAME = "Bit Vector"


class UniqReads(object):
    """ Collection of bit vectors of unique reads. """

    @classmethod
    def from_dataset(cls, dataset: MaskMerger):
        (muts_per_pos,
         batch_to_uniq,
         count_per_uniq) = get_uniq_reads(dataset.section.unmasked_int,
                                          dataset.pattern,
                                          dataset.iter_batches())
        return cls(dataset.sample,
                   dataset.ref,
                   dataset.section,
                   dataset.min_mut_gap,
                   muts_per_pos,
                   batch_to_uniq,
                   count_per_uniq)

    def __init__(self,
                 sample: str,
                 ref: str,
                 section: Section,
                 min_mut_gap: int,
                 muts_per_pos: list[np.ndarray],
                 batch_to_uniq: list[pd.Series],
                 counts_per_uniq: np.ndarray):
        self.sample = sample
        self.ref = ref
        self.section = section
        self.min_mut_gap = min_mut_gap
        if len(muts_per_pos) != self.num_pos:
            raise ValueError(f"Expected {self.num_pos} positions, "
                             f"but got {len(muts_per_pos)}")
        self.muts_per_pos = muts_per_pos
        self.batch_to_uniq = batch_to_uniq
        self.counts_per_uniq = counts_per_uniq

    @cached_property
    def pos_nums(self):
        return self.section.unmasked_int

    @cached_property
    def num_pos(self):
        """ Number of positions in each bit vector. """
        return get_length(self.pos_nums, "pos_nums")

    @cached_property
    def num_batches(self):
        """ Number of batches. """
        return len(self.batch_to_uniq)

    @cached_property
    def num_uniq(self):
        """ Number of unique reads. """
        return get_length(self.counts_per_uniq, "counts")

    @cached_property
    def num_reads(self):
        """ Number of total reads (including duplicates). """
        return np.sum(self.counts_per_uniq)

    def get_matrix(self):
        """ Full boolean matrix of the unique bit vectors. """
        # Initialize an all-False matrix with one row for each unique
        # bit vector and one column for each position.
        full = np.zeros((self.num_uniq, self.num_pos), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in enumerate(self.muts_per_pos):
            full[indexes, j] = True
        return full

    def get_uniq_names(self):
        """ Unique bit vectors as byte strings. """
        # Get the full boolean matrix of the unique bit vectors and cast
        # the data from boolean to unsigned 8-bit integer type.
        chars = self.get_matrix().astype(REL_TYPE, copy=False)
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
        if isinstance(other, UniqReads):
            return (self.sample == other.sample
                    and self.ref == other.ref
                    and self.section == other.section
                    and self.min_mut_gap == other.min_mut_gap
                    and np.array_equal(self.pos_nums,
                                       other.pos_nums)
                    and all(np.array_equal(m1, m2)
                            for m1, m2 in zip(self.muts_per_pos,
                                              other.muts_per_pos,
                                              strict=True))
                    and self.num_batches == other.num_batches
                    and all(b1.equals(b2)
                            for b1, b2 in zip(self.batch_to_uniq,
                                              other.batch_to_uniq,
                                              strict=True))
                    and np.array_equal(self.counts_per_uniq,
                                       other.counts_per_uniq))
        return NotImplemented


def uniq_reads_to_mutations(uniq_reads: Iterable[tuple],
                            pos_nums: Iterable[int]):
    """ Map each position to the numbers of the unique reads that are
    mutated at the position. """
    mutations = defaultdict(list)
    for uniq_read_num, read_muts_pos in enumerate(uniq_reads):
        for pos in read_muts_pos:
            mutations[pos].append(uniq_read_num)
    return [np.array(mutations[pos]) for pos in pos_nums]


def count_uniq_reads(uniq_read_nums: Iterable[list]):
    """ Count the occurrances of each unique value in the original. """
    return np.array(list(map(len, uniq_read_nums)))


def batch_to_uniq_read_num(read_nums_per_batch: list[np.ndarray],
                           uniq_read_nums: Iterable[list]):
    """ Map each read's number in its own batch to its unique number in
    the pool of all batches. """
    # For each batch, initialize a map from the read numbers of the .
    batch_to_uniq = [pd.Series(-1, index=read_nums)
                     for read_nums in read_nums_per_batch]
    # Fill in the map for each unique read.
    for uniq_read_num, batch_read_nums in enumerate(uniq_read_nums):
        for batch_num, read_num in batch_read_nums:
            batch_to_uniq[batch_num][read_num] = uniq_read_num
    for batch, read_nums in enumerate(batch_to_uniq):
        if read_nums.size > 0 > read_nums.min():
            raise ValueError(f"Got unmapped read numbers in batch {batch}: "
                             f"{read_nums.index[read_nums < 0]}")
    return batch_to_uniq


def get_uniq_reads(pos_nums: Iterable[int],
                   pattern: RelPattern,
                   batches: Iterable[RefseqMutsBatch]):
    uniq_reads = defaultdict(list)
    read_nums_per_batch = list()
    for batch_num, batch in enumerate(batches):
        if batch.batch != batch_num:
            raise ValueError(
                f"Batch {batch} is not in order (expected {batch_num})")
        read_nums_per_batch.append(batch.read_nums)
        for read_num, (e, m) in batch.iter_reads(pattern):
            uniq_reads[m].append((batch_num, read_num))
    muts_per_pos = uniq_reads_to_mutations(uniq_reads, pos_nums)
    batch_to_uniq = batch_to_uniq_read_num(read_nums_per_batch,
                                           uniq_reads.values())
    count_per_uniq = count_uniq_reads(uniq_reads.values())
    return muts_per_pos, batch_to_uniq, count_per_uniq

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
