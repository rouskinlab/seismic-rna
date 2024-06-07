from collections import defaultdict
from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .names import BIT_VECTOR_NAME
from ..core.array import get_length
from ..core.batch import EndCoords, SectionMutsBatch
from ..core.rel import RelPattern
from ..core.seq import Section
from ..mask.data import MaskMutsDataset


class UniqReads(EndCoords):
    """ Collection of bit vectors of unique reads. """

    @classmethod
    def from_dataset(cls, dataset: MaskMutsDataset, **kwargs):
        """ Get unique reads from a dataset. """
        ((seg_end5s, seg_end3s),
         muts_per_pos,
         batch_to_uniq,
         count_per_uniq) = get_uniq_reads(dataset.section.unmasked_int,
                                          dataset.pattern,
                                          dataset.iter_batches(),
                                          **kwargs)
        return cls(dataset.sample,
                   dataset.section,
                   dataset.min_mut_gap,
                   dataset.quick_unbias,
                   dataset.quick_unbias_thresh,
                   muts_per_pos,
                   batch_to_uniq,
                   count_per_uniq,
                   seg_end5s=seg_end5s,
                   seg_end3s=seg_end3s)

    @classmethod
    def from_dataset_contig(cls, dataset: MaskMutsDataset):
        """ Get unique reads from a dataset of contiguous reads. """
        return cls.from_dataset(dataset,
                                only_read_ends=True,
                                require_contiguous=True)

    def __init__(self,
                 sample: str,
                 section: Section,
                 min_mut_gap: int,
                 quick_unbias: bool,
                 quick_unbias_thresh: float,
                 muts_per_pos: list[np.ndarray],
                 batch_to_uniq: list[pd.Series],
                 counts_per_uniq: np.ndarray,
                 **kwargs):
        super().__init__(section=section, **kwargs)
        self.sample = sample
        self.section = section
        self.min_mut_gap = min_mut_gap
        self.quick_unbias = quick_unbias
        self.quick_unbias_thresh = quick_unbias_thresh
        if len(muts_per_pos) != (npos := get_length(self.section.unmasked_int,
                                                    "pos_nums")):
            raise ValueError(f"Expected {npos} positions, "
                             f"but got {len(muts_per_pos)}")
        self.muts_per_pos = muts_per_pos
        self.batch_to_uniq = batch_to_uniq
        self.counts_per_uniq = counts_per_uniq

    @property
    def ref(self):
        """ Reference name. """
        return self.section.ref

    @cached_property
    def seg_end5s_zero(self):
        """ 5' end of every segment (0-indexed in the section). """
        return self.seg_end5s - self.section.end5

    @cached_property
    def seg_end3s_zero(self):
        """ 3' end of every segment (0-indexed in the section). """
        return self.seg_end3s - self.section.end5

    @cached_property
    def read_end5s_zero(self):
        """ 5' end coordinates (0-indexed in the section). """
        return self.read_end5s - self.section.end5

    @cached_property
    def read_end3s_zero(self):
        """ 3' end coordinates (0-indexed in the section). """
        return self.read_end3s - self.section.end5

    @cached_property
    def num_batches(self):
        """ Number of batches. """
        return len(self.batch_to_uniq)

    @cached_property
    def num_uniq(self):
        """ Number of unique reads. """
        return get_length(self.counts_per_uniq, "counts")

    @cached_property
    def num_nonuniq(self) -> int:
        """ Number of total reads (including non-unique reads). """
        return self.counts_per_uniq.sum()

    def get_mut_matrix(self):
        """ Full boolean matrix of the mutations. """
        # Initialize an all-False matrix with one row for each unique
        # read and one column for each position.
        muts = np.zeros((self.num_uniq, self.section.length), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in zip(self.section.unmasked_zero,
                              self.muts_per_pos,
                              strict=True):
            muts[indexes, j] = True
        return muts

    def get_cov_matrix(self):
        """ Full boolean matrix of the covered positions. """
        # Initialize an all-False matrix with one row for each unique
        # read and one column for each position.
        covs = np.zeros((self.num_uniq, self.section.length), dtype=bool)
        # For each read (i), set the covered elements to True.
        for segment in range(self.num_segments):
            end5s = self.seg_end5s_zero[:, segment]
            end3s = self.seg_end3s_zero[:, segment] + 1
            for read, (end5, end3) in enumerate(zip(end5s, end3s)):
                covs[read, end5: end3] = True
        return covs

    def get_uniq_names(self):
        """ Unique bit vectors as byte strings. """
        # Get the full boolean matrices of the coverage and mutations of
        # the unique reads, cast the data from boolean to 8-bit integer
        # type, and merge them into one matrix via subtraction.
        chars = (self.get_mut_matrix().astype(np.int8, copy=False)
                 - 2 * (~self.get_cov_matrix()).astype(np.int8, copy=False))
        if chars.size > 0:
            # Add ord('0') to transform every 0 into b'0' and every 1
            # into # b'1', and convert each row (bit vector) into a
            # bytes object of b'0' and b'1' characters.
            names = np.char.decode(
                np.apply_along_axis(np.ndarray.tobytes, 1, chars + ord('0'))
            )
        else:
            # If there are no unique bit vectors, then apply_along_axis
            # will fail, so set names to an empty list.
            names = list()
        return pd.Index(names, name=BIT_VECTOR_NAME)

    def __eq__(self, other):
        if not isinstance(other, UniqReads):
            return NotImplemented
        return (self.sample == other.sample
                and self.section == other.section
                and self.min_mut_gap == other.min_mut_gap
                and self.num_batches == other.num_batches
                and np.array_equal(self.seg_end5s, other.seg_end5s)
                and np.array_equal(self.seg_end3s, other.seg_end3s)
                and all(np.array_equal(m1, m2)
                        for m1, m2 in zip(self.muts_per_pos,
                                          other.muts_per_pos,
                                          strict=True))
                and all(b1.equals(b2)
                        for b1, b2 in zip(self.batch_to_uniq,
                                          other.batch_to_uniq,
                                          strict=True))
                and np.array_equal(self.counts_per_uniq,
                                   other.counts_per_uniq))

    def __str__(self):
        return f"{type(self).__name__} of {self.sample} over {self.section}"


def _uniq_reads_to_ends_muts(uniq_reads: Iterable[tuple[tuple, tuple]],
                             pos_nums: Iterable[int]):
    """ Map each position to the numbers of the unique reads that are
    mutated at the position. """
    end5s = list()
    end3s = list()
    muts = defaultdict(list)
    for uniq_read_num, ((end5, end3), read_muts) in enumerate(uniq_reads):
        end5s.append(end5)
        end3s.append(end3)
        for pos in read_muts:
            muts[pos].append(uniq_read_num)
    end5s = np.array(end5s) if end5s else np.empty((0, 1), dtype=int)
    end3s = np.array(end3s) if end3s else np.empty((0, 1), dtype=int)
    muts = [np.array(muts[pos], dtype=int) for pos in pos_nums]
    return (end5s, end3s), muts


def _batch_to_uniq_read_num(read_nums_per_batch: list[np.ndarray],
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


def _count_uniq_reads(uniq_read_nums: Iterable[list]):
    """ Count the occurrances of each unique value in the original. """
    return np.fromiter(map(len, uniq_read_nums), dtype=int)


def get_uniq_reads(pos_nums: Iterable[int],
                   pattern: RelPattern,
                   batches: Iterable[SectionMutsBatch],
                   **kwargs):
    uniq_reads = defaultdict(list)
    read_nums_per_batch = list()
    for batch_num, batch in enumerate(batches):
        if batch.batch != batch_num:
            raise ValueError(
                f"Batch {batch} is not in order (expected {batch_num})"
            )
        # Record the number of reads in the batch.
        read_nums_per_batch.append(batch.read_nums)
        # Find the reads with unique end coordinates and mutations.
        for (read_num, read_data) in batch.iter_reads(pattern, **kwargs):
            # Key each read by its end coordinates and by the positions
            # at which it has mutations, so that redundant reads map to
            # the same key; record each read as a tuple of its batch and
            # its read number within the batch.
            uniq_reads[read_data].append((batch_num, read_num))
    # Pre-process the unique reads to extract necessary information.
    reads_ends, muts_per_pos = _uniq_reads_to_ends_muts(uniq_reads,
                                                        pos_nums)
    batch_to_uniq = _batch_to_uniq_read_num(read_nums_per_batch,
                                            uniq_reads.values())
    count_per_uniq = _count_uniq_reads(uniq_reads.values())
    return reads_ends, muts_per_pos, batch_to_uniq, count_per_uniq

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
