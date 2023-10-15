from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
from functools import cache, cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .read import ReadBatch, MaskReadBatch
from .util import (POS_INDEX,
                   get_coverage_matrix,
                   get_length,
                   sanitize_ends,
                   sanitize_pos)
from ..rel import MATCH, NOCOV, REL_TYPE, RelPattern
from ..seq import BASE_NAME, DNA, seq_pos_to_index
from ..types import fit_uint_type


class MutsBatch(ReadBatch, ABC):

    def __init__(self, *,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 seqlen: int,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 sanitize: bool = True,
                 **kwargs):
        super().__init__(**kwargs)
        self.seqlen = seqlen
        (self.end5s,
         self.mid5s,
         self.mid3s,
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
    def max_pos(self):
        """ Maximum allowed position. """
        return self.seqlen + (POS_INDEX - 1)

    @property
    def pos_dtype(self):
        """ Data type for position numbers. """
        return fit_uint_type(self.max_pos)

    @cached_property
    @abstractmethod
    def pos_nums(self) -> np.ndarray:
        """ Positions in use. """

    @cache
    def pos_index(self, refseq: DNA):
        return seq_pos_to_index(refseq, self.pos_nums, POS_INDEX)

    @cached_property
    @abstractmethod
    def pattern(self) -> RelPattern | None:
        """ Relationship pattern. """

    @cache
    def fits(self, base: str, rel: int):
        return self.pattern is None or all(self.pattern.fits(base, rel))

    @cache
    def count_base_types(self, refseq: DNA):
        bases = self.pos_index(refseq).get_level_values(BASE_NAME)
        base_types, counts = np.unique(bases, return_counts=True)
        return pd.Series(counts, base_types)

    def iter_base_types(self, refseq: DNA):
        index = self.pos_index(refseq)
        bases = index.get_level_values(BASE_NAME)
        for base in refseq.alph():
            index_base = index[bases == base]
            if index_base.size > 0:
                yield base, index_base

    def iter_windows(self, size: int):
        """ Yield the positions in each window of size positions of the
        section. """
        if size < 1:
            raise ValueError(f"size must be ≥ 1, but got {size}")
        if self.num_pos > 0:
            # Create a Series with the positions as its index.
            positions = pd.Series(True, self.pos_nums)
            min_pos = positions.index[0]
            max_pos = positions.index[-1]
            # Define the 5' and 3' ends of the window.
            win5 = min_pos
            win3 = min(win5 + (size - 1), max_pos)
            # Yield the positions in each window.
            while win3 <= max_pos:
                yield (win5, win3), positions.loc[win5: win3].index.values
                win5 += 1
                win3 += 1

    @cache
    def coverage_matrix(self, refseq: DNA):
        return pd.DataFrame(np.logical_or(get_coverage_matrix(self.pos_nums,
                                                              self.end5s,
                                                              self.mid3s),
                                          get_coverage_matrix(self.pos_nums,
                                                              self.mid5s,
                                                              self.end3s)),
                            index=self.read_nums,
                            columns=self.pos_index(refseq),
                            copy=False)

    @cache
    def cover_per_pos(self, refseq: DNA):
        """ Number of reads covering each position. """
        cov_mat = self.coverage_matrix(refseq)
        return pd.Series(np.count_nonzero(cov_mat.values, axis=0),
                         index=cov_mat.columns,
                         dtype=self.read_dtype)

    @cache
    def cover_per_read(self, refseq: DNA):
        """ Number of positions covered by each read. """
        cov_mat = self.coverage_matrix(refseq)
        return pd.DataFrame.from_dict(
            {base: pd.Series(np.count_nonzero(cov_mat.loc[:, index].values,
                                              axis=1),
                             index=cov_mat.index)
             for base, index in self.iter_base_types(refseq)},
            orient="columns",
            dtype=self.pos_dtype
        )

    @cache
    def rels_per_pos(self, refseq: DNA):
        """ For each relationship that matches the pattern of the batch,
        the number of reads at each position with that relationship. """
        cov_per_pos = self.cover_per_pos(refseq)
        counts = defaultdict(lambda: pd.Series(self.read_dtype(0),
                                               cov_per_pos.index))
        for key, coverage in cov_per_pos.items():
            pos, base = key
            num_reads_pos = 0
            for mut, reads in self.muts.get(pos, dict()).items():
                num_reads_pos_mut = get_length(reads, "read numbers")
                num_reads_pos += num_reads_pos_mut
                if self.fits(base, mut):
                    counts[mut].loc[key] = num_reads_pos_mut
            # The number of matches is the coverage minus the number of
            # reads with another kind of relationship that is not the
            # no-coverage relationship (no coverage is counted later).
            if self.fits(base, MATCH):
                counts[MATCH].loc[key] = coverage - num_reads_pos
            # The number of non-covered positions is the number of reads
            # minus the number that cover the position.
            if self.fits(base, NOCOV):
                counts[NOCOV].loc[key] = self.num_reads - coverage
        return dict(counts)

    @cache
    def rels_per_read(self, refseq: DNA):
        """ For each relationship that matches the pattern of the batch,
        the number of positions in each read with that relationship. """
        cov_per_read = self.cover_per_read(refseq)
        bases = list(cov_per_read.columns)
        counts = defaultdict(lambda: pd.DataFrame(self.pos_dtype(0),
                                                  index=cov_per_read.index,
                                                  columns=cov_per_read.columns))
        if self.pattern is None:
            counts[NOCOV] = self.count_base_types(refseq) - cov_per_read
            counts[MATCH] = cov_per_read.copy()
            implicit_per_pos = False
        else:
            implicit_per_pos = True
        for pos, base in self.pos_index(refseq):
            column = bases.index(base)
            if implicit_per_pos:
                covered = self.coverage_matrix(refseq)[pos]
                if self.fits(base, NOCOV):
                    counts[NOCOV].values[:, column] += np.logical_not(covered)
                if self.fits(base, MATCH):
                    counts[MATCH].values[:, column] += covered
            for mut, reads in self.muts.get(pos, dict()).items():
                rows = self.read_idx[reads]
                if self.fits(base, MATCH):
                    counts[MATCH].values[rows, column] -= 1
                if self.fits(base, mut):
                    counts[mut].values[rows, column] += 1
        return dict(counts)

    @cache
    def reads_per_pos(self, refseq: DNA, pattern: RelPattern):
        """ For each position, find all reads matching a relationship
        pattern. """
        reads = dict()
        for pos, base in seq_pos_to_index(refseq, self.pos_nums, POS_INDEX):
            pos_reads = [pos_mut_reads for mut, pos_mut_reads
                         in self.muts.get(pos, dict()).items()
                         if all(pattern.fits(base, mut))]
            reads[pos] = (np.hstack(pos_reads) if pos_reads
                          else np.array([], self.read_dtype))
        return reads

    def count_per_pos(self, refseq: DNA, pattern: RelPattern):
        """ Count the reads that fit a relationship pattern at each
        position in a section. """
        cov_per_pos = self.cover_per_pos(refseq)
        info = pd.Series(self.read_dtype(0), cov_per_pos.index)
        fits = pd.Series(self.read_dtype(0), cov_per_pos.index)
        rels_per_pos = self.rels_per_pos(refseq)
        for base, index in self.iter_base_types(refseq):
            for rel, counts in rels_per_pos.items():
                is_info, is_fits = pattern.fits(base, rel)
                if is_info:
                    pos_counts = counts.loc[index]
                    info.loc[index] += pos_counts
                    if is_fits:
                        fits.loc[index] += pos_counts
        return info, fits

    def count_per_read(self, refseq: DNA, pattern: RelPattern):
        """ Count the positions in a section that fit a relationship
        pattern in each read. """
        cov_per_read = self.cover_per_read(refseq)
        info = pd.Series(self.pos_dtype(0), cov_per_read.index)
        fits = pd.Series(self.pos_dtype(0), cov_per_read.index)
        for rel, rel_counts in self.rels_per_read(refseq).items():
            for base, base_counts in rel_counts.items():
                is_info, is_fits = pattern.fits(base, rel)
                if is_info:
                    info += base_counts
                    if is_fits:
                        fits += base_counts
        return info, fits

    def nonprox_muts(self, refseq: DNA, pattern: RelPattern, min_gap: int):
        """ List the reads with non-proximal mutations. """
        if min_gap < 0:
            raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
        if min_gap == 0:
            # No reads can have proximal mutations.
            return self.read_nums
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(refseq, pattern)
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
                read_counts[self.read_idx[reads_per_pos[pos]]] += 1
                pos_counted.add(pos)
            for pos in pos_counted - win_pos_set:
                # These positions have exited the window: drop them.
                read_counts[self.read_idx[reads_per_pos[pos]]] -= 1
                pos_counted.remove(pos)
            if len(pos_counted) > 1:
                # Check for reads counted more than once in the window,
                # which means that they have proximal mutations.
                nonprox &= read_counts <= 1
        return self.read_nums[nonprox]

    def mask(self,
             pattern: RelPattern | None = None,
             positions: Iterable[int] | None = None,
             reads: Iterable[int] | None = None):
        # Clean and validate the selection.
        if positions is not None:
            positions = sanitize_pos(positions, self.max_pos)
        # Select mutations at each position.
        muts = dict()
        for pos in positions if positions is not None else self.pos_nums:
            muts[pos] = dict()
            # Select reads with each type of mutation at this position.
            for mut, pos_mut_reads in self.muts.get(pos, dict()).items():
                muts[pos][mut] = (np.intersect1d(pos_mut_reads,
                                                 reads,
                                                 assume_unique=True)
                                  if reads is not None
                                  else pos_mut_reads)
        if reads is not None:
            read_nums = np.asarray(reads, dtype=self.read_dtype)
            read_indexes = self.read_idx[read_nums]
            end5s = self.end5s[read_indexes]
            mid5s = self.mid5s[read_indexes]
            mid3s = self.mid3s[read_indexes]
            end3s = self.end3s[read_indexes]
        else:
            read_nums = self.read_nums
            end5s = self.end5s
            mid5s = self.mid5s
            mid3s = self.mid3s
            end3s = self.end3s
        return MaskMutsBatch(batch=self.batch,
                             muts=muts,
                             seqlen=self.seqlen,
                             end5s=end5s,
                             mid5s=mid5s,
                             mid3s=mid3s,
                             end3s=end3s,
                             read_nums=read_nums,
                             max_read=self.max_read,
                             pattern=pattern,
                             sanitize=False)


class MaskMutsBatch(MaskReadBatch, MutsBatch, ABC):

    def __init__(self, *, pattern: RelPattern | None, **kwargs):
        super().__init__(**kwargs)
        self._pattern = pattern

    @cached_property
    def pos_nums(self):
        return np.array(list(self.muts))

    @property
    def pattern(self):
        return self._pattern


def accumulate(positions: np.ndarray,
               refseq: DNA,
               patterns: dict[str, RelPattern],
               batches: Iterable[MutsBatch],
               per_pos: bool = True,
               per_read: bool = True):
    columns = pd.Index(list(patterns))
    # Initialize the counts per position.
    if per_pos:
        index_per_pos = seq_pos_to_index(refseq, positions, POS_INDEX)
        info_per_pos = pd.DataFrame(0, index_per_pos, columns)
        fits_per_pos = pd.DataFrame(0, index_per_pos, columns)
    else:
        info_per_pos = None
        fits_per_pos = None
    # Initialize the counts per read.
    if per_read:
        info_per_read_per_batch = list()
        fits_per_read_per_batch = list()
    else:
        info_per_read_per_batch = None
        fits_per_read_per_batch = None
    # Accumulate the counts from the batches.
    for batch in batches:
        # Confirm that the positions of every batch match those given.
        if not np.array_equal(batch.pos_nums, positions):
            raise ValueError(f"Positions of {batch} ({batch.pos_nums}) do not "
                             f"match the given positions ({positions})")
        if info_per_read_per_batch is not None:
            # Make two DataFrames for the per-read counts of this batch.
            info_per_read_per_batch.append(pd.DataFrame(batch.pos_dtype(0),
                                                        batch.multiindex,
                                                        columns))
            fits_per_read_per_batch.append(pd.DataFrame(batch.pos_dtype(0),
                                                        batch.multiindex,
                                                        columns))
        # Count the positions and/or reads matching each pattern.
        for column, pattern in patterns.items():
            if info_per_pos is not None:
                # Count the matching reads per position.
                ipp, fpp = batch.count_per_pos(refseq, pattern)
                info_per_pos.loc[:, column] += ipp
                fits_per_pos.loc[:, column] += fpp
            if info_per_read_per_batch is not None:
                # Count the matching positions per read.
                ipr, fpr = batch.count_per_read(refseq, pattern)
                info_per_read_per_batch[-1].loc[:, column] = ipr
                fits_per_read_per_batch[-1].loc[:, column] = fpr
    if info_per_read_per_batch is not None:
        # Produce two DataFrames of per-read counts.
        if info_per_read_per_batch:
            # Concatenate the per-read counts for the batches.
            info_per_read = pd.concat(info_per_read_per_batch, axis=0)
            fits_per_read = pd.concat(fits_per_read_per_batch, axis=0)
        else:
            # There are no batches.
            info_per_read = pd.DataFrame(columns=columns, dtype=int)
            fits_per_read = pd.DataFrame(columns=columns, dtype=int)
    else:
        # Per-read information is not being counted.
        info_per_read = None
        fits_per_read = None
    return (info_per_pos, fits_per_pos), (info_per_read, fits_per_read)


def count_per_pos(positions: np.ndarray,
                  refseq: DNA,
                  patterns: dict[str, RelPattern],
                  batches: Iterable[MutsBatch]
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    per_pos, per_read = accumulate(positions,
                                   refseq,
                                   patterns,
                                   batches,
                                   per_read=False)
    return per_pos


def count_per_read(positions: np.ndarray,
                   refseq: DNA,
                   patterns: dict[str, RelPattern],
                   batches: Iterable[MutsBatch]
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    per_pos, per_read = accumulate(positions,
                                   refseq,
                                   patterns,
                                   batches,
                                   per_pos=False)
    return per_read

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
