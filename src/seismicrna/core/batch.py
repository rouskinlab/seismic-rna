import re
from abc import ABC, abstractmethod
from collections import defaultdict
from functools import cache, cached_property
from logging import getLogger
from typing import Iterable

import numpy as np
import pandas as pd

from .pattern import RelPattern
from .rel import MATCH, NOCOV, REL_TYPE
from .sect import POS_NAME, Section
from .types import fit_uint_type, get_uint_type, UINT_NBYTES

logger = getLogger(__name__)

READ_NUM = "Read Number"
BATCH_NUM = "Batch Number"
INDEX_NAMES = BATCH_NUM, READ_NUM

# Whether numbers are 0-indexed or 1-indexed.
BATCH_INDEX = 0
PATTERN_INDEX = 0
POS_INDEX = 1


def _zero_indexed(indexes: np.ndarray, basis: int):
    """ Adjust indexes from basis-indexed to 0-indexed. """
    return indexes if basis == 0 else indexes - basis


def list_batch_nums(num_batches: int):
    """ List the batch numbers. """
    return list(range(BATCH_INDEX, BATCH_INDEX + num_batches))


def get_length(array: np.ndarray, what: str = "array") -> int:
    if array.ndim != 1:
        raise ValueError(f"{what} must have 1 dimension, but got {array.ndim}")
    length, = array.shape
    return length


def ensure_same_length(arr1: np.ndarray,
                       arr2: np.ndarray,
                       what1: str,
                       what2: str):
    if (len1 := get_length(arr1, what1)) != (len2 := get_length(arr2, what2)):
        raise ValueError(
            f"Lengths differ between {what1} ({len1}) and {what2} ({len2})")
    return len1


def ensure_order(array1: np.ndarray,
                 array2: np.ndarray,
                 what1: str = "array1",
                 what2: str = "array2",
                 gt_eq: bool = False):
    num_reads = ensure_same_length(array1, array2, what1, what2)
    ineq_func, ineq_sign = (np.less, '<') if gt_eq else (np.greater, '>')
    if np.any(is_err := ineq_func(array1, array2)):
        index = pd.Index(np.arange(num_reads)[is_err], name=READ_NUM)
        errors = pd.DataFrame.from_dict({what1: pd.Series(array1[is_err],
                                                          index=index),
                                         what2: pd.Series(array2[is_err],
                                                          index=index)})
        raise ValueError(f"Got {what1} {ineq_sign} {what2}:\n{errors}")
    return num_reads


def sanitize_values(values: Iterable[int],
                    lower_limit: int,
                    upper_limit: int,
                    whats: str = "values"):
    """ Validate and sort values, and return them as an array. """
    # Convert the values to an array and ensure it is one-dimensional.
    if not isinstance(values, (np.ndarray, list)):
        values = list(values)
    array = np.asarray(values)
    n_values = get_length(array, whats)
    if n_values == 0:
        # The array is empty.
        return np.array([], dtype=get_uint_type(min(UINT_NBYTES)))
    # Find and sort the unique values.
    array = np.unique(array)
    if array.size != n_values:
        raise ValueError(f"Got non-unique {whats}")
    # Validate the minimum and maximum values.
    min_value = array[0]
    max_value = array[-1]
    if min_value < lower_limit:
        raise ValueError(f"All {whats} must be ≥ {lower_limit}, but got "
                         f"{array[array < lower_limit]}")
    if max_value > upper_limit:
        raise ValueError(f"All {whats} must be ≤ {upper_limit}, but got "
                         f"{array[array > upper_limit]}")
    # Return the array as the smallest data type that will fit the data.
    return np.asarray(array, dtype=fit_uint_type(max_value))


def sanitize_pos(positions: Iterable[int], seq_length: int):
    """ Validate and sort positions, and return them as an array. """
    return sanitize_values(positions, POS_INDEX, seq_length, "positions")


def sanitize_ends(max_pos: int,
                  end5s: list[int] | np.ndarray,
                  mid5s: list[int] | np.ndarray,
                  mid3s: list[int] | np.ndarray,
                  end3s: list[int] | np.ndarray):
    pos_dtype = fit_uint_type(max_pos)
    end5s = np.asarray(end5s, pos_dtype)
    mid5s = np.asarray(mid5s, pos_dtype)
    mid3s = np.asarray(mid3s, pos_dtype)
    end3s = np.asarray(end3s, pos_dtype)
    # Verify 5' end ≥ min position
    ensure_order(end5s,
                 np.broadcast_to(POS_INDEX, end5s.shape),
                 "5' end positions",
                 f"minimum position ({POS_INDEX})",
                 gt_eq=True)
    # Verify 5' end ≤ 5' mid
    ensure_order(end5s, mid5s, "5' end positions", "5' middle positions")
    # Verify 5' end ≤ 3' mid
    ensure_order(end5s, mid3s, "5' end positions", "3' middle positions")
    # Verify 5' mid ≤ 3' end
    ensure_order(mid5s, end3s, "5' middle positions", "3' end positions")
    # Verify 3' mid ≤ 3' end
    ensure_order(mid3s, end3s, "3' middle positions", "3' end positions")
    # Verify 3' end ≤ max position
    ensure_order(end3s,
                 np.broadcast_to(max_pos, end3s.shape),
                 "3' end positions",
                 f"maximum position ({max_pos})")
    return end5s, mid5s, mid3s, end3s


class Batch(ABC):
    """ Batch of reads. """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z]*)batch", cls.__name__.lower()).groups()
        return btype

    def __init__(self, *, batch: int, **kwargs):
        super().__init__(**kwargs)
        self.batch = batch

    @cached_property
    @abstractmethod
    def read_nums(self) -> np.ndarray:
        """ Read numbers in use. """

    @property
    @abstractmethod
    def num_reads(self) -> int:
        """ Number of reads that are actually in use. """

    @property
    @abstractmethod
    def max_read(self) -> int:
        """ Maximum possible value for a read index. """

    @property
    def read_dtype(self):
        """ Data type for read numbers. """
        return fit_uint_type(self.max_read)

    @cached_property
    @abstractmethod
    def read_idx(self) -> np.ndarray:
        """ Map each read number to its index in self.read_nums. """


def _get_coverage_matrix_single(positions: np.ndarray,
                                end5s: np.ndarray,
                                end3s: np.ndarray):
    # Determine and validate the dimensions.
    npos = get_length(positions)
    nreads = get_length(end5s, "5' end positions")
    # Reshape the positions and 5'/3' ends to row and column vectors.
    pos_row = positions.reshape((1, npos))
    end5s_col = end5s.reshape((nreads, 1))
    end3s_col = end3s.reshape((nreads, 1))
    # Generate a boolean matrix where each element indicates whether
    # the read (row) covers the position (column).
    return np.logical_and(end5s_col <= pos_row, pos_row <= end3s_col)


def _get_coverage_matrix_pair(positions: np.ndarray,
                              end5s: np.ndarray,
                              mid5s: np.ndarray,
                              mid3s: np.ndarray,
                              end3s: np.ndarray):
    return np.logical_or(_get_coverage_matrix_single(positions, end5s, mid3s),
                         _get_coverage_matrix_single(positions, mid5s, end3s))


class RelsBatch(Batch, ABC):

    def __init__(self, *,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 seqlen: int,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 **kwargs):
        super().__init__(**kwargs)
        self.seqlen = seqlen
        (self.end5s,
         self.mid5s,
         self.mid3s,
         self.end3s) = sanitize_ends(self.max_pos, end5s, mid5s, mid3s, end3s)
        self.muts = {pos: {REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                           for rel, reads in muts[pos].items()}
                     for pos in sanitize_pos(muts, self.max_pos)}

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
    def coverage_matrix(self, section: Section):
        index = section.unmasked
        positions = index.get_level_values(POS_NAME).values
        if not np.all(isin := np.isin(positions,
                                      self.pos_nums,
                                      assume_unique=True)):
            raise ValueError(
                f"Positions in {section} are not in {self}: {positions[~isin]}")
        return pd.DataFrame(_get_coverage_matrix_pair(positions,
                                                      self.end5s,
                                                      self.mid5s,
                                                      self.mid3s,
                                                      self.end3s),
                            index=self.read_nums,
                            columns=index,
                            copy=False)

    @cache
    def cov_per_pos(self, section: Section):
        """ Number of reads covering each position. """
        cov_mat = self.coverage_matrix(section)
        return pd.Series(np.count_nonzero(cov_mat.values, axis=0),
                         index=cov_mat.columns,
                         dtype=self.read_dtype)

    @cache
    def cov_per_read(self, section: Section):
        """ Number of positions covered by each read. """
        cov_mat = self.coverage_matrix(section)
        return pd.DataFrame.from_dict(
            {base: pd.Series(np.count_nonzero(cov_mat.loc[:, index].values,
                                              axis=1),
                             index=cov_mat.index)
             for base, index in section.iter_bases()},
            orient="columns",
            dtype=self.pos_dtype
        )

    @cache
    def rels_per_pos(self, section: Section):
        """ For each type of relationship, the number of reads at each
        position with that relationship. """
        cov_per_pos = self.cov_per_pos(section)
        counts = defaultdict(lambda: pd.Series(self.read_dtype(0),
                                               cov_per_pos.index))
        for key, coverage in cov_per_pos.items():
            pos, base = key
            num_reads_pos = 0
            for mut, reads in self.muts.get(pos, dict()).items():
                num_reads_pos_mut = get_length(reads, "read numbers")
                num_reads_pos += num_reads_pos_mut
                counts[mut].loc[key] = num_reads_pos_mut
            # The number of matches is the coverage minus the number of
            # reads with another kind of relationship that is not the
            # no-coverage relationship (no coverage is counted later).
            counts[MATCH].loc[key] = coverage - num_reads_pos
            # The number of non-covered positions is the number of reads
            # minus the number that cover the position.
            counts[NOCOV].loc[key] = self.num_reads - coverage
        return dict(counts)

    @cache
    def rels_per_read(self, section: Section):
        """ For each type of relationship, the number of positions in
        each read with that relationship. """
        cov_per_read = self.cov_per_read(section)
        bases = list(cov_per_read.columns)
        counts = defaultdict(lambda: pd.DataFrame(self.pos_dtype(0),
                                                  index=cov_per_read.index,
                                                  columns=cov_per_read.columns))
        counts[NOCOV] = section.base_counts - cov_per_read
        counts[MATCH] = cov_per_read.copy()
        for pos, base in section.unmasked:
            column = bases.index(base)
            for mut, reads in self.muts.get(pos, dict()).items():
                rows = self.read_idx[reads]
                counts[mut].values[rows, column] += 1
                counts[MATCH].values[rows, column] -= 1
        return dict(counts)

    @cache
    def reads_per_pos(self, section: Section, pattern: RelPattern):
        """ For each position, find all reads matching a relationship
        pattern. """
        reads = dict()
        for pos, base in section.unmasked:
            pos_reads = [pos_mut_reads for mut, pos_mut_reads
                         in self.muts.get(pos, dict()).items()
                         if all(pattern.fits(base, mut))]
            reads[pos] = (np.hstack(pos_reads) if pos_reads
                          else np.array([], self.read_dtype))
        return reads

    def count_per_pos(self, section: Section, pattern: RelPattern):
        """ Count the reads that fit a relationship pattern at each
        position in a section. """
        cov_per_pos = self.cov_per_pos(section)
        info = pd.Series(self.read_dtype(0), cov_per_pos.index)
        fits = pd.Series(self.read_dtype(0), cov_per_pos.index)
        rels_per_pos = self.rels_per_pos(section)
        for base, index in section.iter_bases():
            for rel, counts in rels_per_pos.items():
                is_info, is_fits = pattern.fits(base, rel)
                if is_info:
                    pos_counts = counts.loc[index]
                    info.loc[index] += pos_counts
                    if is_fits:
                        fits.loc[index] += pos_counts
        return info, fits

    def count_per_read(self, section: Section, pattern: RelPattern):
        """ Count the positions in a section that fit a relationship
        pattern in each read. """
        cov_per_read = self.cov_per_read(section)
        info = pd.Series(self.pos_dtype(0), cov_per_read.index)
        fits = pd.Series(self.pos_dtype(0), cov_per_read.index)
        for rel, rel_counts in self.rels_per_read(section).items():
            for base, base_counts in rel_counts.items():
                is_info, is_fits = pattern.fits(base, rel)
                if is_info:
                    info += base_counts
                    if is_fits:
                        fits += base_counts
        return info, fits

    def nonprox_muts(self, section: Section, pattern: RelPattern, min_gap: int):
        """ List the reads with non-proximal mutations. """
        if min_gap < 0:
            raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
        if min_gap == 0:
            # No reads can have proximal mutations.
            return self.read_nums
        # For each position, list the reads with a mutation.
        reads_per_pos = self.reads_per_pos(section, pattern)
        # Track which reads have no proximal mutations.
        nonprox = np.ones(self.num_reads, dtype=bool)
        # Count the number of times each read occurs in a window.
        read_counts = np.zeros_like(nonprox, dtype=self.pos_dtype)
        # Track which positions are being counted.
        pos_counted = set()
        # Iterate over all windows in the section.
        for _, win_pos in section.windows(min_gap + 1):
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


def count_per_pos(section: Section,
                  patterns: RelPattern | list[RelPattern],
                  batches: Iterable[RelsBatch]):
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    if isinstance(patterns, RelPattern):
        # If just one pattern is given, then return two Series.
        info, fits = count_per_pos(section, [patterns], batches)
        return info[PATTERN_INDEX], fits[PATTERN_INDEX]
    # Initialize a count for each position and relationship pattern.
    index = section.unmasked
    columns = pd.RangeIndex(PATTERN_INDEX, PATTERN_INDEX + len(patterns))
    info = pd.DataFrame(0, index=index, columns=columns)
    fits = pd.DataFrame(0, index=index, columns=columns)
    # Accumulate the counts from the batches.
    for batch in batches:
        for column, pattern in enumerate(patterns, start=PATTERN_INDEX):
            i, f = batch.count_per_pos(section, pattern)
            info.loc[:, column] += i
            fits.loc[:, column] += f
    return info, fits


def count_per_read(section: Section,
                   patterns: RelPattern | list[RelPattern],
                   batches: Iterable[RelsBatch]):
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    if isinstance(patterns, RelPattern):
        # If just one pattern is given, then return two Series.
        info, fits = count_per_read(section, [patterns], batches)
        return info[PATTERN_INDEX], fits[PATTERN_INDEX]
    # Initialize a count for each relationship pattern.
    columns = pd.RangeIndex(PATTERN_INDEX, PATTERN_INDEX + len(patterns))
    info = list()
    fits = list()
    # Accumulate the counts from the batches.
    for batch in batches:
        # Create two new DataFrames to hold the counts for each pattern
        # from the batch.
        index = pd.MultiIndex.from_arrays([np.broadcast_to(batch.batch,
                                                           batch.num_reads),
                                           batch.read_nums],
                                          names=INDEX_NAMES)
        info.append(pd.DataFrame(batch.pos_dtype(0), index, columns))
        fits.append(pd.DataFrame(batch.pos_dtype(0), index, columns))
        # Count each pattern and add it to the counts for the batch.
        for column, pattern in enumerate(patterns, start=PATTERN_INDEX):
            i, f = batch.count_per_read(section, pattern)
            info[-1].loc[:, column] = i
            fits[-1].loc[:, column] = f
    if info and fits:
        # Concatenate the batches into one DataFrame.
        return pd.concat(info, axis=0), pd.concat(fits, axis=0)
    # There are no batches.
    return (pd.DataFrame(columns=columns, dtype=int),
            pd.DataFrame(columns=columns, dtype=int))

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
