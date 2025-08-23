from collections import defaultdict
from functools import partial

import numpy as np
import pandas as pd

from .ends import (END_COORDS,
                   match_reads_segments,
                   merge_read_ends,
                   sort_segment_ends)
from .index import count_base_types, iter_base_types
from ..array import find_dims, get_length
from ..logs import logger
from ..rel import MATCH, NOCOV, RelPattern
from ..seq import DNA, POS_NAME

POSITIONS = "positions"
READS = "reads"
CLUSTERS = "clusters"
SEGMENTS = "segments"
PRECISION = 3


def _count_end_coords(end5s: np.ndarray,
                      end3s: np.ndarray):
    """ Count each pair of 5' and 3' end coordinates. """
    uniq_ends, counts = np.unique(merge_read_ends(end5s, end3s),
                                  return_counts=True,
                                  axis=0)
    return pd.Series(counts,
                     pd.MultiIndex.from_arrays(uniq_ends.T, names=END_COORDS))


def _count_end_coords_weights(end5s: np.ndarray,
                              end3s: np.ndarray,
                              weights: pd.DataFrame):
    """ Count each pair of 5' and 3' end coordinates, with weights. """
    # Make a MultiIndex of all 5' and 3' coordinates.
    index = pd.MultiIndex.from_frame(pd.DataFrame(merge_read_ends(end5s,
                                                                  end3s),
                                                  columns=END_COORDS,
                                                  copy=False))
    # Convert the read weights into a DataFrame with that index.
    weights = pd.DataFrame(weights.values, index, weights.columns)
    # Sum the weights for each unique pair of 5'/3' coordinates.
    return weights.groupby(level=list(range(weights.index.nlevels))).sum()


def count_end_coords(end5s: np.ndarray,
                     end3s: np.ndarray,
                     weights: pd.DataFrame | None = None):
    """ Count each pair of 5' and 3' end coordinates. """
    if weights is not None:
        return _count_end_coords_weights(end5s, end3s, weights)
    return _count_end_coords(end5s, end3s)


def _calc_coverage(ends_sorted: np.ndarray,
                   is_contig_end5: np.ndarray,
                   is_contig_end3: np.ndarray,
                   read_weights: int | float | np.ndarray,
                   base_count: np.ndarray):
    """ Calculate the coverage per position and per read.

    Parameters
    ----------
    ends_sorted: np.ndarray
        Array (reads x ends) of the 5' and 3' ends of the segments in
        each read; must be sorted ascendingly over axis 1, with 5' and
        3' ends intermixed; 5' ends must be 0-indexed.
    is_contig_end5: np.ndarray
        Array (reads x ends) indicating whether each coordinate in
        `ends_sorted` is the 5' end of a contiguous segment.
    is_contig_end3: np.ndarray
        Array (reads x ends) indicating whether each coordinate in
        `ends_sorted` is the 3' end of a contiguous segment.
    read_weights: int | float | np.ndarray
        Array (reads x clusters) of the weight of each read per cluster.
    base_count: np.ndarray
        Array ((positions + 1) x bases) of the cumulative count of this
        kind of base up to each position.
    """
    num_reads, _ = ends_sorted.shape
    if isinstance(read_weights, np.ndarray):
        assert read_weights.ndim == 2
        _, num_clusters = read_weights.shape
        assert num_clusters >= 1
    else:
        assert isinstance(read_weights, (int, float))
        num_clusters = 1
    inc_pos, num_bases = base_count.shape
    num_pos = inc_pos - 1
    # Initialize coverage arrays.
    per_pos = np.zeros((num_pos, num_clusters))
    per_read = np.zeros((num_reads, num_bases), dtype=np.int64)
    # Get row indices and column indices for contiguous ends.
    row_indices_5, col_indices_5 = np.nonzero(is_contig_end5)
    row_indices_3, col_indices_3 = np.nonzero(is_contig_end3)
    assert np.array_equal(row_indices_5, row_indices_3)
    # Get the actual end coordinates.
    end5_coords = ends_sorted[row_indices_5, col_indices_5]
    end3_coords = ends_sorted[row_indices_3, col_indices_3]
    assert np.all(end5_coords <= end3_coords)
    assert np.all(end3_coords <= num_pos)
    # Calculate coverate per read for each type of base.
    np.add.at(per_read,
              row_indices_5,
              base_count[end3_coords] - base_count[end5_coords])
    # Process each cluster separately.
    for cluster in range(num_clusters):
        # Get weights for this cluster.
        if isinstance(read_weights, np.ndarray):
            weights_cluster = read_weights[row_indices_5, cluster]
        else:
            weights_cluster = read_weights
        # Create a difference array for this cluster.
        pos_diff = np.zeros(num_pos + 1)
        # Add weights at 5' positions
        np.add.at(pos_diff, end5_coords, weights_cluster)
        # Subtract weights at 3' positions.
        np.add.at(pos_diff, end3_coords, -weights_cluster)
        # Compute cumulative sum to get coverage.
        per_pos[:, cluster] = np.cumsum(pos_diff[:-1])
    return per_pos, per_read


def calc_coverage(pos_index: pd.Index,
                  read_nums: np.ndarray,
                  seg_end5s: np.ndarray,
                  seg_end3s: np.ndarray,
                  seg_ends_mask: np.ndarray | None,
                  read_weights: pd.DataFrame | None = None):
    """ Calculate the coverage per position and per read. """
    logger.routine("Began calculating coverage per position and per read")
    match_reads_segments(seg_end5s, seg_end3s, seg_ends_mask)
    # Find the positions in use.
    positions = pos_index.get_level_values(POS_NAME).values
    logger.detail(f"There are {positions.size} position(s) in use")
    if positions.size == 0:
        # If there are no positions in use, return empty arrays.
        cover_per_pos = (pd.DataFrame(0., pos_index, read_weights.columns)
                         if read_weights is not None
                         else pd.Series(0., pos_index))
        cover_per_read = pd.DataFrame.from_dict({base: pd.Series(0, read_nums)
                                                 for base in DNA.alph()})
        logger.routine("Ended calculating coverage per position and per read")
        return cover_per_pos, cover_per_read
    if positions.size > 1 and np.diff(positions).min() <= 0:
        raise ValueError("positions must increase monotonically, "
                         f"but got {positions}")
    min_pos = positions[0]
    max_pos = positions[-1]
    # Validate the dimensions.
    dim_names = [(POSITIONS,), (READS,), (READS, SEGMENTS), (READS, SEGMENTS)]
    arrays = [positions, read_nums, seg_end5s, seg_end3s]
    names = ["positions", "read_nums", "seg_end5s", "seg_end3s"]
    if read_weights is not None:
        if not isinstance(read_weights, pd.DataFrame):
            raise TypeError("If given, read_weights must be DataFrame, "
                            f"but got {type(read_weights).__name__}")
        read_weights = read_weights.loc[read_nums]
        dim_names.append((READS, CLUSTERS))
        arrays.append(read_weights.values)
        names.append("read_weights")
    find_dims(dim_names, arrays, names)
    # Clip the end coordinates to the minimum and maximum positions.
    # Sort the end coordinates and label the 3' end of each contiguous
    # segment.
    ends_sorted, is_contig_end5, is_contig_end3 = sort_segment_ends(
        seg_end5s.clip(min_pos, max_pos + 1),
        seg_end3s.clip(min_pos - 1, max_pos),
        seg_ends_mask
    )
    # Find the cumulative count of each base up to each position.
    bases = list()
    base_count = list()
    for base, base_pos in iter_base_types(pos_index):
        bases.append(base)
        is_base = np.zeros(max_pos + 1, dtype=bool)
        is_base[base_pos.get_level_values(POS_NAME)] = True
        base_count.append(np.cumsum(is_base))
    base_count = np.stack(base_count, axis=1)
    # Compute the coverage per position and per read.
    cover_per_pos, cover_per_read = _calc_coverage(ends_sorted,
                                                   is_contig_end5,
                                                   is_contig_end3,
                                                   (read_weights.values
                                                    if read_weights is not None
                                                    else 1),
                                                   base_count)
    cover_per_pos = cover_per_pos.round(PRECISION)
    # Reformat the coverage into pandas objects.
    cover_per_pos = cover_per_pos[positions - 1]
    if read_weights is not None:
        cover_per_pos = pd.DataFrame(cover_per_pos,
                                     index=pos_index,
                                     columns=read_weights.columns)
    else:
        cover_per_pos = pd.Series(cover_per_pos.reshape(positions.size),
                                  index=pos_index)
    cover_per_read = pd.DataFrame.from_dict(
        {base: pd.Series((cover_per_read[:, bases.index(base)]
                          if base in bases
                          else 0),
                         index=read_nums)
         for base in DNA.alph()},
    )
    logger.routine("Ended calculating coverage per position and per read")
    return cover_per_pos, cover_per_read


def calc_rels_per_pos(mutations: dict[int, dict[int, np.ndarray]],
                      num_reads: int | pd.Series,
                      cover_per_pos: pd.Series | pd.DataFrame,
                      read_indexes: np.ndarray | None = None,
                      read_weights: pd.DataFrame | None = None):
    """ For each relationship, the number of reads at each position. """
    slice_type = type(num_reads)
    array_type = type(cover_per_pos)
    pos_index = cover_per_pos.index
    if read_weights is not None:
        zero = 0.
        if not isinstance(read_weights, array_type):
            raise TypeError(f"Expected read_weights to be {array_type}, "
                            f"but got {type(read_weights)}")
        clusters = read_weights.columns
        if slice_type is not pd.Series:
            raise TypeError(f"Expected num_reads to be {pd.Series}, "
                            f"but got {slice_type}")
        slice_indexes = dict(index=clusters)
        if array_type is not pd.DataFrame:
            raise TypeError(f"Expected cover_per_pos to be {pd.DataFrame}, "
                            f"but got {array_type}")
        array_indexes = dict(index=pos_index, columns=clusters)
        if isinstance(read_indexes, np.ndarray):
            if read_indexes.ndim != 1:
                raise ValueError(f"Expected read_indexes to have 1 dimension, "
                                 f"but got {read_indexes.ndim}")
        else:
            raise TypeError(f"Expected read_indexes to be {np.ndarray}, "
                            f"but got {type(read_indexes)}")
        if not clusters.equals(num_reads.index):
            raise ValueError(f"Clusters differ between the number of reads "
                             f"({num_reads.index}) and the weights "
                             f"({clusters})")
        if not clusters.equals(cover_per_pos.columns):
            raise ValueError(f"Clusters differ between the coverage matrix "
                             f"({cover_per_pos.columns}) and the weights "
                             f"({clusters})")
    else:
        zero = 0
        if slice_type is not int:
            raise TypeError(f"Expected num_reads to be {int}, "
                            f"but got {slice_type}")
        slice_indexes = dict()
        if array_type is not pd.Series:
            raise TypeError(f"Expected cover_per_pos to be {pd.Series}, "
                            f"but got {array_type}")
        array_indexes = dict(index=pos_index)
    counts = defaultdict(partial(array_type, zero, **array_indexes))
    # Determine covered positions.
    if isinstance(cover_per_pos, pd.DataFrame):
        # A row is covered if its sum (across all clusters) > 0.
        covered = cover_per_pos.sum(axis=1) > 0
    else:
        covered = cover_per_pos > 0
    cov_indices = cover_per_pos[covered].index
    nocov_indices = cover_per_pos[~covered].index
    for pos_base in cov_indices:
        pos, base = pos_base
        num_reads_pos = slice_type(zero, **slice_indexes)
        for mut, reads in mutations.get(pos, dict()).items():
            if read_weights is not None:
                rows = read_indexes[reads]
                num_reads_pos_mut = read_weights.values[rows].sum(axis=0)
                num_reads_pos_mut = num_reads_pos_mut.round(PRECISION)
            else:
                num_reads_pos_mut = get_length(reads, "read numbers")
            num_reads_pos += num_reads_pos_mut
            counts[mut].loc[pos_base] = num_reads_pos_mut
        # The number of matches is the coverage minus the number of
        # reads with another kind of relationship that is not the
        # no-coverage relationship (no coverage is counted later).
        counts[MATCH].loc[pos_base] = (cover_per_pos.loc[pos_base]
                                       - num_reads_pos).round(PRECISION)
        if np.atleast_1d(counts[MATCH].loc[pos_base]).min(initial=0) < 0:
            raise ValueError("Number of matches must be ≥ 0, "
                             f"but got {counts[MATCH].loc[pos_base]} "
                             f"at position {pos}")
        # The number of non-covered positions is the number of reads
        # minus the number that cover the position.
        counts[NOCOV].loc[pos_base] = (
                num_reads - cover_per_pos.loc[pos_base]
        ).round(PRECISION)
        if np.atleast_1d(counts[NOCOV].loc[pos_base]).min(initial=0) < 0:
            raise ValueError("Number of non-covered positions must be ≥ 0, "
                             f"but got {counts[NOCOV].loc[pos_base]} "
                             f"at position {pos}")
    counts[MATCH].loc[nocov_indices] = 0
    if isinstance(cover_per_pos, pd.DataFrame):
        tot_reads = pd.DataFrame(
            np.repeat(num_reads.values[np.newaxis, :], len(nocov_indices), axis=0),
            index=nocov_indices,
            columns=num_reads.index
        )
    else:
        tot_reads = num_reads
    counts[NOCOV].loc[nocov_indices] = tot_reads
    return dict(counts)


def calc_rels_per_read(mutations: dict[int, dict[int, np.ndarray]],
                       pos_index: pd.Index,
                       cover_per_read: pd.DataFrame,
                       read_indexes: np.ndarray):
    """ For each relationship, the number of positions in each read. """
    counts = defaultdict(partial(pd.DataFrame,
                                 0,
                                 index=cover_per_read.index,
                                 columns=cover_per_read.columns))
    counts[NOCOV] = count_base_types(pos_index) - cover_per_read
    match = cover_per_read.copy()
    counts[MATCH] = match
    for pos, base in pos_index:
        for mut, reads in mutations[pos].items():
            if reads.size > 0:
                rows = read_indexes[reads]
                match[base].values[rows] -= 1
                counts[mut][base].values[rows] += 1
    errors = [(f"Number of matches for base {repr(base)} must be ≥ 0 "
               f"for every read, but got\n{values[values < 0]}")
              for base, values in match.items()
              if values.min() < 0]
    if errors:
        raise ValueError("\n".join(errors))
    return dict(counts)


def calc_reads_per_pos(pattern: RelPattern,
                       mutations: dict[int, dict[int, np.ndarray]],
                       pos_index: pd.Index):
    """ For each position, find all reads matching a pattern. """
    reads = dict()
    for pos, base in pos_index:
        pos_reads = [pos_mut_reads for mut, pos_mut_reads
                     in mutations.get(pos, dict()).items()
                     if all(pattern.fits(base, mut))]
        reads[pos] = (np.unique(np.hstack(pos_reads))
                      if pos_reads
                      else np.array([], int))
    return reads


def calc_covered_reads_per_pos(pos_index: pd.Index,
                               read_nums: np.ndarray,
                               seg_end5s: np.ndarray,
                               seg_end3s: np.ndarray,
                               seg_ends_mask: np.ndarray | None):
    """ For each position, find all reads covering it. """
    covering_reads = dict()
    for pos in pos_index.get_level_values(POS_NAME):
        covering_segs = np.logical_and(seg_end5s <= pos,
                                       seg_end3s >= pos)
        if seg_ends_mask is not None:
            covering_segs &= ~seg_ends_mask
        covering_reads[pos] = read_nums[covering_segs.any(axis=1)]
    return covering_reads


def calc_count_per_pos(pattern: RelPattern,
                       cover_per_pos: pd.Series | pd.DataFrame,
                       rels_per_pos: dict[int, pd.Series | pd.DataFrame]):
    """ Count the reads that fit a pattern at each position. """
    array_type = type(cover_per_pos)
    pos_index = cover_per_pos.index
    if array_type is pd.Series:
        zero = 0
        cov_pos_index = cover_per_pos[cover_per_pos > 0].index
        indexes = dict(index=pos_index)
    elif array_type is pd.DataFrame:
        zero = 0.
        cov_pos_index = cover_per_pos.index[cover_per_pos.sum(axis=1) > 0]
        indexes = dict(index=pos_index, columns=cover_per_pos.columns)
    else:
        raise TypeError(f"Expected cover_per_pos to be {pd.Series} or "
                        f"{pd.DataFrame}, but got {array_type}")
    info = array_type(zero, **indexes)
    fits = array_type(zero, **indexes)
    for base, index in iter_base_types(cov_pos_index):
        for rel, counts in rels_per_pos.items():
            is_info, is_fits = pattern.fits(base, rel)
            if is_info:
                pos_counts = counts.loc[index]
                info.loc[index] += pos_counts
                if is_fits:
                    fits.loc[index] += pos_counts
    return info, fits


def calc_count_per_read(pattern: RelPattern,
                        cover_per_read: pd.DataFrame,
                        rels_per_read: dict[int, pd.DataFrame]):
    """ Count the positions that fit a pattern in each read. """
    read_nums = cover_per_read.index
    info = pd.Series(0, index=read_nums)
    fits = pd.Series(0, index=read_nums)
    for rel, rel_counts in rels_per_read.items():
        for base, base_counts in rel_counts.items():
            is_info, is_fits = pattern.fits(str(base), rel)
            if is_info:
                read_counts = base_counts.values
                info += read_counts
                if is_fits:
                    fits += read_counts
    return info, fits
