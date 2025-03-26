from functools import cached_property

import numpy as np

from ..array import ensure_order, ensure_same_length, find_dims
from ..logs import logger
from ..seq import Region
from ..types import fit_uint_type
from ..validate import require_isinstance, require_equal, require_atleast

rng = np.random.default_rng()

# Indexes of read end coordinates.
END5_COORD = "5' End"
END3_COORD = "3' End"
END_COORDS = [END5_COORD, END3_COORD]

NUM_UNIQ = "uniq"
NUM_SEGS = "segments"


class BadSegmentEndsError(ValueError):
    """ Segment 5' and 3' ends are not valid. """


def count_reads_segments(seg_ends: np.ndarray,
                         what: str = "seg_ends") -> tuple[int, int]:
    require_isinstance(what, seg_ends, np.ndarray)
    require_equal(f"{what}.ndim", seg_ends.ndim, 2,
                  error_type=BadSegmentEndsError)
    num_reads, num_segs = seg_ends.shape
    return num_reads, num_segs


def match_reads_segments(seg_end5s: np.ndarray,
                         seg_end3s: np.ndarray,
                         seg_ends_mask: np.ndarray | None):
    """ Number of segments for the given end coordinates. """
    num_reads_5, num_segs_5 = count_reads_segments(seg_end5s, "seg_end5s")
    num_reads_3, num_segs_3 = count_reads_segments(seg_end3s, "seg_end3s")
    require_equal("num_reads_5", num_reads_5, num_reads_3, "num_reads_3",
                  error_type=BadSegmentEndsError)
    require_equal("num_segs_5", num_segs_5, num_segs_3, "num_segs_3",
                  error_type=BadSegmentEndsError)
    if seg_ends_mask is not None:
        num_reads_mask, num_segs_mask = count_reads_segments(seg_ends_mask,
                                                             "seg_ends_mask")
        require_equal("num_reads_mask",
                      num_reads_mask,
                      num_reads_5,
                      "num_reads_5",
                      error_type=BadSegmentEndsError)
        require_equal("num_segs_mask",
                      num_segs_mask,
                      num_segs_5,
                      "num_segs_5",
                      error_type=BadSegmentEndsError)
    return num_reads_5, num_segs_5


def find_read_end5s(seg_end5s: np.ndarray, seg_ends_mask: np.ndarray | None):
    """ Find the 5' end of each read. """
    num_reads, num_segs = match_reads_segments(seg_end5s,
                                               seg_end5s,
                                               seg_ends_mask)
    if num_segs == 0:
        # There are 0 segments, so default to 5' end = 1.
        return np.ones(num_reads, dtype=seg_end5s.dtype)
    if seg_ends_mask is not None:
        mask_i, mask_j = np.nonzero(seg_ends_mask)
        assert mask_i.shape == mask_j.shape
        # Avoid modifying the original values when applying the mask.
        seg_end5s = seg_end5s.copy()
    else:
        mask_i = None
        mask_j = None
    if num_segs == 1:
        # There is exactly 1 segment, so use it.
        read_end5s = seg_end5s[:, 0]
    else:
        # Take the minimum over all segments.
        if seg_ends_mask is not None:
            # Replace masked values with the maximum possible value for
            # their data type so that they won't affect the minimum.
            seg_end5s[mask_i, mask_j] = np.iinfo(seg_end5s.dtype).max
        read_end5s = seg_end5s.min(axis=1)
    if seg_ends_mask is not None:
        # Fill reads where every 5' end is masked with 5' end = 1.
        read_end5s[seg_ends_mask.all(axis=1)] = 1
    return read_end5s


def find_read_end3s(seg_end3s: np.ndarray, seg_ends_mask: np.ndarray | None):
    """ Find the 3' end of each read. """
    num_reads, num_segs = match_reads_segments(seg_end3s,
                                               seg_end3s,
                                               seg_ends_mask)
    if num_segs == 0:
        # There are 0 segments, so default to 3' end = 0.
        return np.zeros(num_reads, dtype=seg_end3s.dtype)
    if seg_ends_mask is not None:
        mask_i, mask_j = np.nonzero(seg_ends_mask)
        assert mask_i.shape == mask_j.shape
        # Avoid modifying the original values when applying the mask.
        seg_end3s = seg_end3s.copy()
    else:
        mask_i = None
        mask_j = None
    if num_segs == 1:
        # There is exactly 1 segment, so use it.
        read_end3s = seg_end3s[:, 0]
    else:
        # Take the maximum over all segments.
        if seg_ends_mask is not None:
            # Replace masked values with the minimum possible value for
            # their data type so that they won't affect the maximum.
            seg_end3s[mask_i, mask_j] = np.iinfo(seg_end3s.dtype).min
        read_end3s = seg_end3s.max(axis=1)
    if seg_ends_mask is not None:
        # Fill reads where every 5' end is masked with 3' end = 0.
        read_end3s[seg_ends_mask.all(axis=1)] = 0
    return read_end3s


def merge_read_ends(read_end5s: np.ndarray, read_end3s: np.ndarray):
    """ Return the 5' and 3' ends as one 2D array. """
    ensure_same_length(read_end5s, read_end3s, "read_end5s", "read_end3s")
    return np.column_stack([read_end5s, read_end3s])


def sanitize_segment_ends(seg_end5s: np.ndarray,
                          seg_end3s: np.ndarray,
                          min_pos: int,
                          max_pos: int,
                          check_values: bool = True):
    """ Sanitize end coordinates.

    Parameters
    ----------
    seg_end5s: np.ndarray
        5' end coordinate of each segment in each read.
    seg_end3s: np.ndarray
        3' end coordinate of each segment in each read.
    min_pos: int
        Minimum allowed value of a position.
    max_pos: int
        Maximum allowed value of a position.
    check_values: bool = True
        Whether to check the bounds of the values, which is the most
        expensive operation in this function. Can be set to False if the
        only desired effect is to ensure the output is a positive, even
        number of arrays in the proper data type.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Sanitized end coordinates: encoded in the most efficient data
        type, and if `check_values` is True then all between `min_pos`
        and `max_pos` (inclusive).
    """
    num_reads, num_segs = match_reads_segments(seg_end5s, seg_end3s, None)
    if check_values:
        for i in range(num_segs):
            # All coordinates must be ≥ 1.
            ensure_order(seg_end5s[:, i],
                         np.broadcast_to(min_pos, num_reads),
                         f"segment {i} 5' coordinates",
                         f"minimum position ({min_pos})",
                         gt_eq=True)
            ensure_order(seg_end3s[:, i],
                         np.broadcast_to(min_pos - 1, num_reads),
                         f"segment {i} 3' coordinates",
                         f"minimum position - 1 ({min_pos - 1})",
                         gt_eq=True)
            # All coordinates must be ≤ max position
            ensure_order(seg_end5s[:, i],
                         np.broadcast_to(max_pos + 1, num_reads),
                         f"segment {i} 5' coordinates",
                         f"maximum position + 1 ({max_pos + 1})")
            ensure_order(seg_end3s[:, i],
                         np.broadcast_to(max_pos, num_reads),
                         f"segment {i} 3' coordinates",
                         f"maximum position ({max_pos})")
    # Convert all ends into arrays of the most efficient data type.
    # Conversion must be done after checking the values against max_pos.
    # Otherwise, values greater than the data type can accomodate would
    # be converted to their modulo relative to the maximum value of the
    # data type, which would mean checking the wrong values.
    dtype = fit_uint_type(max_pos)
    return np.asarray(seg_end5s, dtype), np.asarray(seg_end3s, dtype)


def _roll_unmasked(matrix: np.ndarray, mask: np.ndarray):
    assert isinstance(matrix, np.ndarray)
    assert isinstance(mask, np.ndarray)
    assert matrix.ndim == 2
    assert mask.ndim == 2
    assert matrix.shape == mask.shape
    # Get the row and column indices for valid entries.
    valid_rows, valid_cols = np.nonzero(np.logical_not(mask))
    # Count the valid entries in each row.
    valid_counts = np.bincount(valid_rows, minlength=matrix.shape[0])
    # Calculate what index the first valid value in each row would have
    # if all valid values were put in order into a flattened (1D) array.
    row_offsets = np.concatenate([[0], np.cumsum(valid_counts)[:-1]])
    # Compute the "global" (flattened) indices for the valid values.
    global_idx = np.arange(valid_rows.size)
    # Compute the "local" index within each row for every valid value.
    local_idx = global_idx - row_offsets[valid_rows]
    # For each valid value, compute its new local index after a roll by
    # one, which is (local_idx - 1) modulo the number of valid entries.
    rolled_local_idx = (local_idx - 1) % valid_counts[valid_rows]
    # Map the new local index back to a global index, which is the row’s
    # offset plus rolled_local_idx.
    rolled_global_idx = row_offsets[valid_rows] + rolled_local_idx
    # Get the rolled valid values using rolled_global_idx.
    assert np.array_equal(valid_rows, valid_rows[rolled_global_idx])
    rolled_values = matrix[valid_rows, valid_cols[rolled_global_idx]]
    # Place the rolled values into the result matrix at the positions of
    # valid entries.
    rolled_matrix = matrix.copy()
    rolled_matrix[valid_rows, valid_cols] = rolled_values
    return rolled_matrix


def sort_segment_ends(seg_end5s: np.ndarray,
                      seg_end3s: np.ndarray,
                      seg_ends_mask: np.ndarray | None):
    """ Sort the segment end coordinates and label the 3' end of each
    contiguous set of segments.

    Parameters
    ----------
    seg_end5s: np.ndarray
        5' end of each segment in each read.
    seg_end3s: np.ndarray
        3' end of each segment in each read.
    seg_ends_mask: np.ndarray | None
        Whether each pair of 5' and 3' ends is masked.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - Sorted 5' and 3' coordinates of the segments in each read
        - Labels of whether each coordinate is a 5' end of a segment
        - Labels of whether each coordinate is a 5' end of a contiguous
          segment
        - Labels of whether each coordinate is a 3' end of a contiguous
          segment
    """
    _, num_segs = match_reads_segments(seg_end5s, seg_end3s, seg_ends_mask)
    # Make the 5' ends 0-indexed and merge them with the 3' ends.
    if seg_end5s.size > 0 and seg_end5s.min() == 0:
        # Prevent overflow when subtracting 1.
        raise ValueError("All 5' ends must be ≥ 1, but got 0")
    seg_ends = np.hstack([seg_end5s - 1, seg_end3s])
    if seg_ends_mask is not None:
        seg_ends_mask = np.hstack([seg_ends_mask, seg_ends_mask])
        # Fill the masked 5' and 3' ends with 0.
        seg_ends[np.nonzero(seg_ends_mask)] = 0
    # Sort the coordinates in ascending order within each read.
    sort_order = seg_ends.argsort(axis=1, kind="stable")
    seg_ends_sorted = np.take_along_axis(seg_ends, sort_order, axis=1)
    # Label the sorted coordinates that correspond to segment 5' ends.
    is_seg_end5 = sort_order < num_segs
    # Label the 3' end of each contiguous set of segments, which occurs
    # when the cumulative numbers of 5' and 3' segment ends are equal.
    is_contig_end3 = np.logical_not(np.cumsum(np.where(is_seg_end5, 1, -1),
                                              axis=1))
    if seg_ends_mask is not None:
        # Apply the mask that seg_ends_sorted has to is_contig_end3
        # and is_contig_end5.
        mask_sorted = np.take_along_axis(seg_ends_mask, sort_order, axis=1)
        is_contig_end3[np.nonzero(mask_sorted)] = False
        is_contig_end5 = _roll_unmasked(is_contig_end3, mask_sorted)
    else:
        is_contig_end5 = np.roll(is_contig_end3, 1, axis=1)
    # In every read, the number of 5' and 3' ends of contiguous segments
    # must be equal.
    assert np.array_equal(np.count_nonzero(is_contig_end5, axis=1),
                          np.count_nonzero(is_contig_end3, axis=1))
    return seg_ends_sorted, is_contig_end5, is_contig_end3


def find_contiguous_reads(seg_end5s: np.ndarray,
                          seg_end3s: np.ndarray,
                          seg_ends_mask: np.ndarray | None):
    """ Whether the segments of each read are contiguous. """
    num_reads, num_segs = match_reads_segments(seg_end5s,
                                               seg_end3s,
                                               seg_ends_mask)
    if num_segs <= 1:
        # Reads with fewer than 2 segments cannot be discontiguous.
        return np.ones(num_reads, dtype=bool)
    # For contiguous reads, only the last coordinate (when sorted) will
    # be the end of a contiguous segment.
    _, _, is_contig_end3 = sort_segment_ends(seg_end5s,
                                             seg_end3s,
                                             seg_ends_mask)
    return np.logical_not(np.count_nonzero(is_contig_end3[:, :-1], axis=1))


def simulate_segment_ends(uniq_end5s: np.ndarray,
                          uniq_end3s: np.ndarray,
                          p_ends: np.ndarray,
                          num_reads: int,
                          read_length: int = 0,
                          p_rev: float = 0.5):
    """ Simulate segment end coordinates from their probabilities.

    Parameters
    ----------
    uniq_end5s: np.ndarray
        Unique read 5' end coordinates.
    uniq_end3s: np.ndarray
        Unique read 3' end coordinates.
    p_ends: np.ndarray
        Probability of each set of unique end coordinates.
    num_reads: int
        Number of reads to simulate.
    read_length: int = 0
        If == 0, then generate single-end reads (1 segment per read);
        if > 0, then generate paired-end reads (2 segments per read)
        with at most this number of base calls in each segment.
    p_rev: float = 0.5
        For paired-end reads, the probability that mate 1 aligns in the

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        5' and 3' segment end coordinates of the simulated reads.
    """
    require_atleast("num_reads", num_reads, 0)
    if num_reads == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    find_dims([(NUM_UNIQ,), (NUM_UNIQ,), (NUM_UNIQ,)],
              [uniq_end5s, uniq_end3s, p_ends],
              ["end5s", "end3s", "p_ends"],
              nonzero=[NUM_UNIQ])
    # Drop pairs of 5'/3' ends where the 5' is greater than the 3'.
    valid_ends = np.less_equal(uniq_end5s, uniq_end3s)
    num_ends = np.count_nonzero(valid_ends)
    if num_ends == 0:
        raise ValueError("Got 0 pairs of 5'/3' ends for which 5' > 3'")
    elif num_ends < valid_ends.size:
        logger.warning(f"Got {valid_ends.size - num_ends} pairs of 5'/3' ends "
                       f"for which 5' > 3'")
        uniq_end5s = uniq_end5s[valid_ends]
        uniq_end3s = uniq_end3s[valid_ends]
        p_ends = p_ends[valid_ends]
        total_p_ends = p_ends.sum()
        if total_p_ends > 0:
            p_ends = p_ends / total_p_ends
        else:
            # num_ends is not 0.
            p_ends = np.full(num_ends, 1. / num_ends)
    elif not np.isclose(p_ends.sum(), 1.):
        raise ValueError(f"p_ends must sum to 1, but got {p_ends.sum()}")
    # Choose reads based on their probabilities.
    indexes = rng.choice(num_ends, num_reads, p=p_ends)
    read_end5s = uniq_end5s[indexes]
    read_end3s = uniq_end3s[indexes]
    require_atleast("read_length", read_length, 0)
    if read_length > 0:
        diff = read_length - 1
        seg_end5s = np.stack([read_end5s,
                              np.maximum(read_end3s - diff, read_end5s)],
                             axis=1)
        seg_end3s = np.stack([np.minimum(read_end5s + diff, read_end3s),
                              read_end3s],
                             axis=1)
        # For reverse reads, swap mates 1 and 2.
        reverse = rng.random(num_reads) < p_rev
        indices = np.stack([reverse, ~reverse], axis=1, dtype=int)
        seg_end5s = np.take_along_axis(seg_end5s, indices, axis=1)
        seg_end3s = np.take_along_axis(seg_end3s, indices, axis=1)
    else:
        seg_end5s = read_end5s[:, np.newaxis]
        seg_end3s = read_end3s[:, np.newaxis]
    return seg_end5s, seg_end3s


class EndCoords(object):
    """ Collection of 5' and 3' segment end coordinates. """

    def __init__(self, *,
                 region: Region,
                 seg_end5s: np.ndarray,
                 seg_end3s: np.ndarray,
                 sanitize: bool = True,
                 **kwargs):
        super().__init__(**kwargs)
        # Validate and store the segment end coordinates.
        self.seg_end5s, self.seg_end3s = sanitize_segment_ends(seg_end5s,
                                                               seg_end3s,
                                                               region.end5,
                                                               region.end3,
                                                               sanitize)

    @cached_property
    def num_reads(self):
        """ Number of reads. """
        num_reads, _ = match_reads_segments(self.seg_end5s,
                                            self.seg_end3s,
                                            self.seg_ends_mask)
        return num_reads

    @cached_property
    def num_segments(self):
        """ Number of segments in each read. """
        _, num_segs = match_reads_segments(self.seg_end5s,
                                           self.seg_end3s,
                                           self.seg_ends_mask)
        return num_segs

    @cached_property
    def seg_ends_mask(self) -> np.ndarray | None:
        """ Whether each pair of 5' and 3' ends is masked out. """
        # Store seg_ends_mask as a cached_property, not an instance
        # attribute, so that it never consumes space in the file system.
        seg_ends_mask = self.seg_end5s > self.seg_end3s
        return seg_ends_mask if np.any(seg_ends_mask) else None

    @cached_property
    def read_end5s(self):
        """ 5' end of each read. """
        return find_read_end5s(self.seg_end5s, self.seg_ends_mask)

    @cached_property
    def read_end3s(self):
        """ 3' end of each read. """
        return find_read_end3s(self.seg_end3s, self.seg_ends_mask)

    @cached_property
    def read_lengths(self):
        """ Length of each read. """
        return self.read_end3s - self.read_end5s + 1

    @property
    def pos_dtype(self):
        """ Data type for positions. """
        if self.seg_end5s.dtype is not self.seg_end3s.dtype:
            raise ValueError("Data types differ for 5' and 3' ends")
        return self.seg_end5s.dtype

    @cached_property
    def contiguous(self):
        """ Whether the segments of each read are contiguous. """
        return find_contiguous_reads(self.seg_end5s,
                                     self.seg_end3s,
                                     self.seg_ends_mask)

    @cached_property
    def num_contiguous(self):
        """ Number of contiguous reads. """
        return np.count_nonzero(self.contiguous)

    @cached_property
    def num_discontiguous(self):
        """ Number of discontiguous reads. """
        return self.num_reads - self.num_contiguous
