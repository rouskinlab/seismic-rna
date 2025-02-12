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
    require_atleast("num_segs", num_segs, 1,
                    error_type=BadSegmentEndsError)
    return num_reads, num_segs


def match_reads_segments(seg_end5s: np.ndarray, seg_end3s: np.ndarray):
    """ Number of segments for the given end coordinates. """
    num_reads_5, num_segs_5 = count_reads_segments(seg_end5s, "seg_end5s")
    num_reads_3, num_segs_3 = count_reads_segments(seg_end3s, "seg_end3s")
    require_equal("num_reads_5", num_reads_5, num_reads_3, "num_reads_3",
                  error_type=BadSegmentEndsError)
    require_equal("num_segs_5", num_segs_5, num_segs_3, "num_segs_3",
                  error_type=BadSegmentEndsError)
    return num_reads_5, num_segs_5


def mask_segment_ends(seg_end5s: np.ndarray, seg_end3s: np.ndarray):
    """ Mask segments with no coverage (5' end > 3' end). """
    match_reads_segments(seg_end5s, seg_end3s)
    no_coverage = seg_end5s > seg_end3s
    if np.any(no_coverage):
        seg_end5s = np.ma.masked_array(seg_end5s, no_coverage)
        seg_end3s = np.ma.masked_array(seg_end3s, no_coverage)
    return seg_end5s, seg_end3s


def merge_segment_ends(seg_end5s: np.ndarray,
                       seg_end3s: np.ndarray,
                       fill_value: int | None = None):
    """ """
    match_reads_segments(seg_end5s, seg_end3s)
    if np.ma.isarray(seg_end5s) or np.ma.isarray(seg_end3s):
        seg_ends = np.ma.hstack([seg_end5s, seg_end3s])
        if fill_value is not None:
            seg_ends.fill_value = 0
    else:
        seg_ends = np.hstack([seg_end5s, seg_end3s])
    return seg_ends


def _check_no_coverage_reads(seg_ends: np.ndarray, what: str = "seg_ends"):
    num_reads, num_segs = count_reads_segments(seg_ends, what)
    if np.ma.isarray(seg_ends):
        # Check if any reads have no coverage (are completely masked).
        no_coverage = np.ma.getmaskarray(seg_ends).all(axis=1)
        if np.any(no_coverage):
            raise BadSegmentEndsError(
                f"{np.count_nonzero(no_coverage)} reads have no coverage "
            )
    return num_reads, num_segs


def find_read_end5s(seg_end5s: np.ndarray):
    """ Find the 5' end of each read. """
    num_reads, num_segs = _check_no_coverage_reads(seg_end5s, "seg_end5s")
    return np.asarray(seg_end5s.min(axis=1)
                      if num_segs > 1
                      else seg_end5s[:, 0])


def find_read_end3s(seg_end3s: np.ndarray):
    """ Find the 3' end of each read. """
    num_reads, num_segs = _check_no_coverage_reads(seg_end3s, "seg_end3s")
    return np.asarray(seg_end3s.max(axis=1)
                      if num_segs > 1
                      else seg_end3s[:, 0])


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
    num_reads, num_segs = match_reads_segments(seg_end5s, seg_end3s)
    require_atleast("num_segs", num_segs, 1, error_type=BadSegmentEndsError)
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
    return np.asanyarray(seg_end5s, dtype), np.asanyarray(seg_end3s, dtype)


def sort_segment_ends(seg_end5s: np.ndarray,
                      seg_end3s: np.ndarray,
                      zero_indexed: bool = True,
                      fill_mask: bool = False):
    """ Sort the segment end coordinates and label the 3' end of each
    contiguous set of segments.

    Parameters
    ----------
    seg_end5s: np.ndarray
        5' end of each segment in each read; may be masked.
    seg_end3s: np.ndarray
        3' end of each segment in each read; may be masked.
    zero_indexed: bool = True
        In the return array, make the 5' ends 0-indexed; if False, then
        they will be 1-indexed (like the input).
    fill_mask: bool = False
        If `seg_end5s` or `seg_end3s` is a masked array, then return a
        regular array with all masked coordinates set to 0 (or 1 for
        5' ends if `one_indexed` is True) rather than a masked array.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        - Sorted 5' and 3' coordinates of the segments in each read
        - Labels of whether each coordinate is a 5' end of a segment
        - Labels of whether each coordinate is a 3' end of a contiguous
          segment
    """
    _, num_segs = match_reads_segments(seg_end5s, seg_end3s)
    # Make the 5' ends 0-indexed and join the coordinates horizontally.
    seg_ends = merge_segment_ends(seg_end5s - 1, seg_end3s)
    # Sort the coordinates in ascending order within each read.
    kwargs = dict(fill_value=0) if np.ma.isarray(seg_ends) else dict()
    sort_order = seg_ends.argsort(axis=1, kind="stable", **kwargs)
    if np.ma.isarray(seg_ends) and fill_mask:
        # Fill the masked values with zero and remove the mask.
        mask = seg_ends.mask
        seg_ends = seg_ends.data
        seg_ends[mask] = 0
    seg_ends_sorted = np.take_along_axis(seg_ends, sort_order, axis=1)
    # Label the sorted coordinates that correspond to segment 5' ends.
    is_seg_end5 = sort_order < num_segs
    # Label the 3' end of each contiguous set of segments, which occurs
    # when the cumulative numbers of 5' and 3' segment ends are equal.
    is_contig_end3 = np.logical_not(np.cumsum(np.where(is_seg_end5, 1, -1),
                                              axis=1))
    if np.ma.isarray(seg_ends_sorted):
        # Mask is_seg_end5 and is_contig_end3 like seg_ends_sorted.
        is_seg_end5 = np.ma.masked_array(is_seg_end5,
                                         seg_ends_sorted.mask)
        is_contig_end3 = np.ma.masked_array(is_contig_end3,
                                            seg_ends_sorted.mask)
    if not zero_indexed:
        # Make the 5' coordinates 1-indexed again.
        seg_ends_sorted[np.nonzero(is_seg_end5)] += 1
    return seg_ends_sorted, is_seg_end5, is_contig_end3


def find_contiguous_reads(seg_end5s: np.ndarray, seg_end3s: np.ndarray):
    """ Whether the segments of each read are contiguous. """
    num_reads, num_segs = match_reads_segments(seg_end5s, seg_end3s)
    if num_segs == 1:
        # Single-segment reads must be contiguous.
        return np.ones(num_reads, dtype=bool)
    # For contiguous reads, only the last coordinate (when sorted) will
    # be the end of a contiguous segment.
    _, _, is_contig_end3 = sort_segment_ends(seg_end5s,
                                             seg_end3s,
                                             zero_indexed=True)
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
        self._end5s, self._end3s = sanitize_segment_ends(seg_end5s,
                                                         seg_end3s,
                                                         region.end5,
                                                         region.end3,
                                                         sanitize)

    @cached_property
    def num_reads(self):
        """ Number of reads. """
        num_reads, _ = match_reads_segments(self._end5s, self._end3s)
        return num_reads

    @cached_property
    def num_segments(self):
        """ Number of segments in each read. """
        _, num_segs = match_reads_segments(self._end5s, self._end3s)
        return num_segs

    @cached_property
    def _seg_ends(self):
        """ 5' and 3' ends of each segment in each read. """
        return mask_segment_ends(self._end5s, self._end3s)

    @property
    def seg_end5s(self):
        """ 5' end of each segment in each read. """
        seg_end5s, _ = self._seg_ends
        return seg_end5s

    @property
    def seg_end3s(self):
        """ 3' end of each segment in each read. """
        _, seg_end3s = self._seg_ends
        return seg_end3s

    @cached_property
    def read_end5s(self):
        """ 5' end of each read. """
        return find_read_end5s(self.seg_end5s)

    @cached_property
    def read_end3s(self):
        """ 3' end of each read. """
        return find_read_end3s(self.seg_end3s)

    @cached_property
    def read_lengths(self):
        """ Length of each read. """
        return self.read_end3s - self.read_end5s + 1

    @property
    def pos_dtype(self):
        """ Data type for positions. """
        if self._end5s.dtype is not self._end3s.dtype:
            raise ValueError("Data types differ for 5' and 3' ends")
        return self._end5s.dtype

    @cached_property
    def contiguous(self):
        """ Whether the segments of each read are contiguous. """
        return find_contiguous_reads(self.seg_end5s, self.seg_end3s)

    @cached_property
    def num_contiguous(self):
        """ Number of contiguous reads. """
        return np.count_nonzero(self.contiguous)

    @cached_property
    def num_discontiguous(self):
        """ Number of discontiguous reads. """
        return self.num_reads - self.num_contiguous
