import re
from collections import defaultdict, namedtuple
from functools import cached_property, reduce
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .refs import RefSeqs
from .xna import BASEA, BASEC, BASEN, DNA
from ..logs import logger

# Names of the region index levels.
POS_NAME = "Position"
BASE_NAME = "Base"
SEQ_INDEX_NAMES = POS_NAME, BASE_NAME
FULL_NAME = "full"

# Fields in the regions file.
FIELD_REF = "Reference"
FIELD_REG = "Region"
FIELD_END5 = "5' End"
FIELD_END3 = "3' End"
FIELD_PFWD = "Forward Primer"
FIELD_PREV = "Reverse Primer"

RegionTuple = namedtuple("RegionTuple", ("pos5", "pos3"))


def get_reg_coords_primers(regs_file: Path):
    """
    Parse a file defining each region by the name of its reference and
    either its 5' and 3' coordinates or its forward and reverse primer
    sequences. Return one map from each reference and 5'/3' coordinate
    pair to the name of the corresponding region, and another from each
    reference and primer pair to the name of the corresponding region.

    Parameters
    ----------
    regs_file: Path
        CSV file of a table that defines the regions. The table must
        have columns labeled "Reference", "Region", "5' End", "3' End",
        "Forward Primer", and "Reverse Primer". Others are ignored.

    Returns
    -------
    tuple[dict[tuple[str, int, int], str],
          dict[tuple[str, DNA, DNA], str]]
        Two mappings, the first from (ref name, 5' coord, 3' coord) to
        each region, the second from (ref name, fwd primer, rev primer)
        to each region. If the region is named in the "Region" column
        of the table, then that name will be used as the region name.
        Otherwise, the region name will be an empty string.
    """

    # Initialize dictionaries mapping references and coordinates/primers
    # to region names.
    coords: dict[tuple[str, int, int], str] = dict()
    primers: dict[tuple[str, DNA, DNA], str] = dict()

    def map_reg(mapping: dict[tuple, str], key: tuple, value: str):
        """ Add one region to the map if not already present. """
        # Check whether the mapping already contains the key.
        try:
            prev = mapping[key]
        except KeyError:
            # The mapping does not already contain the key: add it.
            mapping[key] = value
        else:
            # Check whether the value already mapped by the key matches
            # the value given currently.
            if prev == value:
                # If so, then warn about it.
                logger.warning(
                    f"Key {key} mapped to {repr(value)} multiple times"
                )
            else:
                # If not, then raise an error because it is ambiguous
                # which value to use.
                raise ValueError(
                    f"Key {key} mapped to {repr(prev)} and to {repr(value)}"
                )

    # Read every row of the regions file.
    regions = pd.read_csv(regs_file)
    lines = zip(regions[FIELD_REF], regions[FIELD_REG],
                regions[FIELD_END5], regions[FIELD_END3],
                regions[FIELD_PFWD], regions[FIELD_PREV])
    for i, (ref, reg, end5, end3, fwd, rev) in enumerate(lines, start=1):
        try:
            # The reference name must have a value.
            if pd.isnull(ref):
                raise ValueError(f"Missing {FIELD_REF}")
            else:
                ref = str(ref)
            # The region name may be left blank.
            reg = "" if pd.isnull(reg) else str(reg)
            if reg.lower() == FULL_NAME.lower():
                raise ValueError(
                    f"A region cannot be given the name {repr(FULL_NAME)}, "
                    "which reserved for when a reference is automatically "
                    "given a full region in the absence of coordinates/primers"
                )
            # Check whether coordinates or primers were given.
            has_coords = not (pd.isnull(end5) or pd.isnull(end3))
            has_primers = not (pd.isnull(fwd) or pd.isnull(rev))
            if has_coords and has_primers:
                raise ValueError(f"Got both coordinates ({end5}, {end3}) "
                                 f"and primers ({fwd}, {rev})")
            elif has_coords:
                # Map the reference and coordinates to the region.
                map_reg(coords, (ref, int(end5), int(end3)), reg)
            elif has_primers:
                # Map the reference and primers to the region.
                map_reg(primers, (ref, DNA(fwd), DNA(rev)), reg)
            else:
                raise ValueError("Got neither coordinates nor primers")
        except Exception as error:
            logger.error(error)
    return coords, primers


def seq_pos_to_index(seq: DNA, positions: Sequence[int], start: int):
    """
    Convert a sequence and positions to indexes, where each index is a
    tuple of (position, base).

    Parameters
    ----------
    seq: DNA
        DNA sequence.
    positions: Sequence[int]
        Positions of the sequence from which to build the index. Every
        position must be an integer ≥ `start`.
    start: int
        Numerical position to assign to the first base in the sequence.
        Must be a positive integer.

    Returns
    -------
    pd.MultiIndex
        MultiIndex of the same length as positions where each index is a
        tuple of (position, base).
    """
    if start < 1:
        raise ValueError(f"The start position must be ≥ 1, but got {start}")
    # Cast positions to a NumPy integer array.
    pos = np.asarray(positions, dtype=int)
    # Validate the positions.
    if pos.size > 0 and np.min(pos) < start:
        raise ValueError(
            f"All positions must be ≥ start ({start}), but got {positions}")
    end = start + len(seq) - 1
    if pos.size > 0 and np.max(pos) > end:
        raise ValueError(
            f"All positions must be ≤ end ({end}), but got {positions}")
    # Create a 2-level MultiIndex from the positions and the bases in
    # the sequence at those positions.
    index = pd.MultiIndex.from_arrays([pos, seq.array[pos - start]],
                                      names=SEQ_INDEX_NAMES)
    if index.has_duplicates:
        raise ValueError(f"Duplicated positions: {positions}")
    if not index.is_monotonic_increasing:
        raise ValueError(f"Unsorted positions: {positions}")
    return index


def verify_index_names(index: pd.MultiIndex):
    """ Verify that the names of the index are correct. """
    if tuple(index.names) != SEQ_INDEX_NAMES:
        raise ValueError(f"Expected index with names {SEQ_INDEX_NAMES}, "
                         f"but got {index.names}")


def index_to_pos(index: pd.MultiIndex):
    """ Get the positions from a MultiIndex of (pos, base) pairs. """
    verify_index_names(index)
    positions = index.get_level_values(POS_NAME)
    if positions.has_duplicates:
        raise ValueError(f"Index has duplicate positions:\n{positions}")
    if not positions.is_monotonic_increasing:
        raise ValueError(f"Positions in index are not sorted:\n{positions}")
    return positions.values


def index_to_seq(index: pd.MultiIndex, allow_gaps: bool = False):
    """ Get the DNA sequence from a MultiIndex of (pos, base) pairs. """
    # Get the numeric positions and verify that there is at least one.
    if index.size == 0:
        # Checks for sorted and contiguous positions will fail if there
        # are no positions. Just return an empty sequence.
        return DNA("")
    pos = index_to_pos(index)
    # Verify that the positions are sorted and contiguous.
    if not (allow_gaps or np.array_equal(pos, np.arange(pos[0], pos[-1] + 1))):
        raise ValueError("A sequence cannot be assembled from an index with "
                         f"missing positions:\n{pos}")
    # Join the bases in the index and convert them to a DNA sequence.
    return DNA("".join(index.get_level_values(BASE_NAME)))


def get_shared_index(indexes: Iterable[pd.MultiIndex], empty_ok: bool = False):
    """ Get the shared index among all those given, as follows:

    - If indexes contains no elements and empty_ok is True, then return
      an empty MultiIndex with levels named 'Positions' and 'Base'.
    - If indexes contains one element or multiple identical elements,
      and each has two levels named 'Positions' and 'Base', then return
      the first element.
    - Otherwise, raise an error.

    Parameters
    ----------
    indexes: Iterable[pandas.MultiIndex]
        Indexes to compare.
    empty_ok: bool = False
        If given no indexes, then default to an empty index (if True)
        or raise a ValueError (if False).

    Returns
    -------
    pandas.MultiIndex
        The shared index.
    """
    # Ensure indexes is a list-like object.
    indexes = [index.copy(deep=True).remove_unused_levels()
               for index in indexes]
    try:
        # Get the first index.
        index = indexes[0]
    except IndexError:
        if empty_ok:
            # Return an empty MultiIndex.
            return pd.MultiIndex.from_arrays([np.array([], dtype=int),
                                              np.array([], dtype=str)],
                                             names=SEQ_INDEX_NAMES)
        raise ValueError("No indexes were given")
    # Ensure the first index is a MultiIndex with levels 'Positions' and
    # 'Base'.
    if tuple(index.names) != SEQ_INDEX_NAMES:
        raise ValueError(
            f"Expected levels named {SEQ_INDEX_NAMES}, but got {index.names}"
        )
    # Ensure all indexes are identical.
    for i, other in enumerate(indexes[1:], start=1):
        if not index.equals(other):
            raise ValueError(f"Indexes 0 and {i} differ: {index} ≠ {other}")
        if index.nlevels != other.nlevels:
            raise ValueError(f"Indexes 0 and {i} have {index.nlevels} and "
                             f"{other.nlevels} levels, respectively")
        for level in range(index.nlevels):
            # Compare sets because order does not matter.
            ilevel = index.levels[level]
            olevel = other.levels[level]
            if set(ilevel) != set(olevel):
                raise ValueError(f"Indexes 0 and {i} level {level} is "
                                 f"{ilevel} and {olevel}, respectively")
    return index


def window_to_margins(window: int):
    """ Compute the 5' and 3' margins from the size of the window. """
    if not isinstance(window, int):
        raise TypeError(f"window must be int, but got {type(window).__name__}")
    if window <= 0:
        raise ValueError(f"window must be ≥ 1, but got {window}")
    margin3 = window // 2
    margin5 = window - (margin3 + 1)
    return margin5, margin3


def iter_windows(*series: pd.Series,
                 size: int,
                 min_count: int = 1,
                 include_nan: bool = False):
    # Determine the index; the index of all series must match.
    index = get_shared_index((s.index for s in series), empty_ok=True)
    # Calculate the 5' and 3' margins.
    margin5, margin3 = window_to_margins(size)
    if min_count > size:
        logger.warning(f"min_count ({min_count}) is > window size ({size}), "
                       "so no positions will have enough data")
    # Yield each window from each series.
    for center in index_to_pos(index):
        # Determine the 5' and 3' ends of the window.
        end5 = center - margin5
        end3 = center + margin3
        # Slice the window from each series.
        windows = tuple(s.loc[end5: end3] for s in series)
        if include_nan:
            if min_count > 0 and windows[0].size < min_count:
                # If there are insufficient positions in this window,
                # then skip it.
                continue
        else:
            # Determine which positions have no NaN value in any series.
            no_nan = np.logical_not(reduce(np.logical_or,
                                           map(np.isnan, windows)))
            if min_count > 0 and np.count_nonzero(no_nan) < min_count:
                # If there are insufficient positions where no series
                # has a NaN value, then skip this window.
                continue
            if not np.all(no_nan):
                # Keep only positions with no NaN value in any series.
                windows = tuple(w[no_nan] for w in windows)
        # Yield the window from each series.
        yield center, windows


def hyphenate_ends(end5: int, end3: int):
    """ Return the 5' and 3' ends as a hyphenated string.

    Parameters
    ----------
    end5: int
        5' end (1-indexed)
    end3: int
        3' end (1-indexed)

    Returns
    -------
    str
        Hyphenated 5' and 3' ends
    """
    return f"{end5}-{end3}"


class Region(object):
    """ Region of a reference sequence between two coordinates. """

    MASK_POLYA = "pos-polya"
    MASK_GU = "pos-gu"
    MASK_LIST = "pos-list"

    def __init__(self,
                 ref: str,
                 seq: DNA, *,
                 seq5: int = 1,
                 reflen: int | None = None,
                 end5: int | None = None,
                 end3: int | None = None,
                 name: str | None = None):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence.
        seq: DNA
            The full reference sequence or a part of it.
        seq5: int = 1
            Positional number to assign the 5' end of the given part of
            the reference sequence. Must be ≥ 1.
        reflen: int | None = None
            Length of the full reference sequence. Must be ≥ 0. If None,
            defaults to the 3' end position in `seq`.
        end5: int | None = None
            Coordinate of the reference sequence at which the region's
            5' end is located.
        end3: int | None = None
            Coordinate of the reference sequence at which the region's
            3' end is located.
        name: str | None = None
            Name of the region. If None, defaults to `self.range`.
        """
        self.ref = ref
        if seq5 < 1:
            raise ValueError(f"seq5 must be ≥ 1, but got {seq5}")
        # Compute the 3' end position of the given sequence.
        seq3 = seq5 + len(seq) - 1
        if reflen is None:
            # Default to the 3' end position of the given sequence.
            reflen = seq3
        elif reflen < 0:
            raise ValueError(f"reflen must be ≥ 0, but got {reflen}")
        elif reflen < seq3:
            raise ValueError(f"The 3' end of the given sequence is {seq3}, "
                             f"but the full reference is only {reflen} nt")
        # Validate that the region lies within the given sequence.
        if end5 is not None:
            if end5 < seq5:
                raise ValueError(f"Need end5 ≥ seq5, but got {end5} < {seq5}")
            self.end5 = end5
        else:
            self.end5 = seq5
        if end3 is not None:
            if end3 > seq3:
                raise ValueError(f"Need end3 ≤ seq3, but got {end3} > {seq3}")
            self.end3 = end3
        else:
            self.end3 = seq3
        if self.end5 > self.end3 + 1:
            raise ValueError("Need end5 ≤ end3 + 1, "
                             f"but got {self.end5} > {self.end3 + 1}")
        # Determine the sequence of the region and whether it is the
        # full reference sequence.
        self.seq = DNA(seq[self.end5 - seq5: self.end3 - (seq5 - 1)])
        self.full = self.end5 == 1 and self.end3 == reflen
        # Assign the name of the region.
        if name is None:
            # Default to "full" if the region spans the full reference
            # sequence and neither the 5' nor 3' coordinate was given.
            # Otherwise, default to the hyphenated coordinates.
            self.name = (FULL_NAME
                         if self.full and end5 is None and end3 is None
                         else self.hyphen)
        elif isinstance(name, str):
            # Use the given name unless it is an empty string, in which
            # case default to the hyphenated coordinates.
            self.name = name if name else self.hyphen
        else:
            raise TypeError(f"name must be a str, but got {repr(name)}")
        # Initialize an empty set of masks.
        self._masks: dict[str, np.ndarray] = dict()

    @cached_property
    def length(self):
        """ Length of the entire region. """
        return self.end3 - self.end5 + 1

    @cached_property
    def coord(self):
        """ Tuple of the 5' and 3' coordinates. """
        return self.end5, self.end3

    @cached_property
    def hyphen(self):
        return hyphenate_ends(self.end5, self.end3)

    @cached_property
    def ref_reg(self):
        return f"{self.ref}__{self.name}"

    def to_dict(self):
        return dict(ref=self.ref,
                    seq=self.seq,
                    reg=self.name,
                    end5=self.end5,
                    end3=self.end3)

    @property
    def mask_names(self):
        """ Names of the masks. """
        return list(self._masks)

    @property
    def masked_int(self) -> np.ndarray:
        """ Masked positions as integers. """
        # Do not cache this method since self._masks can change.
        return reduce(np.union1d, self._masks.values(), np.array([], int))

    @property
    def masked_zero(self) -> np.ndarray:
        """ Masked positions as integers (0-indexed with respect to the
        first position in the region). """
        # Do not cache this method since self.masked_int can change.
        return self.masked_int - self.end5

    @property
    def masked_bool(self) -> np.ndarray:
        """ Masked positions as a boolean array. """
        # Do not cache this method since self.masked_int can change.
        return np.isin(self.range_int, self.masked_int)

    @property
    def unmasked_bool(self) -> np.ndarray:
        """ Unmasked positions as a boolean array. """
        # Do not cache this method since self.masked_bool can change.
        return np.logical_not(self.masked_bool)

    @property
    def unmasked_int(self) -> np.ndarray:
        """ Unmasked positions as integers (1-indexed). """
        # Do not cache this method since self.unmasked_bool can change.
        return self.range_int[self.unmasked_bool]

    @property
    def unmasked_zero(self) -> np.ndarray:
        """ Unmasked positions as integers (0-indexed with respect to
        the first position in the region). """
        # Do not cache this method since self.unmasked_int can change.
        return self.unmasked_int - self.end5

    @property
    def unmasked(self):
        """ Index of unmasked positions in the region. """
        # Do not cache this method since self.unmasked_int can change.
        return seq_pos_to_index(self.seq, self.unmasked_int, self.end5)

    @cached_property
    def range_int(self):
        """ All positions in the region as integers. """
        return np.arange(self.end5, self.end3 + 1, dtype=int)

    @cached_property
    def range_one(self):
        """ All 1-indexed positions in the region as integers. """
        return np.arange(1, self.length + 1, dtype=int)

    @cached_property
    def range(self):
        """ Index of all positions in the region. """
        return seq_pos_to_index(self.seq, self.range_int, self.end5)

    @property
    def size(self):
        """ Number of relevant positions in the region. """
        return self.length - self.masked_int.size

    def copy(self, masks: bool = True):
        """ Return an identical region. """
        copied = self.__class__(ref=self.ref,
                                seq=self.seq,
                                seq5=self.end5,
                                reflen=(self.end3 if self.full
                                        else self.end3 + 1),
                                end5=self.end5,
                                end3=self.end3,
                                name=self.name)
        if masks:
            for name, masked in self._masks.items():
                copied.add_mask(name, masked)
        return copied

    def get_mask(self, name: str):
        """ Get the positions masked under the given name. """
        return self._masks[name]

    def add_mask(self,
                 name: str,
                 positions: Iterable[int],
                 complement: bool = False):
        """ Mask the integer positions in the array `positions`.

        Parameters
        ----------
        name: str
            Name of the mask.
        positions: Iterable[int]
            Positions to mask (1-indexed).
        complement: bool = False
            If True, then leave only positions in `positions` unmasked.
        """
        # Convert positions to a NumPy integer array.
        p = np.unique(np.asarray(list(positions), dtype=int))
        # Check for positions outside the region.
        if np.any(p < self.end5) or np.any(p > self.end3):
            out = p[np.logical_or(p < self.end5, p > self.end3)]
            raise ValueError(f"Got positions to mask outside of {self}: {out}")
        if complement:
            # Mask all positions except those listed.
            p = np.setdiff1d(self.range_int, p, assume_unique=True)
        # Record the positions that have not already been masked.
        p = np.setdiff1d(p, self.masked_int, assume_unique=True)
        # Combine the positions with any existing positions.
        existing = self._masks.get(name)
        if existing is not None:
            p = np.union1d(p, existing)
        self._masks[name] = p
        # Return self to allow chaining repeated calls to add_mask(),
        # e.g. region.add_mask("a", [1]).add_mask("b", [2], True).
        return self

    def remove_mask(self, name: str, missing_ok: bool = False):
        """ Remove the specified mask from the region. """
        try:
            self._masks.pop(name)
        except KeyError:
            if not missing_ok:
                raise

    def _find_gu(self) -> np.ndarray:
        """ Array of each position whose base is neither A nor C. """
        # Mark whether each position is neither A nor C.
        gu_pos = np.logical_and(self.seq.array != BASEA,
                                self.seq.array != BASEC)
        # Return the integer positions.
        return self.range_int[gu_pos]

    def mask_gu(self):
        """ Mask positions whose base is neither A nor C. """
        return self.add_mask(self.MASK_GU, self._find_gu())

    def _find_polya(self, min_length: int) -> np.ndarray:
        """ Array of each position within a stretch of `min_length` or
        more consecutive adenines. """
        if min_length < 0:
            raise ValueError(f"min_length must be ≥ 0, but got {min_length}")
        # Initialize a list of 0-indexed positions in poly(A) sequences.
        polya_pos = list()
        if min_length > 0:
            # Generate a pattern that matches stretches of consecutive
            # adenines that are at least as long as min_length.
            polya_pattern = "%c{%d,}" % (BASEA, min_length)
            # Add the 0-indexed positions in every poly(A) sequence.
            for polya in re.finditer(polya_pattern, str(self.seq)):
                polya_pos.extend(range(polya.start(), polya.end()))
        # Convert the positions to an array with natural indexing.
        return np.array(polya_pos, dtype=int) + self.end5

    def mask_polya(self, min_length: int):
        """ Mask poly(A) stretches with length ≥ `min_length`. """
        return self.add_mask(self.MASK_POLYA, self._find_polya(min_length))

    def mask_list(self, pos: Iterable[int]):
        """ Mask a list of positions. """
        return self.add_mask(self.MASK_LIST, pos)

    def subregion(self,
                  end5: int | None = None,
                  end3: int | None = None,
                  name: str | None = None):
        """ Return a new region from part of this region. """
        return self.__class__(self.ref,
                              self.seq,
                              seq5=self.end5,
                              end5=end5 if end5 is not None else self.end5,
                              end3=end3 if end3 is not None else self.end3,
                              name=(FULL_NAME if (self.full
                                                  and end5 is None
                                                  and end3 is None
                                                  and name is None)
                                    else name))

    def renumber_from(self, seq5: int, name: str | None = None):
        """ Return a new region renumbered starting from a position.

        Parameters
        ----------
        seq5: int
            Position from which to start the new numbering system.
        name: str | None = None
            Name of the renumbered region.

        Returns
        -------
        Region
            Region with renumbered positions.
        """
        renumbered = self.__class__(self.ref,
                                    self.seq,
                                    seq5=seq5,
                                    name=name if name is not None else "")
        # Copy any masked positions from this region, offseting them by
        # the difference between the numbering systems.
        offset = renumbered.end5 - self.end5
        for mask_name, pos in self._masks.items():
            renumbered.add_mask(mask_name, pos + offset)
        return renumbered

    def __str__(self):
        return f"{type(self).__name__}({self.ref_reg}:{self.hyphen})"

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, Region):
            return NotImplemented
        # Compare the regions' sequences, positions, and names.
        if any([self.ref != other.ref,
                self.seq != other.seq,
                self.full != other.full,
                self.end5 != other.end5,
                self.end3 != other.end3,
                self.name != other.name]):
            return False
        # If that comparison passed, then compare their mask names.
        if sorted(self.mask_names) != sorted(other.mask_names):
            return False
        # Compare all mask values and return False if any differ.
        for name, mask in self._masks.items():
            if not np.array_equal(mask, other.get_mask(name)):
                return False
        # All checks for equality passed.
        return True

    def __ne__(self, other):
        return not self == other


def intersect(*regions: Region, name: str | None = None):
    """ Intersect one or more regions.

    Parameters
    ----------
    *regions: Region
        Regions to intersect.
    name: str | None = None
        Name for the region to return.

    Returns
    -------
    Region
        Intersection of all given regions.
    """
    if not regions:
        raise ValueError("Cannot intersect zero regions")
    # Confirm that all reference names match.
    refs = list(set(region.ref for region in regions))
    if len(refs) != 1:
        raise ValueError(f"Expected exactly one reference, but got {refs}")
    ref = refs[0]
    # Compute the 5' and 3' coordinates of the intersection.
    end5 = max(region.end5 for region in regions)
    end3 = min(region.end3 for region in regions)
    if end5 <= end3:
        # Confirm that the sequences match over the intersection.
        seqs = list(set(region.seq[end5 - region.end5:
                                   end3 - region.end5 + 1]
                        for region in regions))
        if len(seqs) != 1:
            raise ValueError(f"Expected exactly one sequence, but got {seqs}")
        seq = seqs[0]
    else:
        # The intersection among the given regions is empty.
        seq = DNA("")
        end3 = end5 - 1
    if all(region.full for region in regions):
        # If all regions are full-length, then make the intersection
        # full-length as well.
        seq5 = 1
        end5 = None
        end3 = None
    else:
        seq5 = end5
    # Create the intersection.
    intersection = Region(ref, seq, seq5=seq5, end5=end5, end3=end3, name=name)
    # Keep a position only if unmasked in all regions; equivalently,
    # mask a position if masked in any region.
    masked = regions[0].masked_int
    for region in regions[1:]:
        masked = np.union1d(region.masked_int, masked)
    # The masked positions must be limited to those in the intersection.
    masked = np.intersect1d(masked, intersection.range_int, assume_unique=True)
    if masked.size > 0:
        intersection.add_mask("intersection", masked)
    return intersection


def unite(regions: Iterable[Region],
          name: str | None = None,
          refseq: DNA | None = None):
    """ Unite one or more regions.

    Parameters
    ----------
    regions: Iterable[Region]
        Regions to unite.
    name: str | None = None
        Name for the region to return.
    refseq: DNA | None = None
        Reference sequence (optional) for filling any gaps in the union
        of the regions. If given, then it must match every region at
        the corresponding positions. If omitted, then any positions not
        covered by at least one region will be filled with N.

    Returns
    -------
    Region
        Union of all given regions.
    """
    regions = list(regions)
    if not regions:
        raise ValueError("Cannot unite 0 regions")
    # Confirm that all reference names match.
    refs = list(set(region.ref for region in regions))
    if len(refs) != 1:
        raise ValueError(f"Expected exactly one reference, but got {refs}")
    ref = refs[0]
    # Compute the 5' and 3' coordinates of the union.
    end5 = min(region.end5 for region in regions)
    end3 = max(region.end3 for region in regions)
    # Determine the coverage and sequence of the union.
    if refseq is not None:
        # Create a region from the given reference sequence.
        ref_reg = Region(ref, refseq, end5=end5, end3=end3)
        # Verify that every region matches the reference sequence.
        for reg in regions:
            # This will succeed only if the sequences match.
            intersect(reg, ref_reg)
        seq = ref_reg.seq
    else:
        # Determine a consensus sequence by overlapping the regions.
        seq_array = pd.Series(BASEN, index=np.arange(end5, end3 + 1))
        for reg in regions:
            # Verify that the region matches the sequence.
            if np.any((reg.seq.array != seq_array.loc[reg.end5: reg.end3])
                      & (seq_array.loc[reg.end5: reg.end3] != BASEN)):
                seq = DNA("".join(seq_array.loc[reg.end5: reg.end3]))
                raise ValueError(f"Sequences differ: {reg.seq} ≠ {seq}")
            # Fill in the sequence based on the region.
            seq_array.loc[reg.end5: reg.end3] = reg.seq.array
        seq = DNA("".join(seq_array))
    if all(region.full for region in regions):
        # If all regions are full-length, then also make the union
        # full-length.
        seq5 = 1
        end5 = None
        end3 = None
    else:
        seq5 = end5
    # Create the union.
    union = Region(ref, seq, seq5=seq5, end5=end5, end3=end3, name=name)
    # Keep a position only if unmasked in at least one region.
    unmasked = regions[0].unmasked_int
    for reg in regions[1:]:
        unmasked = np.union1d(reg.unmasked_int, unmasked)
    if unmasked.size < union.size:
        union.add_mask("gaps", unmasked, complement=True)
    return union


class RegionFinder(Region):
    """
    The 5' and 3' ends of a region can be given explicitly as integers,
    but if the sample is of an amplicon (i.e. generated by RT-PCR using
    site-specific primers), then it is often more convenient to enter
    the sequences of the PCR primers and have the software determine the
    coordinates. RegionFinder accepts 5' and 3' coordinates given as
    integers or primers, validates them, and stores the coordinates as
    integers, as follows:

    end5 = end5 if end5 is given, else the 3' end of the forward primer
           + (primer_gap + 1) if fwd is given, else 1
    end3 = end3 if end3 is given, else the 5' end of the reverse primer
       - (primer_gap + 1) if rev is given, else the length of refseq
    """

    def __init__(self,
                 ref: str,
                 seq: DNA, *,
                 seq5: int = 1,
                 end5: int | None = None,
                 end3: int | None = None,
                 fwd: DNA | None = None,
                 rev: DNA | None = None,
                 primer_gap: int = 0,
                 exclude_primers: bool = False,
                 **kwargs):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence.
        seq: DNA
            The full reference sequence or a part of it.
        seq5: int = 1
            Positional number to assign the 5' end of the given part of
            the reference sequence. Must be ≥ 1.
        end5: int | None = None
            Coordinate of the reference sequence at which the region's
            5' end is located.
        end3: int | None = None
            Coordinate of the reference sequence at which the region's
            3' end is located.
        fwd: DNA | None = None
            (For amplicons only) Sequence of the forward PCR primer
            that was used to generate the amplicon
        rev: DNA | None = None
            (For amplicons only) Sequence of the reverse PCR primer
            that was used to generate the amplicon (the actual sequence,
            not its reverse complement)
        primer_gap: int = 1
            (For coordinates specified by fwd/rev only) Number of
            positions 3' of the forward primer and 5' of the reverse
            primer to exclude from the region. Coordinates within 1 - 2
            nucleotides of each primer may contain DMS reactivity
            artifacts. If primer_gap = 0, then end5 and end3 are set,
            respectively, to the coordinates immediately adjacent to
            (i.e. 1 nucleotide 3' and 5' of) the 3' end of the forward
            and reverse primers.
        exclude_primers: bool = False
            Whether to exclude the primer sequences from the region.
        """
        if primer_gap < 0:
            raise ValueError(f"primer_gap must be ≥ 0, but got {primer_gap}")
        if seq is None:
            raise ValueError(f"No sequence for reference {repr(ref)}. Check "
                             "that you gave the right reference sequence file "
                             "and spelled the name of the reference correctly.")
        if end5 is None:
            # No 5' end coordinate was given.
            if fwd is not None:
                # Locate the forward primer.
                primer_site = self.locate(seq, fwd, seq5)
                if exclude_primers:
                    # Place the 5' end of the region (primer_gap + 1)
                    # positions after the 3' end of the reverse primer.
                    end5 = primer_site.pos3 + (primer_gap + 1)
                else:
                    # Place the 5' end of the region at the 5' end of
                    # the reverse primer.
                    end5 = primer_site.pos5
        if end3 is None:
            # No 3' end coordinate was given.
            if rev is not None:
                # Locate the reverse primer.
                primer_site = self.locate(seq, rev.rc, seq5)
                if exclude_primers:
                    # Place the 3' end of the region (primer_gap + 1)
                    # positions before the 5' end of the reverse primer.
                    end3 = primer_site.pos5 - (primer_gap + 1)
                else:
                    # Place the 3' end of the region at the 3' end of
                    # the reverse primer.
                    end3 = primer_site.pos3
        super().__init__(ref, seq, seq5=seq5, end5=end5, end3=end3, **kwargs)

    @staticmethod
    def locate(seq: DNA, primer: DNA, seq5: int) -> RegionTuple:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        Parameters
        ----------
        seq: DNA
            The full reference sequence or a part of it.
        primer: DNA
            Sequence of the forward PCR primer or the reverse complement
            of the reverse PCR primer
        seq5: int = 1
            Positional number to assign the 5' end of the given part of
            the reference sequence. Must be ≥ 1.

        Returns
        -------
        RegionTuple
            Named tuple of the first and last positions that the primer
            occupies in the reference sequence. Positions are 1-indexed
            and include the first and last coordinates.
        """
        matches = list(re.finditer(str(primer), str(seq)))
        if not matches:
            raise ValueError(f"Primer {primer} is not in ref {seq}")
        if len(matches) > 1:
            raise ValueError(f"Primer {primer} occurs {len(matches)} times "
                             f"in ref {seq}")
        pos5 = matches[0].start() + seq5
        pos3 = matches[0].end() + (seq5 - 1)
        return RegionTuple(pos5, pos3)


def get_coords_by_ref(coords: Iterable[tuple[str, int | DNA, int | DNA]]):
    ref_coords: dict[str, set[tuple[int | DNA, int | DNA]]] = defaultdict(set)
    for ref, end5, end3 in coords:
        coord = end5, end3
        if coord in ref_coords[ref]:
            logger.warning(f"Skipping duplicate coordinates: {coord}")
        else:
            ref_coords[ref].add(coord)
    return ref_coords


class RefRegions(object):
    """ A collection of regions, grouped by reference. """

    def __init__(self,
                 ref_seqs: Iterable[tuple[str, DNA]], *,
                 regs_file: Path | None = None,
                 coords: Iterable[tuple[str, int, int]] = (),
                 primers: Iterable[tuple[str, DNA, DNA]] = (),
                 primer_gap: int = 0,
                 exclude_primers: bool = False,
                 default_full: bool = True):
        # Get the names of the regions from the regions file, if any.
        reg_coords = dict()
        reg_primers = dict()
        if regs_file is not None:
            try:
                reg_coords, reg_primers = get_reg_coords_primers(regs_file)
                # Combines the coordinates from the regs_file and from the
                # coord parameter.
            except Exception as error:
                logger.error(error)
            else:
                coords = list(coords) + list(reg_coords)
                primers = list(primers) + list(reg_primers)
        # Group coordinates and primers by reference.
        ref_coords = get_coords_by_ref(coords)
        ref_primers = get_coords_by_ref(primers)
        # For each reference, generate regions from the coordinates.
        self._regions: dict[str, dict[tuple[int, int], Region]] = dict()
        for ref, refseq in RefSeqs(ref_seqs):
            assert ref not in self._regions
            self._regions[ref] = dict()
            for end5, end3 in ref_coords[ref]:
                # Add a region for each pair of 5' and 3' coordinates.
                reg = reg_coords.get((ref, end5, end3))
                self._add_region(ref,
                                 refseq,
                                 end5=end5,
                                 end3=end3,
                                 name=reg)
            for fwd, rev in ref_primers[ref]:
                # Add a region for each pair of fwd and rev primers.
                reg = reg_primers.get((ref, fwd, rev))
                self._add_region(ref,
                                 refseq,
                                 fwd=fwd,
                                 rev=rev,
                                 primer_gap=primer_gap,
                                 exclude_primers=exclude_primers,
                                 name=reg)
            if default_full and not self._regions[ref]:
                # If no regions were given for the reference, then add
                # a region that spans the full reference.
                self._add_region(ref, refseq)

    def _add_region(self, *args, **kwargs):
        """ Create a region and add it to the object. """
        try:
            region = RegionFinder(*args, **kwargs)
        except Exception as error:
            logger.error(error)
        else:
            # Check if the region was seen already.
            if (seen := self._regions[region.ref].get(region.coord)) is None:
                # The region was not seen already: add it.
                self._regions[region.ref][region.coord] = region
            elif seen.name == region.name:
                # The region was seen already with the same name.
                logger.warning(f"Got duplicate region: {region}")
            else:
                # The region was seen already with a different name.
                logger.error(f"Region {seen} was redefined with as {region}; "
                             f"using the first encountered: {seen}")

    def list(self, ref: str):
        """ List the regions for a given reference. """
        return list(self._regions[ref].values())

    @property
    def refs(self):
        """ Reference names. """
        return list(self._regions)

    @property
    def dict(self):
        """ List the regions for every reference. """
        return {ref: list(regions.values())
                for ref, regions in self._regions.items()}

    @property
    def regions(self):
        """ List all regions. """
        return [region for regions in self._regions.values()
                for region in regions.values()]

    @property
    def count(self):
        """ Total number of regions. """
        return sum(map(len, self._regions.values()))

    def __str__(self):
        return f"{type(self).__name__} ({self.count}): {list(self._regions)}"
