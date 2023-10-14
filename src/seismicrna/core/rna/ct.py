from logging import getLogger
from pathlib import Path
from typing import TextIO

from ..seq import RNA

logger = getLogger(__name__)


NUM_FIELDS = 6


def _parse_int(text: str, name: str, zero: bool = False) -> int:
    """ Try to parse the text into an integer/positive integer. """
    try:
        value = int(text)
    except ValueError:
        value = None
    if value is None or value < 0 or (value == 0 and not zero):
        s = "non-negative" if zero else "positive"
        raise ValueError(f"{name} must be a {s} integer, but got {repr(text)}")
    return value


def _parse_ct_header_line(line: str):
    """ Get the title and sequence length from a CT header line. """
    content = line.strip()
    if not content:
        raise ValueError(f"Got all-whitespace CT header line")
    # Determine the length of the sequence, which is the first item
    # in the header line.
    length_str = content.split()[0]
    length = _parse_int(length_str, "sequence length")
    # Determine the title, which is the part of the line following
    # the sequence length.
    title = content[len(length_str):].lstrip()
    return title, length


def _parse_ct_body_line(line: str, first: bool, last: bool):
    """ Get the position and pairing data from a CT body line. """
    content = line.strip()
    if not content:
        raise ValueError(f"Got all-whitespace CT body line")
    # Parse each whitespace-delimited field of the CT file.
    fields = content.split()
    if len(fields) != NUM_FIELDS:
        raise ValueError(f"CT body line needs {NUM_FIELDS} fields,"
                         f"but got {len(fields)}: '{content}'")
    curr_idx = _parse_int(fields[0], "current index")
    prev_idx = _parse_int(fields[2], "previous index", zero=first)
    if first:
        if prev_idx != 0:
            raise ValueError(f"Expected previous index of first line "
                             f"to be 0, but got {prev_idx}")
    elif prev_idx != curr_idx - 1:
        raise ValueError(f"Previous index ({prev_idx}) does not precede "
                         f"current index ({curr_idx})")
    next_idx = _parse_int(fields[3], "next index", zero=last)
    if last:
        if next_idx != 0:
            raise ValueError(f"Expected next index of last line to be 0, "
                             f"but got {next_idx}")
    elif next_idx != curr_idx + 1:
        raise ValueError(f"Next index ({next_idx}) does not succeed "
                         f"current index ({curr_idx})")
    partner = _parse_int(fields[4], "partner index", zero=True)
    position = _parse_int(fields[5], "natural position")
    base = fields[1]
    return curr_idx, base, partner, position


def _parse_ct_structure(ct_file: TextIO,
                        length: int,
                        index_offset: int | None = None):
    """ Return the sequence and pairs for the current structure. """
    # Initialize the bases, pairs, and position numbers.
    bases: list[str] = list()
    pairs: dict[int, int] = dict()
    reverse_pairs: dict[int, int] = dict()
    unpaired: set[int] = set()
    positions_map: dict[int, int] = dict()
    positions_set: set[int] = set()
    # Read the next length lines from the file: one per position in
    # the structure.
    indexes = list(range(1, length + 1))
    for index, line in zip(indexes, ct_file, strict=False):
        is_first = index == 1
        is_last = index == length
        # Parse the current line in the CT file.
        ct_index, base, partner, position = _parse_ct_body_line(line,
                                                                is_first,
                                                                is_last)
        if ct_index != index:
            raise ValueError(f"Index ({ct_index}) and line # ({index}) differ")
        # Add the current base to the growing sequence.
        bases.append(base)
        # Check that the position has not been encountered, and
        # record it.
        if position in positions_set:
            raise ValueError(f"Position {position} appears twice")
        positions_set.add(position)
        # If an index offset was given, then map the index to the
        # index plus the offset. Otherwise, map the index to the
        # natural position given in the CT file.
        positions_map[index] = (position if index_offset is None
                                else index + index_offset)
        # Check if the current base is unpaired.
        if partner == 0:
            # If the current base is unpaired, then it should not
            # appear as the partner of any other base.
            if (partner0 := reverse_pairs.get(index)) is not None:
                raise ValueError(f"Base {index} is not paired, but base "
                                 f"{partner0} claims to pair with it")
            # Mark this base as unpaired and skip to the next.
            unpaired.add(index)
            continue
        # Validate the pairing partner.
        if partner > length:
            raise ValueError(f"Index of partner ({partner}) exceeds length "
                             f"of sequence ({length})")
        if partner in unpaired:
            raise ValueError(f"Base {index} claims to pair with {partner}, "
                             f"but base {partner} claims to be unpaired")
        # Handle the pair depending on whether the current base
        # comes before or after its partner.
        if index < partner:
            # This is the first time this pair appears: record it.
            pairs[index] = partner
            # Check if another base has already claimed to pair with
            # the partner of the current base.
            if (other := reverse_pairs.get(partner)) is not None:
                raise ValueError(f"Bases {other} and {index} both claim to "
                                 f"pair with base {partner}")
            reverse_pairs[partner] = index
        elif index > partner:
            # This pair should have appeared already. Confirm that
            # it does, and do not add it again.
            if (check := pairs.get(partner)) != index:
                raise ValueError(f"Base {index} claims to pair with "
                                 f"{partner}, but base {partner} claims "
                                 f"to pair with {check}")
            if (check := reverse_pairs.get(index)) != partner:
                raise ValueError(f"Base {index} claims to pair with "
                                 f"{partner}, but base {check} claims "
                                 f"to pair with {index}")
        else:
            # Index and partner are equal.
            raise ValueError(f"Base {index} claims to pair with itself")
    # Confirm that the expected number of lines were read.
    if len(bases) != length:
        raise ValueError(f"Got {len(bases)} lines, but expected {length}")
    # Confirm that all paired bases appear in both maps of pairs.
    expect_paired = set(indexes) - unpaired
    if (is_paired := set(pairs).union(set(reverse_pairs))) != expect_paired:
        raise ValueError(f"Paired bases {is_paired} do not match bases "
                         f"expected to be paired {expect_paired}")
    # Assemble the list of bases into an RNA sequence.
    seq = RNA("".join(bases))
    # Map all the indexes to their corresponding positions.
    pairs_list = [(positions_map[index1], positions_map[index2])
                  for index1, index2 in pairs.items()]
    return seq, pairs_list


def parse_ct(ct_path: Path, start: int | None = None):
    """ Yield the sequence and list of base pairs for each structure
    from a connectivity table (CT) file. """
    # If a starting position was given, then offset the index by that
    # position minus 1 (because CT files are 1-indexed).
    offset = start if start is None else start - 1
    # Parse each structure in the CT file.
    with open(ct_path) as file:
        while header_line := file.readline():
            # Get the title and length of the current structure.
            curr_title, curr_length = _parse_ct_header_line(header_line)
            # Determine the sequence and base pairs.
            curr_seq, curr_pairs = _parse_ct_structure(file, curr_length, offset)
            # Yield the title, sequence, and base pairs.
            yield curr_title, curr_seq, curr_pairs

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
