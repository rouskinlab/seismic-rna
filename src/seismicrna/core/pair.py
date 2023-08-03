from logging import getLogger
from pathlib import Path
from typing import BinaryIO, Iterable

import pandas as pd

from .sect import Section, POS_NAME
from .seq import RNA

logger = getLogger(__name__)


def pairs_to_partners(pairs: Iterable[tuple[int, int]], section: Section):
    """ Return a Series of every position in the section and the base to
    which it pairs, or 0 if it does not pair. """
    # Initialize the series of pairs to 0 for every position.
    partners = pd.Series(0, index=pd.Index(section.range_int, name=POS_NAME))

    def add_pair(at: int, to: int):
        """ Add a base pair at position `at` to position `to`. """
        # Find the current pairing partner at this position.
        try:
            to2 = partners.loc[at]
        except KeyError:
            raise ValueError(f"Position {at} is not in {section}")
        if to2 == 0:
            # There is no pairing partner at this position: add it.
            partners.loc[at] = to
        elif to2 == to:
            # A previous partner matches the current partner.
            logger.warning(f"Pair {at}-{to} was given twice")
        else:
            # A previous partner conflicts with the current partner.
            raise ValueError(f"Position {at} was given pairs with both "
                             f"{to2} and {to}")

    # Add all base pairs (in both directions) to the table.
    for pos1, pos2 in pairs:
        add_pair(pos1, pos2)
        add_pair(pos2, pos1)
    return partners


def parse_ct_pairs(ct_path: Path, start: int | None = None):
    """ Yield the sequence and list of base pairs for each structure
    from a connectivity table file. """
    n_cols = 6

    def parse_int(text: bytes, name: str, zero: bool = False) -> int:
        """ Try to parse the text into an integer/positive integer. """
        try:
            value = int(text)
        except ValueError:
            value = None
        if value is None or value < 0 or (value == 0 and not zero):
            kind = "non-negative" if zero else "positive"
            raise ValueError(f"{name.capitalize()} must be a {kind} integer, "
                             f"but got '{text.decode()}' in {ct_path}")
        return value

    def parse_ct_header_line(line: bytes):
        """ Get the title and sequence length from a CT header line. """
        content = line.strip()
        if not content:
            raise ValueError(f"Got all-whitespace CT header line in {ct_path}")
        # Determine the length of the sequence, which is the first item
        # in the header line.
        length_str = content.split()[0]
        length = parse_int(length_str, "sequence length")
        # Determine the title, which is the part of the line following
        # the sequence length.
        title = content[len(length_str):].lstrip().decode()
        return title, length

    def parse_ct_body_line(line: bytes, first: bool, last: bool):
        """ Get the position and pairing data from a CT body line. """
        content = line.strip()
        if not content:
            raise ValueError(f"Got all-whitespace CT body line in {ct_path}")
        # Parse each whitespace-delimited field of the CT file.
        fields = content.split()
        if len(fields) != n_cols:
            raise ValueError(f"CT body line needs {n_cols} fields,"
                             f"but got {len(fields)}: '{content}'")
        curr_idx = parse_int(fields[0], "current index")
        prev_idx = parse_int(fields[2], "previous index", zero=first)
        if first:
            if prev_idx != 0:
                raise ValueError(f"Expected previous index of first line "
                                 f"to be 0, but got {prev_idx}")
        elif prev_idx != curr_idx - 1:
            raise ValueError(f"Previous index ({prev_idx}) does not precede "
                             f"current index ({curr_idx})")
        next_idx = parse_int(fields[3], "next index", zero=last)
        if last:
            if next_idx != 0:
                raise ValueError(f"Expected next index of last line to be 0, "
                                 f"but got {next_idx}")
        elif next_idx != curr_idx + 1:
            raise ValueError(f"Next index ({next_idx}) does not succeed "
                             f"current index ({curr_idx})")
        partner = parse_int(fields[4], "partner index", zero=True)
        position = parse_int(fields[5], "natural position")
        base = fields[1]
        return curr_idx, base, partner, position

    def parse_ct_structure(ct_file: BinaryIO, length: int,
                           index_offset: int | None = None):
        """ Return the sequence and pairs for the current structure. """
        # Initialize the bases, pairs, and position numbers.
        bases: list[bytes] = list()
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
            ct_index, base, partner, position = parse_ct_body_line(line,
                                                                   is_first,
                                                                   is_last)
            if ct_index != index:
                raise ValueError(f"Index ({ct_index}) does not match line "
                                 f"number ({index}) in {ct_path}")
            # Add the current base to the growing sequence.
            bases.append(base)
            # Check that the position has not been encountered, and
            # record it.
            if position in positions_set:
                raise ValueError(f"Position {position} appears twice "
                                 f"in {ct_path}")
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
                                 f"of sequence ({length}) in {ct_path}")
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
                                     f"pair with base {partner} in {ct_path}")
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
            raise ValueError(f"Got {len(bases)} lines in {ct_path}, "
                             f"but expected {length}")
        # Confirm that all paired bases appear in both maps of pairs.
        expect_paired = set(indexes) - unpaired
        if (is_paired := set(pairs).union(set(reverse_pairs))) != expect_paired:
            raise ValueError(f"Paired bases {is_paired} do not match bases "
                             f"expected to be paired {expect_paired}")
        # Assemble the list of bases into an RNA sequence.
        seq = RNA(b"".join(bases))
        # Map all the indexes to their corresponding positions.
        pairs_list = [(positions_map[index1], positions_map[index2])
                      for index1, index2 in pairs.items()]
        return seq, pairs_list

    # If a starting position was given, then offset the index by that
    # position minus 1 (because CT files are 1-indexed).
    offset = start if start is None else start - 1
    # Parse each structure in the CT file.
    with open(ct_path, "rb") as file:
        while header_line := file.readline():
            # Get the title and length of the current structure.
            curr_title, curr_length = parse_ct_header_line(header_line)
            # Determine the sequence and base pairs.
            curr_seq, curr_pairs = parse_ct_structure(file, curr_length, offset)
            # Yield the title, sequence, and base pairs.
            yield curr_title, curr_seq, curr_pairs
