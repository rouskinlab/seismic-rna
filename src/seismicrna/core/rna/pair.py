from logging import getLogger
from typing import Generator, Iterable

import pandas as pd

from ..seq import FIELD_END5, FIELD_END3, POS_NAME, Section

logger = getLogger(__name__)

UNPAIRED = 0


def pairs_to_dict(pairs: Iterable[tuple[int, int]]):
    """ Return a dictionary that maps each position to the base to which
    it pairs and contains no key for unpaired positions. """
    # Initialize the series of pairs to 0 for every position.
    pair_dict: dict[int, int] = dict()

    def add_pair(a: int, b: int):
        """ Add a base pair at position `a` to position `b`. """
        if a < 1:
            raise ValueError(f"Position must be ≥ 1, but got {a}")
        # Find the current pairing partner at this position.
        if (b2 := pair_dict.get(a)) is None:
            # There is no pairing partner at this position: add it.
            pair_dict[a] = b
        elif b2 == b:
            logger.warning(f"Pair {a, b} was given multiple times")
        else:
            # A previous partner conflicts with the current partner.
            raise ValueError(f"Position {a} was given pairs with both "
                             f"{b2} and {b}")

    # Add all base pairs (in both directions) to the table.
    for pos1, pos2 in pairs:
        if pos1 == pos2:
            raise ValueError(f"Position {pos1} is paired with itself")
        add_pair(pos1, pos2)
        add_pair(pos2, pos1)
    return pair_dict


def dict_to_pairs(pair_dict: dict[int, int]):
    """ Tuples of the 5' and 3' position in each pair. """
    for a, b in pair_dict.items():
        if a < 1:
            raise ValueError(f"Position must be ≥ 1, but got {a}")
        if pair_dict.get(b) != a:
            raise ValueError(f"Pair {a, b} is missing its reverse {b, a}")
        if a == b:
            raise ValueError(f"Position {a} is paired with itself")
        if a < b:
            # Yield only the pairs in which "a" is 5' of "b".
            yield a, b


def pairs_to_table(pairs: Iterable[tuple[int, int]], section: Section):
    """ Series of every position in the section and the base to which it
    pairs, or 0 if it does not pair. """
    table = pd.Series(UNPAIRED, index=section.range)
    seq = str(section.seq)

    def add_pair(a: int, b: int):
        """ Add a base pair at position `a` to position `b`. """
        if not section.end5 <= a <= section.end3:
            raise ValueError(f"Position {a} is not in {section}")
        # Find the current pairing partner at this position.
        index = a, seq[a - section.end5]
        if (b2 := table[index]) == UNPAIRED:
            # There is no pairing partner at this position: add it.
            table[index] = b
        elif b2 == b:
            logger.warning(f"Pair {a, b} was given multiple times")
        else:
            # A previous partner conflicts with the current partner.
            raise ValueError(f"Position {a} was given pairs with both "
                             f"{b2} and {b}")

    # Add all base pairs (in both directions) to the table.
    for pos1, pos2 in pairs:
        if pos1 == pos2:
            raise ValueError(f"Position {pos1} is paired with itself")
        add_pair(pos1, pos2)
        add_pair(pos2, pos1)
    return table


def dict_to_table(pair_dict: dict[int, int], section: Section):
    """ Series of every position in the section and the base to which it
    pairs, or 0 if it does not pair. """
    return pairs_to_table(dict_to_pairs(pair_dict), section)


def table_to_pairs(table: pd.Series):
    """ Tuples of the 5' and 3' position in each pair. """
    pairs = table[table != UNPAIRED]
    return dict_to_pairs(dict(zip(pairs.index.get_level_values(POS_NAME),
                                  pairs,
                                  strict=True)))


def table_to_dict(table: pd.Series):
    """ Dictionary of the 5' and 3' position in each pair. """
    return pairs_to_dict(table_to_pairs(table))


def renumber_pairs(pairs: Iterable[tuple[int, int]], offset: int):
    """ Renumber pairs by offsetting each number.

    Parameters
    ----------
    pairs: Iterable[tuple[int, int]]
        Pairs to renumber.
    offset: int
        Offset by which to chage the numbering.

    Returns
    -------
    Generator[tuple[int, int], Any, None]
        Renumbered pairs, in the same order as given.
    """
    for p1, p2 in pairs:
        if min(p1, p2) < 0:
            raise ValueError(f"Positions must be ≥ 1, but got {p1, p2}")
        r1, r2 = p1 + offset, p2 + offset
        if min(r1, r2) < 0:
            raise ValueError(f"Positions must be ≥ 1, but got {r1, r2}")
        yield r1, r2


def find_enclosing_pairs(table: pd.Series):
    """ Find the base pair that encloses each position. """
    enclosing = pd.DataFrame(UNPAIRED, table.index, [FIELD_END5, FIELD_END3])
    stack = list()
    for (position, _), partner in table.items():
        if partner != UNPAIRED:
            if partner == position:
                raise ValueError(f"Base {position} is paired with itself")
            if partner > position:
                # Create a new pair.
                pair = position, partner
                stack.append(pair)
                enclosing.loc[position] = pair
            else:
                # Remove the last pair.
                pair = stack.pop()
                if pair != (partner, position):
                    raise ValueError(f"Pairs {partner, position} and {pair} "
                                     f"are not nested")
                enclosing.loc[position] = pair
        elif stack:
            enclosing.loc[position] = stack[-1]
    if stack:
        raise ValueError(f"Remaining base pairs: {stack}")
    return enclosing


def map_nested(pairs: Iterable[tuple[int, int]]):
    """ Map each pair to the pair in which it is nested. """
    pair_dict = pairs_to_dict(pairs)
    min_pos = min(pair_dict)
    max_pos = max(pair_dict)
    root_pos = set()

    def find_nested(pair: tuple[int, int]):
        """ Find the pair in which one pair is nested. """
        p5, p3 = pair
        if p5 >= p3:
            raise ValueError(f"Invalid pair: {p5, p3}")
        # Find the nearest pair on the 5' side.
        while (((c5 := pair_dict.get(p5 := p5 - 1)) is None or c5 < p5)
               and p5 not in root_pos
               and p5 >= min_pos):
            pass
        # Find the nearest pair on the 3' side.
        while (((c3 := pair_dict.get(p3 := p3 + 1)) is None or c3 > p3)
               and p3 not in root_pos
               and p3 <= max_pos):
            pass
        if p5 < min_pos and p3 > max_pos:
            # No other pair contains this pair.
            root_pos.update(pair)
            return None
        if p5 != c3 or p3 != c5:
            raise ValueError(f"{pair} and {c3} are not nested")
        if pair_dict.get(p5) != p3:
            raise ValueError(f"{pair} and")
        if p3 <= max_pos:
            raise ValueError("Found a pseudoknot: {}")
        return None

    return {(a, b): find_nested((a, b)) for a, b in pair_dict.items() if a < b}


def _find_root_pairs_nested(pair_dict: dict[int, int]):
    """ Return all pairs that are not contained in any other pair, while
    assuming that all pairs are nested (i.e. no pseudoknots). """
    if pair_dict:
        max_pos = max(pair_dict)
        pos = min(pair_dict)
        while pos <= max_pos:
            if partner := pair_dict.get(pos, UNPAIRED):
                if partner == pos:
                    raise ValueError(f"Base {pos} is paired with itself")
                if partner < pos:
                    raise ValueError(f"Pair {pos, partner} is not nested")
                yield pos, partner
                pos = partner
            pos += 1


def find_root_pairs(pairs: Iterable[tuple[int, int]],
                    assume_nested: bool = False):
    """ Return all pairs which are not contained any other pair. """
    pair_dict = pairs_to_dict(pairs)
    return list(_find_root_pairs_nested(pair_dict))

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
