from logging import getLogger
from typing import Iterable

import pandas as pd

from ..seq import POS_INDEX, Section

logger = getLogger(__name__)


def pairs_to_dict(pairs: Iterable[tuple[int, int]]):
    """ Return a dictionary that maps each position to the base to which
    it pairs and contains no key for unpaired positions. """
    # Initialize the series of pairs to 0 for every position.
    pair_dict: dict[int, int] = dict()

    def add_pair(a: int, b: int):
        """ Add a base pair at position `at` to position `to`. """
        if a < POS_INDEX:
            raise ValueError(f"Position must be ≥ {POS_INDEX}, but got {a}")
        # Find the current pairing partner at this position.
        if (to2 := pair_dict.get(a)) is None:
            # There is no pairing partner at this position: add it.
            pair_dict[a] = b
        elif to2 == b:
            logger.warning(f"Pair {a, b} was given multiple times")
        else:
            # A previous partner conflicts with the current partner.
            raise ValueError(f"Position {a} was given pairs with both "
                             f"{to2} and {b}")

    # Add all base pairs (in both directions) to the table.
    for pos1, pos2 in pairs:
        if pos1 == pos2:
            raise ValueError(f"Base {pos1} is paired with itself")
        add_pair(pos1, pos2)
        add_pair(pos2, pos1)
    return pair_dict


def dict_to_pairs(pair_dict: dict[int, int]):
    """ Tuples of the 5' and 3' position in each pair. """
    for a, b in pair_dict.items():
        if pair_dict.get(b) != a:
            raise ValueError(f"Pair {a, b} is missing its reverse {b, a}")
        if a == b:
            raise ValueError(f"Base {a} is paired with itself")
        if a < b:
            # Yield only the pairs in which at is 5' of to.
            yield a, b


def pairs_to_table(pairs: Iterable[tuple[int, int]], section: Section):
    """ Return a Series of every position in the section and the base to
    which it pairs, or 0 if it does not pair. """
    table = pd.Series(0, index=section.range_int)
    for a, b in pairs_to_dict(pairs).items():
        try:
            table.loc[a] = b
        except KeyError:
            raise ValueError(f"Position {a} is not in {section}")
    return table


def table_to_pairs(table: pd.Series):
    """ Tuples of the 5' and 3' position in each pair. """
    return dict_to_pairs(table[table != 0].to_dict())


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
            if partner := pair_dict.get(pos, 0):
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
