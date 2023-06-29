import re
from typing import Sequence

import pandas as pd

from ..core.sect import (index_to_pos, index_to_seq, seq_pos_to_index,
                         INDEX_NAMES)
from ..core.seq import BASES, DNA


def format_seq_pos(seq: DNA, positions: Sequence[int], start: int):
    """
    Convert sequence and positions to indexes, where each index is a
    string of the base followed by the position.

    Parameters
    ----------
    seq: DNA
        Sequence of DNA
    positions: Sequence[int]
        Positions of the sequence from which to build the index. Every
        position must be an integer ≥ `start`.
    start: int
        Numerical position to assign to the first base in the sequence.
        Must be a positive integer.

    Returns
    -------
    pd.Index
        Index of the same length as positions where each element is a
        string of {base}{position}.
    """
    return pd.Index([f"{base}{pos}" for pos, base
                     in seq_pos_to_index(seq, positions, start)])


def _parse_seq_pos(index: Sequence[str]):
    """
    Convert indexes to a sequence and positions, where each index is a
    string of the base followed by the position.

    Parameters
    ----------
    index: Sequence[str]
        Index where each element is a string of {base}{position}.

    Returns
    -------
    pd.MultiIndex
        MultiIndex of the same length as positions where each index is a
        tuple of (position, base).
    """
    if len(index) == 0:
        raise ValueError("Got empty index")
    # Regex pattern "^([ACGT])([0-9]+)$" finds the base and position
    # in each index name.
    pattern = re.compile(f"^([{BASES.decode()}])([0-9]+)$")
    try:
        # Match each index name using the pattern.
        matches = list(map(pattern.match, index))
    except TypeError:
        raise ValueError(f"Invalid indexes: {index}")
    try:
        # Obtain the two groups (base and position) from each match
        # and unzip them into two tuples.
        bases, pos_strs = zip(*map(re.Match.groups, matches))
    except TypeError:
        # TypeError is raised if any match is None, which happens if
        # a column fails to match the pattern.
        invalid_idx = [idx for idx, match in zip(index, matches, strict=True)
                       if match is None]
        raise ValueError(f"Invalid indexes: {invalid_idx}")
    # Assemble the positions and bases into a NON-VALIDATED MultiIndex.
    return pd.MultiIndex.from_arrays([list(map(int, pos_strs)), list(bases)],
                                     names=INDEX_NAMES)


def parse_seq(index: Sequence[str]):
    """ Convert indexes to a sequence, where each index is a string of
    the base followed by the position. """
    return index_to_seq(_parse_seq_pos(index))


def parse_pos(index: Sequence[str]):
    """ Convert indexes to an array of positions, where each index is a
    string of the base followed by the position. """
    return index_to_pos(_parse_seq_pos(index))
