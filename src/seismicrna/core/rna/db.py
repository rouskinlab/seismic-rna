from collections import defaultdict
from logging import getLogger
from pathlib import Path
from typing import TextIO

from .. import path
from ..seq import RNA, Section

logger = getLogger(__name__)

DB_NAME_MARK = ">"
UNPAIRED_MARK = "."
PAIRED_MARKS = {")": "(",
                "]": "[",
                "}": "{",
                ">": "<"}


def _parse_db_header(header_line: str):
    if not header_line.startswith(DB_NAME_MARK):
        raise ValueError(f"Header line {repr(header_line)} does not start with "
                         f"{repr(DB_NAME_MARK)}")
    return header_line[len(DB_NAME_MARK):].rstrip("\n")


def _parse_db_string(db_file: TextIO, seq: RNA | None):
    if seq is None:
        seq = RNA(db_file.readline().rstrip("\n"))
    struct = db_file.readline().rstrip("\n")
    if len(struct) != len(seq):
        raise ValueError(f"Lengths of structure {repr(struct)} ({len(struct)}) "
                         f"and sequence {seq} ({len(seq)}) do not match")
    return seq, struct


def parse_db_strings(db_path: Path):
    """ Return the sequence and structures from a dot-bracket file. """
    seq = None
    structs = dict()
    with open(db_path) as file:
        while header_line := file.readline():
            # Get the header of the current structure.
            header = _parse_db_header(header_line)
            # Determine the sequence and base pairs.
            seq, struct = _parse_db_string(file, seq)
            structs[header] = struct
    if seq is None:
        raise ValueError(f"File {db_path} contained no sequence")
    return seq, structs


def _parse_db_structure(struct: str, seq5: int):
    stacks: dict[str, list[int]] = defaultdict(list)
    pairs = list()
    opening_marks = "".join(PAIRED_MARKS.values())
    for pos, mark in enumerate(struct, start=seq5):
        if mark != UNPAIRED_MARK:
            if mark in opening_marks:
                stacks[mark].append(pos)
            elif omark := PAIRED_MARKS.get(mark):
                try:
                    pairs.append((stacks[omark].pop(), pos))
                except IndexError:
                    raise ValueError(
                        f"Position {pos} has an unmatched {repr(mark)}"
                    ) from None
            else:
                raise ValueError(
                    f"Position {pos} has an invalid mark: {repr(mark)}"
                )
    for mark, positions in stacks.items():
        if positions:
            raise ValueError(
                f"Position {positions[0]} has an unmatched {repr(mark)}"
            )
    return sorted(pairs)


def _parse_db_record(db_file: TextIO, seq: RNA | None, seq5: int):
    seq, struct = _parse_db_string(db_file, seq)
    return seq, _parse_db_structure(struct, seq5)


def parse_db(db_path: Path, seq5: int = 1):
    """ Yield the title, section, and base pairs for each structure in a
    dot-bracket (DB) file.

    Parameters
    ----------
    db_path: Path
        Path of the DB file.
    seq5: int = 1
        Number to give the 5' position of the sequence.

    Returns
    -------
    Generator[tuple[str, Section, list[tuple[int, int]]], Any, None]
    """
    # Determine the reference and section names from the path.
    fields = path.parse(db_path,
                        path.RefSeg,
                        path.SectSeg,
                        path.DotBracketSeg)
    ref = fields[path.REF]
    sect = fields[path.SECT]
    # Parse each structure in the CT file.
    seq = None
    with open(db_path) as file:
        while header_line := file.readline():
            # Get the header of the current structure.
            title = _parse_db_header(header_line)
            # Determine the sequence and base pairs.
            seq, pairs = _parse_db_record(file, seq, seq5)
            # Make a section from the sequence.
            section = Section(ref, seq.rt(), seq5=seq5, name=sect)
            # Yield the title, section, and base pairs.
            yield title, section, pairs

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
