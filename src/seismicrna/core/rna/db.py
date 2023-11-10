from logging import getLogger
from pathlib import Path
from typing import TextIO

from ..seq import RNA

logger = getLogger(__name__)


DB_NAME_MARK = ">"


def _parse_db_header(header_line: str):
    if not header_line.startswith(DB_NAME_MARK):
        raise ValueError(f"Header line {repr(header_line)} does not start with "
                         f"{repr(DB_NAME_MARK)}")
    return header_line[len(DB_NAME_MARK):]


def _parse_db_structure(db_file: TextIO, seq: RNA | None):
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
            seq, struct = _parse_db_structure(file, seq)
            structs[header] = struct
    if seq is None:
        raise ValueError(f"File {db_path} contained no sequence")
    return seq, structs

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
