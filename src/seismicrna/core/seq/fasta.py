import re
from logging import getLogger
from os import linesep
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

from .xna import XNA
from ..path import STR_CHARS
from ..write import need_write

logger = getLogger(__name__)

# FASTA name line format.
FASTA_NAME_MARK = ">"
FASTA_NAME_CHARS = STR_CHARS
FASTA_NAME_REGEX = re.compile(f"^{FASTA_NAME_MARK}([{FASTA_NAME_CHARS}]*)")


def extract_fasta_seqname(line: str):
    """ Extract the name of a sequence from a line in FASTA format. """
    return mat.groups()[0] if (mat := FASTA_NAME_REGEX.match(line)) else None


def valid_fasta_seqname(line: str) -> str:
    """ Get a valid sequence name from a line in FASTA format. """
    # Confirm that the line matches the pattern for name lines.
    if (name := extract_fasta_seqname(line)) is None:
        # If the pattern failed to match, then the line is misformatted.
        raise ValueError(f"FASTA name line {repr(line)} is misformatted")
    if not name:
        # If pattern matched, then the name can still be blank.
        raise ValueError(f"FASTA name line {repr(line)} has a blank name")
    if name != line[len(FASTA_NAME_MARK):].rstrip():
        # If the name is not blank, then it can have illegal characters.
        raise ValueError(f"FASTA name line {repr(line)} has illegal characters")
    # If none of the above, then the line and the name are valid.
    return name


def format_fasta_name_line(name: str):
    return f"{FASTA_NAME_MARK}{name.rstrip()}{linesep}"


def format_fasta_seq_lines(seq: XNA, wrap: int = 0):
    """ Format a sequence in a FASTA file so that each line has at most
    `wrap` characters, or no limit if `wrap` is 0. """
    if wrap < 0:
        raise ValueError(f"Wrap must be ≥ 0, but got {wrap}")
    if wrap == 0 or wrap >= len(seq):
        return f"{seq}{linesep}"
    return "".join(f"{seq[start: start + wrap]}{linesep}"
                   for start in range(0, len(seq), wrap))


def format_fasta_record(name: str, seq: XNA, wrap: int = 0):
    return f"{format_fasta_name_line(name)}{format_fasta_seq_lines(seq, wrap)}"


def parse_fasta(fasta: Path,
                seq_type: type[XNA] | None,
                only: Iterable[str] | None = None):
    names = set()
    if only is not None and not isinstance(only, set):
        only = set(only)
    with open(fasta) as f:
        line = f.readline()
        # Read to the end of the file.
        while line:
            # Determine the name of the current reference.
            name = valid_fasta_seqname(line)
            if name in names:
                raise ValueError(f"Duplicate name {repr(name)} in {fasta}")
            names.add(name)
            # Advance to the next line to prevent the current line from
            # being read multiple times.
            line = f.readline()
            if only is None or name in only:
                # Yield this record if it is not the case that only some
                # records have been selected or if the record is among
                # those that have been selected.
                if seq_type is None:
                    # In name-only mode, yield just the name of the
                    # reference.
                    yield name
                else:
                    # Otherwise, read the lines until the sequence ends,
                    # then assemble the lines.
                    segments = list()
                    while line and not line.startswith(FASTA_NAME_MARK):
                        segments.append(line.rstrip(linesep))
                        line = f.readline()
                    seq = seq_type("".join(segments))
                    logger.debug(f"Read {seq_type.__name__} {repr(name)} "
                                 f"({len(seq)} nt) from {fasta}")
                    yield name, seq
            # Skip to the next name line if there is one, otherwise to
            # the end of the file; ignore blank lines.
            while line and not line.startswith(FASTA_NAME_MARK):
                line = f.readline()
    logger.info(f"Read {len(names)} sequences(s) from {fasta}")


def get_fasta_seq(fasta: Path, seq_type: type[XNA], name: str):
    """ Get one sequence of a given name from a FASTA file. """
    if not isinstance(seq_type, type) and issubclass(seq_type, XNA):
        raise TypeError("Expected seq_type to be subclass of Seq, but got "
                        f"{repr(seq_type)}")
    try:
        _, seq = next(iter(parse_fasta(fasta, seq_type, (name,))))
    except StopIteration:
        raise ValueError(f"Sequence {repr(name)} is not in {fasta}") from None
    return seq


def write_fasta(fasta: Path,
                refs: Iterable[tuple[str, XNA]],
                wrap: int = 0,
                force: bool = False):
    """ Write an iterable of reference names and DNA sequences to a
    FASTA file. """
    if need_write(fasta, force):
        with NamedTemporaryFile("w",
                                dir=fasta.parent,
                                prefix=fasta.stem,
                                suffix=fasta.suffix,
                                delete=False) as f:
            tmp_fasta = Path(f.file.name)
        try:
            # Write the new FASTA in a temporary file.
            with open(tmp_fasta, "w") as f:
                # Record the names of all the references.
                names = set()
                for name, seq in refs:
                    # Confirm that the name is not blank.
                    if not name:
                        raise ValueError(f"Got blank reference name")
                    # Confirm that the name has no illegal characters.
                    if illegal := (set(name) - set(FASTA_NAME_CHARS)):
                        raise ValueError(f"Reference name {repr(name)} has "
                                         f"illegal characters: {illegal}")
                    if name in names:
                        raise ValueError(f"Duplicate reference: {repr(name)}")
                    f.write(format_fasta_record(name, seq, wrap))
                    logger.debug(f"Wrote reference {repr(name)} "
                                 f"({len(seq)} nt) to {tmp_fasta}")
                    names.add(name)
            # Release the FASTA file.
            tmp_fasta.rename(fasta)
        finally:
            # The temporary FASTA file would have been renamed already
            # if the write operation had succeeded; if not, delete it.
            tmp_fasta.unlink(missing_ok=True)
        logger.info(f"Wrote {len(names)} sequences(s) to {fasta}")

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
