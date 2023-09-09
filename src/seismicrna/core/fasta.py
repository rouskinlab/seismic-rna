import re
from logging import getLogger
from os import linesep
from pathlib import Path
from typing import Iterable

from .path import STR_CHARS
from .seq import Seq

logger = getLogger(__name__)

# FASTA name line format.
FASTA_NAME_MARK = '>'
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


def format_fasta_seq_lines(seq: Seq, wrap: int = 0):
    """ Format a sequence in a FASTA file so that each line has at most
    `wrap` characters, or no limit if `wrap` is 0. """
    if wrap < 0:
        raise ValueError(f"Wrap cannot be negative, but got {wrap}")
    if wrap == 0 or wrap >= len(seq):
        return f"{seq}{linesep}"
    return "".join(f"{seq[start: start + wrap]}{linesep}"
                   for start in range(0, len(seq), wrap))


def format_fasta_record(name: str, seq: Seq, wrap: int = 0):
    return f"{format_fasta_name_line(name)}{format_fasta_seq_lines(seq, wrap)}"


def parse_fasta(fasta: Path, seq_type: type[Seq] | None):
    names = set()
    with open(fasta) as f:
        line = f.readline()
        # Read to the end of the file.
        while line:
            try:
                # Determine the name of the current reference.
                name = valid_fasta_seqname(line)
                if name in names:
                    raise ValueError(f"Duplicate name '{name}' in {fasta}")
                names.add(name)
            except Exception as error:
                # Determining the name failed.
                logger.error(f"Failed to parse name of reference in {fasta}: "
                             f"{error}")
                # Advance to the next line to prevent the current line
                # from being read multiple times.
                line = f.readline()
            else:
                # Advance to the next line to prevent the current line
                # from being read multiple times.
                line = f.readline()
                if seq_type is None:
                    # In name-only mode, yield just the name of the
                    # reference.
                    yield name
                else:
                    # In name-sequence mode, read the lines until the
                    # current sequence ends, then assemble the lines.
                    try:
                        segments = list()
                        while line and not line.startswith(FASTA_NAME_MARK):
                            segments.append(line.rstrip(linesep))
                            line = f.readline()
                        seq = seq_type("".join(segments))
                        logger.debug(f"Read {seq_type.__name__} reference "
                                     f"'{name}' ({len(seq)} nt) from {fasta}")
                        yield name, seq
                    except Exception as error:
                        logger.error(f"Failed to parse sequence of '{name}' "
                                     f"in {fasta}: {error}")
            # Skip to the next name line if there is one, otherwise to
            # the end of the file. Ignore blank lines.
            while line and not line.startswith(FASTA_NAME_MARK):
                line = f.readline()


def write_fasta(fasta: Path, refs: Iterable[tuple[str, Seq]],
                wrap: int = 0, overwrite: bool = False):
    """ Write an iterable of reference names and DNA sequences to a
    FASTA file. """
    logger.info(f"Began writing FASTA file: {fasta}")
    # Record the names of all the references.
    names = set()
    with open(fasta, 'w' if overwrite else 'x') as f:
        for name, seq in refs:
            try:
                # Confirm that the name is not blank.
                if not name:
                    raise ValueError(f"Got blank reference name")
                # Confirm there are no illegal characters in the name.
                if illegal := set(name) - set(FASTA_NAME_CHARS):
                    raise ValueError(f"Reference name '{name}' has illegal "
                                     f"characters: {illegal}")
                # If there are two or more references with the same name,
                # then the sequence of only the first is used.
                if name in names:
                    raise ValueError(f"Duplicate reference name: '{name}'")
                f.write(format_fasta_record(name, seq, wrap))
                logger.debug(f"Wrote reference '{name}' ({len(seq)} nt) "
                             f"to {fasta}")
                names.add(name)
            except Exception as error:
                logger.error(
                    f"Failed to write reference '{name}' to {fasta}: {error}")
    logger.info(f"Wrote {len(names)} sequences(s) to {fasta}")
