import re
from os import linesep
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

from .xna import XNA
from .. import path
from ..logs import logger
from ..write import need_write

# FASTA name line format.
FASTA_NAME_MARK = ">"
FASTA_NAME_CHARS = path.STR_CHARS
FASTA_NAME_REGEX = re.compile(f"^{FASTA_NAME_MARK}([{FASTA_NAME_CHARS}]*)")


class ReferenceNameError(ValueError):
    """Error in the name of a reference sequence."""


class BadReferenceNameError(ReferenceNameError):
    """A reference name is not valid."""


class BadReferenceNameLineError(ReferenceNameError):
    """A line that should contain a reference name is not valid."""


class DuplicateReferenceNameError(ReferenceNameError):
    """A reference name occurs more than once."""


class MissingReferenceNameError(ReferenceNameError):
    """A reference name was expected to appear but is absent."""


def extract_fasta_seqname(line: str):
    """Extract the name of a sequence from a line in FASTA format."""
    return mat.groups()[0] if (mat := FASTA_NAME_REGEX.match(line)) else None


def valid_fasta_seqname(line: str) -> str:
    """Get a valid sequence name from a line in FASTA format."""
    # Confirm that the line matches the pattern for name lines.
    if (name := extract_fasta_seqname(line)) is None:
        # If the pattern failed to match, then the line is misformatted.
        raise BadReferenceNameLineError(f"Misformatted FASTA name line {repr(line)}")
    if not name:
        # If pattern matched, then the name can still be blank.
        raise BadReferenceNameLineError(f"Blank FASTA name line {repr(line)}")
    if name != line[len(FASTA_NAME_MARK) :].rstrip():
        # If the name is not blank, then it can have illegal characters.
        raise BadReferenceNameLineError(
            f"Illegal characters in FASTA name line {repr(line)}"
        )
    # If none of the above, then the line and the name are valid.
    return name


def format_fasta_name_line(name: str):
    return f"{FASTA_NAME_MARK}{name.rstrip()}{linesep}"


def format_fasta_seq_lines(seq: XNA, wrap: int = 0):
    """Format a sequence in a FASTA file so that each line has at most
    `wrap` characters, or no limit if `wrap` is ≤ 0."""
    if not isinstance(seq, XNA):
        raise TypeError(seq)
    if 0 < wrap < len(seq):
        return "".join(
            f"{seq[start : start + wrap]}{linesep}"
            for start in range(0, len(seq), wrap)
        )
    return f"{seq}{linesep}"


def format_fasta_record(name: str, seq: XNA, wrap: int = 0):
    return f"{format_fasta_name_line(name)}{format_fasta_seq_lines(seq, wrap)}"


def parse_fasta(
    fasta: str | Path, seq_type: type[XNA] | None, only: Iterable[str] | None = None
):
    """Parse a FASTA file and yield names and (optionally) sequences.

    Parameters
    ----------
    fasta: str | Path
        Path to the FASTA file.
    seq_type: type[XNA] | None
        Type of sequence to parse (e.g. DNA or RNA); if None, only
        reference names are yielded.
    only: Iterable[str] | None
        If provided, only records whose names are in this collection are
        yielded; all others are skipped.

    Yields
    ------
    tuple[str, XNA] | str
        If `seq_type` is not None, yields (name, sequence) tuples.
        If `seq_type` is None, yields only the name of each reference.
    """
    fasta = path.sanitize(fasta, strict=True)
    path.check_file_extension(fasta, path.FastaExt)
    if seq_type is None:
        item_type = "names of sequence"
    elif issubclass(seq_type, XNA):
        item_type = f"{seq_type.__name__} sequence"
    else:
        raise ValueError(seq_type)
    with logger.debug.begin("parsing {}s in FASTA file {}", item_type, fasta):
        names = set()
        skipped = 0
        if only is not None and not isinstance(only, set):
            only = set(only)
        with open(fasta) as f:
            line = f.readline()
            # Read to the end of the file.
            while line:
                # Determine the name of the current reference.
                name = valid_fasta_seqname(line)
                if name in names:
                    raise DuplicateReferenceNameError(f"{repr(name)} in {fasta}")
                names.add(name)
                # Advance to the next line to prevent the current line from
                # being read multiple times.
                line = f.readline()
                if only is None or name in only:
                    # Yield this record if it is not the case that only some
                    # records have been selected or if the record is among
                    # those that have been selected.
                    if seq_type is not None:
                        # Read lines until the sequence ends, then assemble.
                        segments = list()
                        while line and not line.startswith(FASTA_NAME_MARK):
                            segments.append(line.rstrip(linesep))
                            line = f.readline()
                        seq = seq_type("".join(segments))
                        logger.trace("Parsed {!r} ({} nt {})", name, len(seq), item_type)
                        yield name, seq
                    else:
                        # In name-only mode, yield only the reference name.
                        logger.trace("Parsed {!r}", name)
                        yield name
                else:
                    logger.trace("Skipped {!r}", name)
                    skipped += 1
                # Skip to the next name line if there is one, otherwise to
                # the end of the file; ignore blank lines.
                while line and not line.startswith(FASTA_NAME_MARK):
                    line = f.readline()
        logger.trace("Parsed {} {}s in FASTA file {}", len(names), item_type, fasta)
        logger.trace("Skipped {} {}s in FASTA file {}", skipped, item_type, fasta)


def get_fasta_seq(fasta: str | Path, seq_type: type[XNA], name: str):
    """Get one sequence of a given name from a FASTA file."""
    try:
        _, seq = next(iter(parse_fasta(fasta, seq_type, {name})))
    except StopIteration:
        raise MissingReferenceNameError(f"{repr(name)} is not in {fasta}") from None
    return seq


def write_fasta(
    fasta: str | Path,
    refs: Iterable[tuple[str, XNA]],
    wrap: int = 0,
    force: bool = False,
):
    """Write an iterable of reference names and DNA sequences to a
    FASTA file."""
    fasta = path.sanitize(fasta, strict=False)
    path.check_file_extension(fasta, path.FastaExt)
    if need_write(fasta, force):
        with logger.debug.begin("writing FASTA file {}", fasta):
            with NamedTemporaryFile(
                "w",
                dir=fasta.parent,
                prefix=fasta.stem,
                suffix=fasta.suffix,
                delete=False,
            ) as f:
                tmp_fasta = Path(f.file.name)
            logger.debug("Created temporary FASTA file {}", tmp_fasta)
            try:
                # Write the new FASTA in a temporary file.
                with open(tmp_fasta, "w") as f:
                    # Record the names of all the references.
                    names = set()
                    for name, seq in refs:
                        # Confirm that the name is not blank.
                        if not name:
                            raise BadReferenceNameError("Blank reference name")
                        # Confirm that the name has no illegal characters.
                        if set(name) - set(FASTA_NAME_CHARS):
                            raise BadReferenceNameError(
                                f"Illegal characters in {repr(name)}"
                            )
                        if name in names:
                            raise DuplicateReferenceNameError(name)
                        f.write(format_fasta_record(name, seq, wrap))
                        logger.trace(
                            "Wrote {!r} ({} nt {} sequence)",
                            name,
                            len(seq),
                            type(seq).__name__,
                        )
                        names.add(name)
                # Release the FASTA file.
                tmp_fasta.rename(fasta)
                logger.debug("Released temporary FASTA file {} to {}", tmp_fasta, fasta)
            finally:
                # The temporary FASTA file would have been renamed already
                # if the write operation had succeeded; if not, delete it.
                try:
                    tmp_fasta.unlink()
                except FileNotFoundError:
                    pass
                else:
                    logger.debug("Deleted temporary FASTA file {}", tmp_fasta)
            logger.trace("Wrote {} sequence(s) to FASTA file {}", len(names), fasta)
    else:
        logger.trace("Skipped overwriting FASTA file {}", fasta)
