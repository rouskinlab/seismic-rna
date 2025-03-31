from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO

from .. import path
from ..seq import RNA, Region

DB_NAME_MARK = ">"
UNPAIRED_MARK = "."
PAIRED_MARKS = {")": "(",
                ">": "<",
                "]": "[",
                "}": "{"}


def _parse_db_file_header(header_line: str):
    if not header_line.startswith(DB_NAME_MARK):
        raise ValueError(f"Header line {repr(header_line)} does not start with "
                         f"{repr(DB_NAME_MARK)}")
    return header_line[len(DB_NAME_MARK):].rstrip("\n")


def _parse_db_file_next_record(db_file: TextIO, seq: RNA | None):
    if seq is None:
        seq = RNA.from_any_seq(db_file.readline().rstrip("\n"))
    db_string = db_file.readline().rstrip("\n")
    if len(db_string) != len(seq):
        raise ValueError(
            f"Lengths of structure {repr(db_string)} ({len(db_string)}) "
            f"and sequence {seq} ({len(seq)}) do not match"
        )
    return seq, db_string


def parse_db_file_as_strings(db_path: str | Path):
    """ Return the sequence and dot-bracket strings from a dot-bracket
    file. """
    seq = None
    db_strings = dict()
    with open(db_path) as file:
        while header_line := file.readline():
            # Get the header of the current structure.
            header = _parse_db_file_header(header_line)
            # Determine the sequence and base pairs.
            seq, db_string = _parse_db_file_next_record(file, seq)
            db_strings[header] = db_string
    if seq is None:
        raise ValueError(f"File {db_path} contained no sequence")
    return seq, db_strings


def parse_db_string(db_string: str, seq5: int = 1):
    """ Parse a dot-bracket string into a list of base pairs. """
    stacks: dict[str, list[int]] = defaultdict(list)
    pairs = list()
    opening_marks = "".join(PAIRED_MARKS.values())
    for pos, mark in enumerate(db_string, start=seq5):
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


def parse_db_file_as_pairs(db_path: str | Path, seq5: int = 1):
    """ Yield the title, region, and base pairs for each structure in a
    dot-bracket (DB) file.

    Parameters
    ----------
    db_path: str | Path
        Path of the DB file.
    seq5: int = 1
        Number to give the 5' position of the sequence.

    Returns
    -------
    Generator[tuple[str, Region, list[tuple[int, int]]], Any, None]
    """
    # Determine the reference and region names from the path.
    fields = path.parse(db_path, path.DB_FILE_LAST_SEGS)
    ref = fields[path.REF]
    reg = fields[path.REG]
    # Parse each structure in the CT file.
    seq = None
    with open(db_path) as db_file:
        while header_line := db_file.readline():
            # Get the header of the current structure.
            title = _parse_db_file_header(header_line)
            # Determine the sequence and base pairs.
            seq, db_string = _parse_db_file_next_record(db_file, seq)
            pairs = parse_db_string(db_string, seq5)
            # Make a region from the sequence.
            region = Region(ref, seq.rt(), seq5=seq5, name=reg)
            # Yield the title, region, and base pairs.
            yield title, region, pairs


def format_db_string(pairs: Iterable[tuple[int, int]],
                     length: int,
                     seq5: int = 1):
    """ Create a dot-bracket string from a list of base pairs. """
    db = [UNPAIRED_MARK] * length
    for pos5, pos3 in sorted(pairs):
        i = pos5 - seq5
        j = pos3 - seq5
        if i < 0:
            raise ValueError(f"5' partner must be ≥ {seq5}, but got {pos5}")
        if i >= j:
            raise ValueError("5' partner must be less than 3' partner, "
                             f"but got {pos5} ≥ {pos3}")
        if j >= length:
            raise ValueError(f"3' partner must be ≤ {length + seq5 - 1}, "
                             f"but got {pos3}")
        # Determine which mark to use (it must not be used already).
        used_marks = set(db[i: j])
        for cmark, omark in PAIRED_MARKS.items():
            if omark not in used_marks and cmark not in used_marks:
                if db[i] != UNPAIRED_MARK:
                    raise ValueError(f"Got >1 base pair for position {pos5}")
                if db[j] != UNPAIRED_MARK:
                    raise ValueError(f"Got >1 base pair for position {pos3}")
                db[i] = omark
                db[j] = cmark
                break
        else:
            raise ValueError(f"Cannot write base-pair {pos5, pos3} because all "
                             f"marks are already used in {''.join(db[i: j])}")
    return "".join(db)
