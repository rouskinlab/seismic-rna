from logging import getLogger
from pathlib import Path
from typing import Iterable

from .ct import parse_ct
from .struct import RNAStructure
from ..write import need_write, write_mode

logger = getLogger(__name__)


def from_ct(ct_path: Path, seq5: int | None = None):
    """ Yield an instance of an RNAStructure for each structure in a
    connectivity table (CT) file.

    Parameters
    ----------
    ct_path: Path
        Path of the CT file.
    seq5: int | None = None
        Positional number to assign the 5' end of the given part of the
        reference sequence. Must be ≥ 1. If omitted, then the numbering
        is taken from the last column of the CT file.

    Returns
    -------
    Generator[RNAStructure, Any, None]
        RNA secondary structures from the CT file.
    """
    titles: set[str] = set()
    for title, section, pairs in parse_ct(ct_path, seq5):
        if title in titles:
            logger.warning(f"Title {repr(title)} is repeated in {ct_path}")
        else:
            titles.add(title)
        yield RNAStructure(title=title, section=section, pairs=pairs)


def to_ct(structures: Iterable[RNAStructure],
          ct_path: Path,
          force: bool = False):
    """ Write a connectivity table (CT) file of RNA structures.

    Parameters
    ----------
    structures: Iterable[RNAStructure]
        RNA structures to write to the CT file.
    ct_path: Path
        Path of the CT file.
    force: bool = False
        Overwrite the output CT file if it already exists.
    """
    if need_write(ct_path, force):
        with open(ct_path, write_mode(force)) as f:
            for structure in structures:
                try:
                    f.write(structure.ct_text)
                except Exception as error:
                    logger.error(
                        f"Failed to write {structure} to {ct_path}: {error}"
                    )


def renumber_ct(ct_in: Path, ct_out: Path, start: int, force: bool = False):
    """ Renumber the last column of a connectivity table (CT) file.

    Parameters
    ----------
    ct_in: Path
        Path of the input CT file.
    ct_out: Path
        Path of the output CT file.
    start: int
        Number to give the first position in the renumbered CT file.
    force: bool = False
        Overwrite the output CT file if it already exists.
    """
    to_ct(from_ct(ct_in, start), ct_out, force)
