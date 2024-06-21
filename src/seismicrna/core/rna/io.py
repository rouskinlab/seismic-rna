from functools import partial
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

from .ct import parse_ct
from .db import parse_db
from .struct import RNAStructure
from ..path import CT_EXT, DB_EXT
from ..seq import Section
from ..write import need_write, write_mode

logger = getLogger(__name__)


def _from_file(file: Path, parser: Callable, *args, **kwargs):
    titles: set[str] = set()
    for title, section, pairs in parser(file, *args, **kwargs):
        if title in titles:
            logger.warning(f"Title {repr(title)} is repeated in {file}")
        else:
            titles.add(title)
        yield RNAStructure(title=title, section=section, pairs=pairs)


def from_ct(ct_path: Path):
    """ Yield an instance of an RNAStructure for each structure in a
    connectivity table (CT) file.

    Parameters
    ----------
    ct_path: Path
        Path of the CT file.

    Returns
    -------
    Generator[RNAStructure, Any, None]
        RNA secondary structures from the CT file.
    """
    yield from _from_file(ct_path, parse_ct)


def from_db(db_path: Path, seq5: int = 1):
    """ Yield an instance of an RNAStructure for each structure in a
    dot-bracket (DB) file.

    Parameters
    ----------
    db_path: Path
        Path of the DB file.
    seq5: int = 1
        Number to give the 5' position of the sequence.

    Returns
    -------
    Generator[RNAStructure, Any, None]
        RNA secondary structures from the CT file.
    """
    yield from _from_file(db_path, parse_db, seq5=seq5)


def find_ct_section(ct_path: Path) -> Section:
    """ Section shared among all structures in a CT file. """
    structures = iter(from_ct(ct_path))
    try:
        structure = next(structures)
    except StopIteration:
        raise ValueError(f"No structures in {ct_path}")
    section = structure.section
    for structure in structures:
        if structure.section != section:
            raise ValueError(f"Got > 1 unique section in {ct_path}")
    return section


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
        # Generate the output text for every structure.
        text = "".join(structure.ct_text for structure in structures)
        # Make the output directory, if it does not already exist.
        ct_path.parent.mkdir(parents=True, exist_ok=True)
        # Write the structures to the file.
        with open(ct_path, write_mode(force)) as f:
            f.write(text)
        logger.info(f"Wrote {ct_path}")


def to_db(structures: Iterable[RNAStructure],
          db_path: Path,
          force: bool = False):
    """ Write a dot-bracket (DB) file of RNA structures.

    Parameters
    ----------
    structures: Iterable[RNAStructure]
        RNA structures to write to the CT file.
    db_path: Path
        Path of the DB file.
    force: bool = False
        Overwrite the output DB file if it already exists.
    """
    if need_write(db_path, force):
        # Generate the output text for every structure.
        text = "".join(structure.get_db_text(i == 0)
                       for i, structure in enumerate(structures))
        # Make the output directory, if it does not already exist.
        db_path.parent.mkdir(parents=True, exist_ok=True)
        # Write the structures to the file.
        with open(db_path, write_mode(force)) as f:
            f.write(text)
        logger.info(f"Wrote {db_path}")


def renumber_ct(ct_in: Path, ct_out: Path, seq5: int, force: bool = False):
    """ Renumber the last column of a connectivity table (CT) file.

    Parameters
    ----------
    ct_in: Path
        Path of the input CT file.
    ct_out: Path
        Path of the output CT file.
    seq5: int
        Number to give the 5' position in the renumbered CT file.
    force: bool = False
        Overwrite the output CT file if it already exists.
    """
    to_ct(map(partial(RNAStructure.renumber_from, seq5=seq5),
              from_ct(ct_in)),
          ct_out,
          force)


def ct_to_db(ct_path: Path,
             db_path: Path | None = None,
             force: bool = False):
    """ Write a dot-bracket (DB) file of structures in a connectivity
    table (CT) file. """
    if db_path is None:
        db_path = ct_path.with_suffix(DB_EXT)
    to_db(from_ct(ct_path), db_path, force)
    return db_path


def db_to_ct(db_path: Path,
             ct_path: Path | None = None,
             force: bool = False):
    """ Write a connectivity table (CT) file of structures in a
    dot-bracket (DB) file. """
    if ct_path is None:
        ct_path = db_path.with_suffix(CT_EXT)
    to_ct(from_db(db_path), ct_path, force)
    return ct_path

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
