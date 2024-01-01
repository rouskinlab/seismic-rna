from functools import partial
from logging import getLogger
from pathlib import Path
from typing import Iterable

from .ct import parse_ct
from .struct import RNAStructure
from ..write import need_write, write_mode

logger = getLogger(__name__)


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
    titles: set[str] = set()
    for title, section, pairs in parse_ct(ct_path):
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
        # Generate the output text for every renumbered structure.
        text = "".join(structure.ct_text for structure in structures)
        # Make the output directory, if it does not already exist.
        ct_path.parent.mkdir(parents=True, exist_ok=True)
        # Write the numbered structures to the file.
        with open(ct_path, write_mode(force)) as f:
            f.write(text)


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
