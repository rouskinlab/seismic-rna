from logging import getLogger
from pathlib import Path

logger = getLogger(__name__)


def need_write(query: Path, force: bool = False, warn: bool = True):
    """ Determine whether a file/directory must be written.

    Parameters
    ----------
    query: Path
        File or directory for which to check the need for writing.
    force: bool = False
        Force the query to be written, even if it already exists.
    warn: bool = True
        If the query does not need to be written, then log a warning.

    Returns
    -------
    bool
        Whether the file must be written.
    """
    if force or not query.exists():
        return True
    if warn:
        logger.warning(f"{query} exists: use --force to overwrite")
    return False


def write_mode(force: bool = False, binary: bool = False):
    """ Get the mode in which to open a file for writing.

    Parameters
    ----------
    force: bool = False
        Force the file to be written, truncating the file if it exists.
        If False and the file exists, a FileExistsError will be raised.
    binary: bool = False
        Write the file in binary mode instead of text mode.

    Returns
    -------
    str
        The mode argument for the builtin function `open()`.
    """
    return "".join(["w" if force else "x", "b" if binary else ""])

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
