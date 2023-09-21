from hashlib import md5
from logging import getLogger
from pathlib import Path

logger = getLogger(__name__)


def digest_file(file_path: Path) -> str:
    """
    Compute the checksum of a file.

    Parameters
    ----------
    file_path: Path
        Path of the file on which to compute the checksum. Can be
        any type that the open() function recognizes as a path.

    Returns
    -------
    str
        Checksum of the file (in hexadecimal)
    """
    with open(file_path, "rb") as f:
        digest = md5(f.read()).hexdigest()
    logger.debug(f"Computed MD5 digest of {file_path}: {digest}")
    return digest

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                              #
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
