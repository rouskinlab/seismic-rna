import re
from logging import getLogger

logger = getLogger(__name__)

__version__ = "0.20.0"


def parse_version(version: str = __version__):
    """ Major and minor versions, patch, and pre-release tag. """
    match = re.match(r"^([0-9]+)[.]([0-9]+)[.]([0-9]+)([a-z]*[0-9]*)$", version)
    if not match:
        raise ValueError(f"Malformatted version: {repr(version)}")
    major, minor, patch = map(int, match.groups()[:3])
    prtag = match.groups()[3]
    return major, minor, patch, prtag


MAJOR, MINOR, PATCH, PRTAG = parse_version()


def format_version(major: int = MAJOR,
                   minor: int = MINOR,
                   patch: int = PATCH,
                   prtag: str = PRTAG):
    return f"{major}.{minor}.{patch}{prtag}"

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
