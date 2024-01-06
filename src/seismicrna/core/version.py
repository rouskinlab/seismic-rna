"""
Version information for SEISMIC-RNA
"""

from logging import getLogger

logger = getLogger(__name__)

__version__ = "0.11.3"

VERSION_DELIM = "."


def parse_version(version: str = __version__):
    try:
        major, minor, patch = map(int, version.split(VERSION_DELIM))
    except ValueError:
        raise ValueError(f"Malformatted version number: {version}")
    return major, minor, patch


MAJOR, MINOR, PATCH = parse_version()


def format_version(major: int = MAJOR, minor: int = MINOR, patch: int = PATCH):
    return VERSION_DELIM.join(map(str, (major, minor, patch)))

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
