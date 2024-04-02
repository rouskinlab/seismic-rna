from shutil import which


def dependency_exists(dependency: str) -> bool:
    """ Check whether a dependency exists. """
    return which(dependency) is not None


def require_dependency(dependency: str, module: str = ""):
    """ If a dependency does not exist, return an error message. """
    if dependency_exists(dependency):
        # The dependency exists: no error message.
        return ""
    # The dependency is missing: error message.
    by = f"by '{module}' " if module else ""
    return (f"{repr(dependency)} is required {by}but was not found. "
            f"Please install it (if not yet) and place the executable for "
            f"{repr(dependency)} in your PATH.")

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
