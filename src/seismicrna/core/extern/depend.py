from shutil import which


def dependency_exists(dependency: str) -> bool:
    """ Check whether a dependency exists. """
    return which(dependency) is not None


def require_dependency(dependency: str, module: str = ""):
    """ Check whether a dependency exists and raise an error if not. """
    if not dependency_exists(dependency):
        by = f"by '{module}' " if module else ""
        raise RuntimeError(f"The dependency '{dependency}' is required {by}but "
                           f"was not found. Please install it (if not yet) and "
                           f"place the executable '{dependency}' in your PATH.")

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
