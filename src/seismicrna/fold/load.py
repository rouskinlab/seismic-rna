from pathlib import Path
from typing import Iterable

from ..core import path
from ..core.rna import from_ct


def find_ct_files(files: Iterable[str | Path]):
    """ Yield a file for each given file/directory of a table. """
    yield from path.find_files_chain(map(Path, files), [path.ConnectTableSeg])


def load_ct_structs(files: Iterable[str | Path]):
    """ Yield an RNA structure generator for each CT file. """
    for file in find_ct_files(files):
        yield from_ct(file)

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
