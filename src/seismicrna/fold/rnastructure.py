"""
Struct -- RNAstructure Module

Wrapper around RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

from logging import getLogger
from pathlib import Path

from ..core import path
from ..core.rna.profile import RnaProfile
from ..core.extern import (RNASTRUCTURE_CT2DOT_CMD,
                           RNASTRUCTURE_DOT2CT_CMD,
                           RNASTRUCTURE_FOLD_CMD,
                           args_to_cmd,
                           run_cmd, )

logger = getLogger(__name__)


def fold(rna: RnaProfile, *,
         out_dir: Path, temp_dir: Path, save_temp: bool, force: bool):
    """ Run the 'Fold' program of RNAstructure. """
    ct_file = rna.ct_file(out_dir)
    if force or not ct_file.is_file():
        cmd = [RNASTRUCTURE_FOLD_CMD]
        # Write the DMS reactivities file for the RNA.
        cmd.extend(["--DMS", rna.to_dms(temp_dir)])
        # Write a temporary FASTA file for the RNA.
        cmd.append(fasta := rna.to_fasta(temp_dir))
        try:
            # Get the path of the output CT file.
            cmd.append(ct_file)
            # Run the command.
            run_cmd(args_to_cmd(cmd))
        finally:
            if not save_temp:
                # Delete the temporary files.
                fasta.unlink(missing_ok=True)
    else:
        logger.warning(f"File exists: {ct_file}")
    return ct_file


def ct2dot(ct_file: Path, number: int | str = "all"):
    """ Convert a CT file to a DOT file. """
    dot_file = ct_file.with_suffix(path.DOT_EXT)
    cmd = [RNASTRUCTURE_CT2DOT_CMD, ct_file, number, dot_file]
    run_cmd(args_to_cmd(cmd))
    return dot_file


def dot2ct(dot_file: Path):
    """ Convert a DOT file to a CT file. """
    ct_file = dot_file.with_suffix(path.CT_EXT)
    cmd = [RNASTRUCTURE_DOT2CT_CMD, dot_file, ct_file]
    run_cmd(args_to_cmd(cmd))
    return ct_file

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
