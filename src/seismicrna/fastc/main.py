"""

FASTA Cleaning Module
========================================================================


"""

from logging import getLogger
from pathlib import Path

from click import command

from .fastaclean import FastaCleaner
from ..core import docdef
from ..core.cliparam import arg_fasta, opt_out_dir, opt_rerun
from ..core.clicmd import CMD_FASTC
from ..core.seq import DNA

logger = getLogger(__name__)

params = [
    arg_fasta,
    opt_out_dir,
    opt_rerun,
]


@command(CMD_FASTC, params=params)
def cli(*args, **kwargs):
    """ Clean the names and sequences in a FASTA file. """
    return run(*args, **kwargs)


@docdef.auto()
def run(fasta: str, out_dir: str, rerun: bool):
    """
    Clean a FASTA file.
    """
    try:
        fc = FastaCleaner(DNA)
        fc.run(fasta_path := Path(fasta),
               Path(out_dir).joinpath(fasta_path.name),
               force=rerun)
    except Exception as error:
        logger.critical(error)

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
