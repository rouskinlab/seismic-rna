from sys import stdout

from click import command

from .rnastructure import guess_data_path, DATAPATH
from ..core.arg import (CMD_DATAPATH)
from ..core.run import run_func


@run_func(CMD_DATAPATH, default=None)
def run_datapath():
    """ Guess the DATAPATH for RNAstructure. """
    datapath = guess_data_path()
    # This function should use stdout.write(), not the logger, because
    # the DATAPATH should be printed no matter the logging verbosity.
    stdout.write(f"{DATAPATH}={datapath}\n")
    return datapath


@command(CMD_DATAPATH)
def cli_datapath():
    """ Guess the DATAPATH for RNAstructure. """
    return run_datapath()

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
