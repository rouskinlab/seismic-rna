from click import group

from . import seqbar, seqcorr, seqdiff, scatter
from ..core.arg import CMD_GRAPH


# Group for all graph commands
@group(CMD_GRAPH)
def cli():
    """ Graphing command line interface """


# Add graph commands to the CLI.
cli.add_command(seqbar.cli)
cli.add_command(seqcorr.cli)
cli.add_command(seqdiff.cli)
cli.add_command(scatter.cli)

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
