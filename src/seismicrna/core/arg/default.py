from inspect import getmembers

from click import Argument, Option

from . import cli

# Get every parameter defined for the command line interface.
cli_args = dict(getmembers(cli, lambda member: isinstance(member, Argument)))
cli_opts = dict(getmembers(cli, lambda member: isinstance(member, Option)))

# Get the default value for every parameter.
cli_defaults = {param.name: param.default
                for param in (cli_args | cli_opts).values()
                if param.default is not None}

extra_defaults = dict(mask_sections_file=None,
                      mask_pos_file=None,
                      join_clusts=None,
                      fold_sections_file=None,
                      fold_constraint=None,
                      samples_meta=None,
                      refs_meta=None)

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
