import cProfile
import os

from click import Context, group, pass_context, version_option

from . import (wf,
               demult,
               align,
               relate,
               pool,
               mask,
               cluster,
               join,
               table,
               lists,
               fold,
               graph,
               export,
               draw,
               migrate,
               test,
               sim,
               cleanfa,
               renumct,
               __version__)
from .align import split
from .core import rna
from .core.arg import (opt_log,
                       opt_log_color,
                       opt_profile,
                       opt_quiet,
                       opt_verbose)
from .core.logs import logger, set_config

params = [
    opt_verbose,
    opt_quiet,
    opt_log_color,
    opt_log,
    opt_profile,
]


# Group for main commands
@group(params=params, context_settings={"show_default": True})
@version_option(__version__)
@pass_context
def cli(ctx: Context,
        verbose: int,
        quiet: int,
        log_color: bool,
        log: str,
        profile: str,
        **kwargs):
    """ Command line interface of SEISMIC-RNA. """
    # Configure logging.
    if log:
        log_file_path = os.path.abspath(log)
        os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
    else:
        log_file_path = None
    set_config(verbose - quiet, log_file_path, log_color)
    logger.detail(f"This is SEISMIC-RNA version {__version__}")
    # If no subcommand was given, then run the entire pipeline.
    if ctx.invoked_subcommand is None:
        if profile:
            profile_path = os.path.abspath(profile)
            # Profile the program as it runs and write results to the
            # file given in the parameter profile.
            os.makedirs(os.path.dirname(profile_path), exist_ok=True)
            cProfile.runctx("wf.run(**kwargs)",
                            globals=globals(),
                            locals=locals(),
                            filename=profile_path,
                            sort="time")
        else:
            # Run without profiling.
            wf.run(**kwargs)


# Add all commands to the main CLI command group.
for module in (wf,
               demult,
               align,
               relate,
               pool,
               mask,
               cluster,
               join,
               table,
               fold,
               graph,
               export,
               draw,
               migrate,
               test,
               sim,
               cleanfa,
               renumct):
    cli.add_command(module.cli)
cli.add_command(split.cli)
cli.add_command(lists.listpos.cli)
cli.add_command(rna.convert.cli_ct2db)
cli.add_command(rna.convert.cli_db2ct)
cli.add_command(fold.cli_datapath)

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
