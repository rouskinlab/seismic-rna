import unittest as ut
from os.path import dirname

from click import command

from seismicrna.core.arg import CMD_TEST, opt_verbose
from seismicrna.core.logs import Level, restore_config, set_config
from seismicrna.core.run import run_func


@run_func(CMD_TEST, default=None)
@restore_config
def run(verbose: int):
    """ Run all unit tests. """
    # Write no log file, suppress warnings, and halt on errors.
    set_config(verbosity=Level.ERROR,
               log_file_path=None,
               raise_on_error=True)
    # Discover all unit test modules.
    main_dir = dirname(dirname(__file__))
    # The argument top_level_dir=dirname(main_dir) is needed to make
    # Python treat seismicrna as a package, so relative imports work.
    # Omitting this argument causes an ImportError during every test.
    suite = ut.TestLoader().discover(main_dir,
                                     pattern="*test.py",
                                     top_level_dir=dirname(main_dir))
    # Run all unit tests.
    runner = ut.TextTestRunner(verbosity=verbose)
    runner.run(suite)


# Parameters for command line interface
params = [opt_verbose]


@command(CMD_TEST, params=params)
def cli(**kwargs):
    """ Run all unit tests. """
    return run(**kwargs)


if __name__ == "__main__":
    # Run all unit tests by executing this script on the command line.
    run(verbose=0)

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
