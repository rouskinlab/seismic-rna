"""

Testing Main Module

========================================================================

This module is the entry point for the command line interface. Running

$ seismic [OPTIONS] command [OPTIONS] [ARGS]

calls the function cli() defined in this module.

"""

import unittest as ut
from os.path import dirname

from click import command

from ..core.arg import CMD_TEST, docdef, opt_verbose

# Parameters for command line interface
params = [opt_verbose]


@command(CMD_TEST, params=params)
def cli(**kwargs):
    """ Run all unit tests. """
    return run(**kwargs)


@docdef.auto()
def run(verbose: int):
    """
    Run all unit tests.
    """
    # Discover all unit test modules.
    main_dir = dirname(dirname(__file__))
    # The line top_level_dir=dirname(main_dir) is needed to make Python
    # treat seismicrna as a package, so that the relative imports work.
    # Omitting this line leads to an ImportError during every test.
    suite = ut.TestLoader().discover(main_dir,
                                     pattern="*test.py",
                                     top_level_dir=dirname(main_dir))
    # Run all unit tests.
    runner = ut.TextTestRunner(verbosity=verbose)
    runner.run(suite)

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
