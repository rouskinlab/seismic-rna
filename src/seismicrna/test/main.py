"""

Testing Main Module

========================================================================

This module is the entry point for the command line interface. Running

$ seismic [OPTIONS] command [OPTIONS] [ARGS]

calls the function main_cli() defined in this module.

"""

import unittest as ut
from os.path import dirname

from click import command

from ..core import docdef
from ..core.cli import opt_verbose
from ..core.cmd import CMD_TEST

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
