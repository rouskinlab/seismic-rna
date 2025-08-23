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
    # Write no log file, suppress warnings, and exit on errors.
    set_config(verbosity=Level.ERROR,
               log_file_path=None,
               exit_on_error=True)
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
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise RuntimeError(
            f"Some tests did not succeed ({len(result.failures)} failures, "
            f"{len(result.errors)} errors)"
        )


# Parameters for command line interface
params = [opt_verbose]


@command(CMD_TEST, params=params)
def cli(**kwargs):
    """ Run all unit tests. """
    return run(**kwargs)


if __name__ == "__main__":
    # Run all unit tests by executing this script on the command line.
    run(verbose=2)
