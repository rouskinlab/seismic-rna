from sys import stdout

from click import command

from .rnastructure import guess_data_path, DATAPATH
from ..core.arg import (CMD_DATAPATH)
from ..core.logs import logger
from ..core.run import run_func


@run_func(logger.fatal, default=None)
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
