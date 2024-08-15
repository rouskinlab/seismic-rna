from logging import getLogger

from click import command

from .rnastructure import guess_data_path, DATAPATH
from ..core.arg import (CMD_DATAPATH)
from ..core.run import run_func

logger = getLogger(__name__)


@run_func(logger.critical, default=None)
def run_datapath():
    """ Guess the DATAPATH for RNAstructure. """
    datapath = guess_data_path()
    # This function should use print(), not the logger, because the
    # DATAPATH should be printed regardless of the logging verbosity.
    print(f"{DATAPATH}={datapath}")
    return datapath


@command(CMD_DATAPATH)
def cli_datapath():
    """ Guess the DATAPATH for RNAstructure. """
    return run_datapath()
