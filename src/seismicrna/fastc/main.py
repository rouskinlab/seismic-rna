"""

FASTA Cleaning Module
========================================================================


"""

from logging import getLogger
from pathlib import Path

from click import command

from .fastaclean import FastaCleaner
from ..core import docdef
from ..core.cli import arg_fasta, opt_out_dir, opt_rerun
from ..core.cmd import CMD_FASTC
from ..core.depend import confirm_dependency
from ..core.seq import DNA
from ..core.shell import GREP_CMD

logger = getLogger(__name__)

params = [
    arg_fasta,
    opt_out_dir,
    opt_rerun,
]


@command(CMD_FASTC, params=params)
def cli(*args, **kwargs):
    """ Clean the names and sequences in a FASTA file. """
    return run(*args, **kwargs)


@docdef.auto()
def run(fasta: str, out_dir: str, rerun: bool):
    """
    Clean a FASTA file.
    """
    try:
        confirm_dependency(GREP_CMD, __name__)
        fc = FastaCleaner(DNA)
        fc.run(fasta_path := Path(fasta),
               Path(out_dir).joinpath(fasta_path.name),
               force=rerun)
    except Exception as error:
        logger.critical(error)
