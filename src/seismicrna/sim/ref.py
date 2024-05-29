import os
from logging import getLogger
from pathlib import Path

from click import command

from ..core import path
from ..core.arg import (docdef,
                        opt_sim_dir,
                        opt_ref,
                        opt_refs,
                        opt_reflen,
                        opt_force)
from ..core.seq import DNA, write_fasta
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]

logger = getLogger(__name__)


@docdef.auto()
def run(sim_dir: str,
        refs: str,
        ref: str,
        reflen: int,
        force: bool):
    try:
        fasta = Path(sim_dir).joinpath(refs).with_suffix(path.FASTA_EXTS[0])
        if need_write(fasta, force):
            seq = DNA.random(reflen)
            fasta.parent.mkdir(parents=True, exist_ok=True)
            write_fasta(fasta, [(ref, seq)], force=force)
        return fasta
    except Exception as error:
        logger.critical(f"Failed to create {repr(ref)} ({reflen} nt): {error}")


params = [
    opt_sim_dir,
    opt_refs,
    opt_ref,
    opt_reflen,
    opt_force
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate a FASTA file of a reference sequence. """
    run(*args, **kwargs)
