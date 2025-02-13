import os
from pathlib import Path

from click import command

from ..core import path
from ..core.arg import (opt_sim_dir,
                        opt_ref,
                        opt_refs,
                        opt_reflen,
                        opt_force)
from ..core.logs import logger
from ..core.run import run_func
from ..core.seq import DNA, write_fasta
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]


def get_fasta_path(top: Path, ref: str):
    """ Get the path of a FASTA file. """
    return path.buildpar([path.FastaSeg],
                         {path.TOP: top,
                          path.REF: ref,
                          path.EXT: path.FASTA_EXTS[0]})


@run_func(COMMAND, default=Path)
def run(*,
        sim_dir: str | Path,
        refs: str,
        ref: str,
        reflen: int,
        force: bool):
    top = Path(sim_dir).joinpath(path.SIM_REFS_DIR)
    top.mkdir(parents=True, exist_ok=True)
    fasta = get_fasta_path(top, refs)
    if need_write(fasta, force):
        seq = DNA.random(reflen)
        write_fasta(fasta, [(ref, seq)], force=force)
    return fasta


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
    try:
        run(*args, **kwargs)
    except Exception as error:
        logger.fatal(error)
