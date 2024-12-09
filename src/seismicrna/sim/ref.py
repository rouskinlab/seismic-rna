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
    return path.buildpar(path.FastaSeg,
                         top=top,
                         ref=ref,
                         ext=path.FASTA_EXTS[0])


@run_func(COMMAND, default=Path)
def run(*,
        sim_dir: str,
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
        logger.severe(error)

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
