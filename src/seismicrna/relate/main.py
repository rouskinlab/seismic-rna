"""
Relate -- Main Module
=====================
Auth: Matty

Define the command line interface for the 'relate' command, as well as
its main run function that executes the relate step.
"""

from logging import getLogger
from pathlib import Path

from click import command

from .write import write_all
from ..core import path
from ..core.arg import (CMD_REL,
                        arg_input_path,
                        arg_fasta,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_min_mapq,
                        opt_min_reads,
                        opt_batch_size,
                        opt_phred_enc,
                        opt_min_phred,
                        opt_ambindel,
                        opt_overhangs,
                        opt_clip_end5,
                        opt_clip_end3,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs,
                        opt_force,
                        opt_keep_tmp)
from ..core.run import run_func

logger = getLogger(__name__)


@run_func(logger.critical, with_tmp=True, pass_keep_tmp=True)
def run(fasta: str,
        input_path: tuple[str, ...], *,
        out_dir: str,
        tmp_dir: Path,
        min_reads: int,
        min_mapq: int,
        phred_enc: int,
        min_phred: int,
        batch_size: int,
        ambindel: bool,
        overhangs: bool,
        clip_end5: int,
        clip_end3: int,
        max_procs: int,
        parallel: bool,
        brotli_level: int,
        force: bool,
        keep_tmp: bool):
    """ Compute relationships between references and aligned reads. """
    return write_all(xam_files=path.find_files_chain(map(Path, input_path),
                                                     path.XAM_SEGS),
                     fasta=Path(fasta),
                     out_dir=Path(out_dir),
                     tmp_dir=tmp_dir,
                     min_reads=min_reads,
                     min_mapq=min_mapq,
                     phred_enc=phred_enc,
                     min_phred=min_phred,
                     ambindel=ambindel,
                     overhangs=overhangs,
                     clip_end5=clip_end5,
                     clip_end3=clip_end3,
                     batch_size=batch_size,
                     max_procs=max_procs,
                     parallel=parallel,
                     brotli_level=brotli_level,
                     force=force,
                     keep_tmp=keep_tmp)


# Parameters for command line interface
params = [
    # Input files
    arg_fasta,
    arg_input_path,
    # Output directories
    opt_out_dir,
    opt_tmp_pfx,
    # SAM options
    opt_min_mapq,
    opt_phred_enc,
    opt_min_phred,
    # Relate options
    opt_min_reads,
    opt_batch_size,
    opt_ambindel,
    opt_overhangs,
    opt_clip_end5,
    opt_clip_end3,
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # File generation
    opt_force,
    opt_keep_tmp,
]


@command(CMD_REL, params=params)
def cli(**kwargs):
    """ Compute relationships between references and aligned reads. """
    return run(**kwargs)

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
