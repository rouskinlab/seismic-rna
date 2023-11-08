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
                        docdef,
                        arg_input_path,
                        arg_fasta,
                        opt_out_dir,
                        opt_temp_dir,
                        opt_min_mapq,
                        opt_min_reads,
                        opt_batch_size,
                        opt_phred_enc,
                        opt_min_phred,
                        opt_ambrel,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs,
                        opt_force,
                        opt_keep_temp)
from ..core.parallel import lock_temp_dir

logger = getLogger(__name__)

# Parameters for command line interface
params = [
    # Input files
    arg_fasta,
    arg_input_path,
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # SAM options
    opt_min_mapq,
    opt_phred_enc,
    opt_min_phred,
    # Relate options
    opt_min_reads,
    opt_batch_size,
    opt_ambrel,
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # File generation
    opt_force,
    opt_keep_temp,
]


@command(CMD_REL, params=params)
def cli(**kwargs):
    """ For every read, find the relationship between each read base and
    the reference base to which it aligned. """
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(fasta: str,
        input_path: tuple[str, ...],
        *,
        out_dir: str,
        temp_dir: str,
        min_reads: int,
        min_mapq: int,
        phred_enc: int,
        min_phred: int,
        batch_size: float,
        ambrel: bool,
        max_procs: int,
        parallel: bool,
        brotli_level: int,
        force: bool,
        keep_temp: bool):
    """
    Run the relation step. For each read (or each pair of reads, if the
    reads are paired-end), generate a 'relation vector' that encodes the
    relationship between each base in the read and the corresponding
    base in the reference sequence -- that is, whether the base in the
    read matches the base in the reference, is substituted to another
    base, is deleted, or has one or more extra bases inserted beside it.
    """
    return write_all(xam_files=path.find_files_chain(map(Path, input_path),
                                                     path.XAM_SEGS),
                     fasta=Path(fasta),
                     out_dir=Path(out_dir),
                     temp_dir=Path(temp_dir),
                     min_reads=min_reads,
                     min_mapq=min_mapq,
                     phred_enc=phred_enc,
                     min_phred=min_phred,
                     ambrel=ambrel,
                     batch_size=batch_size,
                     max_procs=max_procs,
                     parallel=parallel,
                     brotli_level=brotli_level,
                     force=force,
                     keep_temp=keep_temp)

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
