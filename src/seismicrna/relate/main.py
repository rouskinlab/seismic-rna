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
from ..align.write import format_ref_minus
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
                        opt_sep_strands,
                        opt_minus_label,
                        opt_parallel,
                        opt_max_procs,
                        opt_force,
                        opt_keep_tmp)
from ..core.run import run_func
from ..core.seq import DNA, parse_fasta, write_fasta

logger = getLogger(__name__)


def generate_both_strands(ref: str, seq: DNA, minus_label: str):
    """ Yield both the plus and minus strand for each sequence. """
    yield ref, seq
    yield format_ref_minus(ref, minus_label), seq.rc


def write_both_strands(fasta_in: Path, fasta_out: Path, minus_label: str):
    """ Write a FASTA file of both plus and minus strands. """
    write_fasta(fasta_out,
                (strand
                 for ref, seq in parse_fasta(fasta_in, DNA)
                 for strand in generate_both_strands(ref, seq, minus_label)))


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
        sep_strands: bool,
        minus_label: str,
        max_procs: int,
        parallel: bool,
        brotli_level: int,
        force: bool,
        keep_tmp: bool):
    """ Compute relationships between references and aligned reads. """
    fasta = Path(fasta)
    if sep_strands:
        # Create a temporary FASTA file of both plus and minus strands.
        fasta_dir = tmp_dir.joinpath("fasta")
        fasta_dir.mkdir(parents=True, exist_ok=False)
        relate_fasta = fasta_dir.joinpath(fasta.name)
        write_both_strands(fasta, relate_fasta, minus_label)
    else:
        relate_fasta = fasta
    return write_all(xam_files=path.find_files_chain(map(Path, input_path),
                                                     path.XAM_SEGS),
                     fasta=relate_fasta,
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
    opt_sep_strands,
    opt_minus_label,
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
