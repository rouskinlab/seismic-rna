from logging import getLogger
from pathlib import Path

from click import command

from .demultiplex import demultiplex_run
from ..align.fqops import FastqUnit
from ..core.arg import (CMD_DEMULT,
                        opt_barcode_end,
                        opt_barcode_start,
                        opt_parallel_demultiplexing,
                        opt_clipped_demultiplexing,
                        opt_mismatch_tolerence,
                        opt_index_tolerence,
                        opt_demulti_overwrite,
                        arg_fasta,
                        opt_fastqx,
                        opt_out_dir,
                        opt_phred_enc,
                        opt_refs_meta,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        extra_defaults)
from ..core.run import run_func

logger = getLogger(__name__)


@run_func(logger.critical,
          with_tmp=True,
          pass_keep_tmp=True,
          extra_defaults=extra_defaults)
def run_dm(fasta: str,
           refs_meta: str,
           out_dir: str,
           tmp_dir: Path,
           fastqx: tuple[str, ...],
           phred_enc: int,
           barcode_start: int,
           barcode_end: int,
           clipped: int,
           index_tolerance: int,
           parallel_demultiplexing: bool,
           mismatch_tolerence: int,
           demulti_overwrite: bool,
           keep_tmp: bool):
    """ Split multiplexed FASTQ files by their barcodes. """
    fq_units = list(FastqUnit.from_paths(fastqx=list(map(Path, fastqx)),
                                         phred_enc=phred_enc))
    return [demultiplex_run(refs_file_csv=refs_meta,
                            overwrite=demulti_overwrite,
                            demulti_workspace=tmp_dir,
                            report_folder=out_dir,
                            fq_unit=fq_unit,
                            barcode_start=barcode_start,
                            barcode_end=barcode_end,
                            clipped=clipped,
                            index_tolerance=index_tolerance,
                            parallel=parallel_demultiplexing,
                            fasta=fasta,
                            mismatch_tolerence=mismatch_tolerence,
                            keep_tmp=keep_tmp)
            for fq_unit in fq_units]


params = [
    # Inputs
    arg_fasta,
    opt_fastqx,
    opt_phred_enc,
    opt_barcode_start,
    opt_barcode_end,
    opt_out_dir,
    # options
    opt_parallel_demultiplexing,
    opt_clipped_demultiplexing,
    opt_mismatch_tolerence,
    opt_index_tolerence,
    opt_demulti_overwrite,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_refs_meta,
]


# Turn into command.

@command(CMD_DEMULT, params=params)
def cli(*args, **kwargs):
    """ Split multiplexed FASTQ files by their barcodes. """
    return run_dm(*args, **kwargs)

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
