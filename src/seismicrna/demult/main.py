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
                        opt_temp_dir,
                        opt_keep_temp, )
from ..core.parallel import lock_temp_dir

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
    opt_temp_dir,
    opt_keep_temp,
    opt_refs_meta,
]


# Turn into DREEM command.

@command(CMD_DEMULT, params=params)
def cli(*args, **kwargs):
    """ Split multiplexed FASTQ files by their barcodes. """
    return run(*args, **kwargs)


@lock_temp_dir
def run(refs_meta: str,
        out_dir: str,
        temp_dir: str,
        fastqx: tuple[str, ...],
        phred_enc: int,
        fasta: str,
        barcode_start=0,
        barcode_end=0,
        clipped: int = 0,
        index_tolerance: int = 0,
        parallel_demultiplexing: bool = False,
        mismatch_tolerence: int = 0,
        demulti_overwrite: bool = False,
        keep_temp: bool = True):
    """ Split multiplexed FASTQ files by their barcodes. """
    fq_units = list(FastqUnit.from_paths(fastqx=list(map(Path, fastqx)),
                                         phred_enc=phred_enc))
    return [demultiplex_run(refs_file_csv=refs_meta,
                            overwrite=demulti_overwrite,
                            demulti_workspace=temp_dir,
                            report_folder=out_dir,
                            fq_unit=fq_unit,
                            barcode_start=barcode_start,
                            barcode_end=barcode_end,
                            clipped=clipped,
                            index_tolerance=index_tolerance,
                            parallel=parallel_demultiplexing,
                            fasta=fasta,
                            mismatch_tolerence=mismatch_tolerence,
                            keep_temp=keep_temp)
            for fq_unit in fq_units]


if __name__ == '__main__':
    pass

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
