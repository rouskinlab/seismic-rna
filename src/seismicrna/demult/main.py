from pathlib import Path
from typing import Iterable

from click import command

from .demultiplex import demultiplex_run
from ..align.fqunit import FastqUnit
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


@run_func(CMD_DEMULT,
          with_tmp=True,
          pass_keep_tmp=True,
          extra_defaults=extra_defaults)
def run_dm(fasta: str | Path,
           refs_meta: str,
           out_dir: str | Path,
           tmp_dir: Path,
           fastqx: Iterable[str | Path],
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
    fq_units = list(FastqUnit.from_paths(fastqx=fastqx,
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
