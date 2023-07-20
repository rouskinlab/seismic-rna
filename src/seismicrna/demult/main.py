from click import command
from pathlib import Path

from ..demult.demultiplex import demultiplex_run
from ..align.fqutil import FastqUnit
from ..core.cli import (
    opt_barcode_length, opt_barcode_start, opt_parallel_demultiplexing, opt_clipped_demultiplexing,
    opt_mismatch_tolerence, opt_index_tolerence, opt_demulti_overwrite, arg_fasta, opt_sections_file, opt_fastqp, opt_out_dir,
    opt_phred_enc)

params = [
    # Inputs
    arg_fasta,
    opt_fastqp,
    opt_phred_enc,
    opt_sections_file,
    opt_barcode_start,
    opt_barcode_length,
    opt_out_dir,

    # options
    opt_parallel_demultiplexing,
    opt_clipped_demultiplexing,
    opt_mismatch_tolerence,
    opt_index_tolerence,
    opt_demulti_overwrite

]


# Turn into DREEM command.

@command("demultiplex", params=params)
def cli(*args, **kwargs):
    """ Split multiplexed FASTQ files by their barcodes. """
    return run(*args, **kwargs)


def run(sections_file: str, out_dir: str, temp_dir: str, fastqp: tuple[str, ...], phred_enc: int, fasta: str, barcode_start=0,
        barcode_length=0, clipped: int = 0, index_tolerance: int = 0, parallel_demultiplexing: bool = False,
        mismatch_tolerence: int = 0, demulti_overwrite: bool = False):
    fq_units = list(FastqUnit.from_paths(fastqp=list(map(Path, fastqp)),
                                         phred_enc=phred_enc))
    return [demultiplex_run(sections_file_csv=sections_file,
                            overwrite=demulti_overwrite,
                            demulti_workspace=temp_dir,
                            report_folder=out_dir,
                            fq_unit=fq_unit,
                            barcode_start=barcode_start,
                            barcode_length=barcode_length,
                            clipped=clipped,
                            index_tolerance=index_tolerance,
                            parallel=parallel_demultiplexing,
                            fasta=fasta,
                            mismatch_tolerence=mismatch_tolerence)
            for fq_unit in fq_units]


if __name__ == '__main__':
    pass
