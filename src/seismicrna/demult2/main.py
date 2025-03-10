from pathlib import Path
from typing import Iterable

from click import command

from ..align.fqunit import FastqUnit
from .write import demult_samples
from ..core.arg import (CMD_DEMULT2,
                        arg_fasta,
                        opt_fastqz,
                        opt_fastqy,
                        opt_fastqx,
                        opt_dmfastqz,
                        opt_dmfastqy,
                        opt_dmfastqx,
                        opt_refs_meta,
                        opt_barcode_start,
                        opt_barcode_end,
                        opt_read_pos,
                        opt_barcode,
                        opt_mismatch_tolerance,
                        opt_index_tolerance,
                        opt_phred_enc,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_branch,
                        opt_force,
                        opt_keep_tmp,
                        opt_max_procs,
                        extra_defaults)
from ..core.extern import (SEQKIT_CMD,
                           require_dependency)
from ..core.run import run_func
from ..core.seq.xna import DNA


@run_func(CMD_DEMULT2,
          with_tmp=True,
          pass_keep_tmp=True,
          extra_defaults=extra_defaults)
def run(fasta: str | Path, *,
        # Inputs
        fastqz: Iterable[str | Path],
        fastqy: Iterable[str | Path],
        fastqx: Iterable[str | Path],
        dmfastqz: Iterable[str | Path],
        dmfastqy: Iterable[str | Path],
        dmfastqx: Iterable[str | Path],
        phred_enc: int,
        # Barcodes
        refs_meta: str | Path,
        barcode_start: int,
        barcode_end: int,
        read_pos: int,
        barcode: tuple[tuple[str, DNA, int]],
        mismatch_tolerance: int,
        index_tolerance: int,
        # Outputs
        out_dir: str | Path,
        tmp_dir: Path,
        keep_tmp: bool,
        branch: str,
        # Parallelization
        max_procs: int,
        force: bool) -> list[Path]:
    """ Demultiplex FASTQ files. """
    # Check for external dependencies.
    require_dependency(SEQKIT_CMD, __name__)
    # FASTQ files of read sequences may come from up to seven different
    # sources (i.e. each argument beginning with "fq_unit"). This step
    # collects all of them into one list (fq_units) and also bundles
    # together pairs of FASTQ files containing mate 1 and mate 2 reads.
    fq_units = list(FastqUnit.from_paths(fastqz=list(map(Path, fastqz)),
                                         fastqy=list(map(Path, fastqy)),
                                         fastqx=list(map(Path, fastqx)),
                                         dmfastqz=list(map(Path, dmfastqz)),
                                         dmfastqy=list(map(Path, dmfastqy)),
                                         dmfastqx=list(map(Path, dmfastqx)),
                                         phred_enc=phred_enc))
    # Generate and return a BAM file for every FASTQ-reference pair.
    return demult_samples(
        fq_units=fq_units,
        fasta=Path(fasta),
        out_dir=Path(out_dir),
        refs_meta=refs_meta,
        barcode_start=barcode_start,
        barcode_end=barcode_end,
        read_pos=read_pos,
        barcode=barcode,
        mismatch_tolerance=mismatch_tolerance,
        index_tolerance=index_tolerance,
        tmp_dir=tmp_dir,
        keep_tmp=keep_tmp,
        branch=branch,
        force=force,
        max_procs=max_procs,
    )


# Parameters for command line interface
params = [
    # Inputs
    arg_fasta,
    opt_fastqx,
    opt_fastqy,
    opt_fastqz,
    opt_dmfastqx,
    opt_dmfastqy,
    opt_dmfastqz,
    opt_phred_enc,
    # Barcodes
    opt_refs_meta,
    opt_barcode_start,
    opt_barcode_end,
    opt_read_pos,
    opt_barcode,
    opt_mismatch_tolerance,
    opt_index_tolerance,
    # Outputs
    opt_out_dir,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_branch,
    # Parallelization
    opt_max_procs,
    opt_force,
]


@command(CMD_DEMULT2, params=params)
def cli(*args, **kwargs):
    """ Trim FASTQ files and align them to reference sequences. """
    return run(*args, **kwargs)
