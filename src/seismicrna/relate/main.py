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

from .write import relate_all, get_relaters
from ..core import docdef, path
from ..core.cli import (arg_input_file, arg_fasta,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_min_reads, opt_batch_size, opt_ambrel,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_save_temp)
from ..core.parallel import lock_temp_dir

logger = getLogger(__name__)

# Parameters for command line interface
params = [
    # Input files
    arg_fasta,
    arg_input_file,
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # SAM options
    opt_phred_enc,
    opt_min_phred,
    # Vectoring options
    opt_min_reads,
    opt_batch_size,
    opt_ambrel,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # File generation
    opt_rerun,
    opt_save_temp,
]


@command(path.MOD_REL, params=params)
def cli(**kwargs):
    """ For every read, find the relationship between each read base and
    the reference base to which it aligned. """
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(fasta: str,
        input_file: tuple[str, ...],
        *,
        out_dir: str,
        temp_dir: str,
        phred_enc: int,
        min_phred: int,
        min_reads: int,
        batch_size: float,
        ambrel: bool,
        max_procs: int,
        parallel: bool,
        rerun: bool,
        save_temp: bool):
    """
    Run the relation step. For each read (or each pair of reads, if the
    reads are paired-end), generate a 'relation vector' that encodes the
    relationship between each base in the read and the corresponding
    base in the reference sequence -- that is, whether the base in the
    read matches the base in the reference, is substituted to another
    base, is deleted, or has one or more extra bases inserted beside it.
    """

    # For each BAM file, create a relation writer.
    relaters = get_relaters(path.find_files_chain(map(Path, input_file),
                                                  [path.SampSeg, path.XamSeg]),
                            Path(fasta),
                            min_reads=min_reads,
                            max_procs=max_procs,
                            parallel=parallel)

    # Compute and write relation vectors for each relation writer.
    return relate_all(relaters=relaters,
                      out_dir=Path(out_dir),
                      temp_dir=Path(temp_dir),
                      phred_enc=phred_enc,
                      min_phred=min_phred,
                      ambrel=ambrel,
                      batch_size=batch_size,
                      max_procs=max_procs,
                      parallel=parallel,
                      rerun=rerun,
                      save_temp=save_temp)
