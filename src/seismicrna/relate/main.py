from pathlib import Path
from typing import Iterable

from click import command

from .strands import write_both_strands
from .write import relate_xam
from ..core import path
from ..core.arg import (CMD_RELATE,
                        arg_input_path,
                        arg_fasta,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_branch,
                        opt_min_mapq,
                        opt_min_reads,
                        opt_batch_size,
                        opt_phred_enc,
                        opt_min_phred,
                        opt_insert3,
                        opt_ambindel,
                        opt_overhangs,
                        opt_clip_end5,
                        opt_clip_end3,
                        opt_brotli_level,
                        opt_sep_strands,
                        opt_rev_label,
                        opt_write_read_names,
                        opt_relate_pos_table,
                        opt_relate_read_table,
                        opt_relate_cx,
                        opt_num_cpus,
                        opt_force,
                        opt_keep_tmp)
from ..core.logs import logger
from ..core.ngs import DuplicateSampleReferenceError
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch


def check_duplicates(xam_files: list[Path]):
    """ Check if any combination of sample, reference, and branches
    occurs more than once. """
    logger.routine("Began checking for duplicates")
    combos = set()
    for xam_file in xam_files:
        fields = path.parse(xam_file, path.XAM_SEGS)
        combo = (fields[path.SAMPLE],
                 fields[path.REF],
                 tuple(fields[path.BRANCHES]))
        logger.detail(f"{xam_file}: {combo}")
        if combo in combos:
            raise DuplicateSampleReferenceError(combo)
        combos.add(combo)
    logger.routine("Ended checking for duplicates")


@run_func(CMD_RELATE, with_tmp=True, pass_keep_tmp=True)
def run(fasta: str | Path,
        input_path: Iterable[str | Path], *,
        out_dir: str | Path,
        tmp_dir: Path,
        branch: str,
        min_reads: int,
        min_mapq: int,
        phred_enc: int,
        min_phred: int,
        batch_size: int,
        insert3: bool,
        ambindel: bool,
        overhangs: bool,
        clip_end5: int,
        clip_end3: int,
        sep_strands: bool,
        rev_label: str,
        write_read_names: bool,
        relate_pos_table: bool,
        relate_read_table: bool,
        relate_cx: bool,
        num_cpus: int,
        brotli_level: int,
        force: bool,
        keep_tmp: bool):
    """ Compute relationships between references and aligned reads. """
    fasta = Path(fasta)
    if sep_strands:
        # Create a temporary FASTA file of forward and reverse strands.
        fasta_dir = tmp_dir.joinpath("fasta")
        fasta_dir.mkdir(parents=True, exist_ok=False)
        relate_fasta = fasta_dir.joinpath(fasta.name)
        write_both_strands(fasta, relate_fasta, rev_label)
    else:
        relate_fasta = fasta
    # List the input XAM files and check for duplicates.
    xam_files = list(path.find_files_chain(input_path, path.XAM_SEGS))
    check_duplicates(xam_files)
    return dispatch(relate_xam,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(xam_files),
                    kwargs=dict(fasta=relate_fasta,
                                out_dir=Path(out_dir),
                                tmp_dir=tmp_dir,
                                branch=branch,
                                min_reads=min_reads,
                                min_mapq=min_mapq,
                                phred_enc=phred_enc,
                                min_phred=min_phred,
                                insert3=insert3,
                                ambindel=ambindel,
                                overhangs=overhangs,
                                clip_end5=clip_end5,
                                clip_end3=clip_end3,
                                batch_size=batch_size,
                                write_read_names=write_read_names,
                                relate_pos_table=relate_pos_table,
                                relate_read_table=relate_read_table,
                                relate_cx=relate_cx,
                                brotli_level=brotli_level,
                                force=force,
                                keep_tmp=keep_tmp))


# Parameters for command line interface
params = [
    # Input files
    arg_fasta,
    arg_input_path,
    opt_sep_strands,
    opt_rev_label,
    # Output directories
    opt_out_dir,
    opt_tmp_pfx,
    opt_branch,
    # SAM options
    opt_min_mapq,
    opt_phred_enc,
    opt_min_phred,
    # Relate options
    opt_min_reads,
    opt_batch_size,
    opt_insert3,
    opt_ambindel,
    opt_overhangs,
    opt_clip_end5,
    opt_clip_end3,
    opt_brotli_level,
    # Output options
    opt_write_read_names,
    opt_relate_pos_table,
    opt_relate_read_table,
    # Implementation options
    opt_relate_cx,
    # Parallelization
    opt_num_cpus,
    # File generation
    opt_force,
    opt_keep_tmp,
]


@command(CMD_RELATE, params=params)
def cli(**kwargs):
    """ Compute relationships between references and aligned reads. """
    return run(**kwargs)
