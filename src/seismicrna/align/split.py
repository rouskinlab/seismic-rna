from pathlib import Path
from typing import Iterable

from click import command

from .fqunit import format_phred_arg
from .write import split_references
from ..core import path
from ..core.arg import (CMD_SPLITBAM,
                        arg_fasta,
                        arg_input_path,
                        opt_phred_enc,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_branch,
                        opt_force,
                        opt_keep_tmp,
                        opt_max_procs,
                        opt_bt2_local,
                        opt_bt2_discordant,
                        opt_bt2_mixed,
                        opt_bt2_dovetail,
                        opt_bt2_contain,
                        opt_bt2_i,
                        opt_bt2_x,
                        opt_bt2_score_min_loc,
                        opt_bt2_score_min_e2e,
                        opt_bt2_s,
                        opt_bt2_l,
                        opt_bt2_d,
                        opt_bt2_r,
                        opt_bt2_gbar,
                        opt_bt2_dpad,
                        opt_bt2_orient,
                        opt_min_mapq,
                        opt_min_reads,
                        opt_sep_strands,
                        opt_rev_label,
                        opt_f1r2_fwd)
from ..core.extern import (BOWTIE2_CMD,
                           BOWTIE2_BUILD_CMD,
                           SAMTOOLS_CMD,
                           require_dependency)
from ..core.ngs import (run_flagstat,
                        run_sort_xam,
                        run_index_xam,
                        xam_paired)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import release_to_out
from ..core.write import need_write


def split_xam_file(xam_file: Path,
                   out_dir: Path,
                   tmp_dir: Path,
                   branch: str,
                   fasta: Path,
                   phred_enc: int,
                   force: bool,
                   n_procs: int,
                   **kwargs):
    # Assume the XAM file is named for the sample.
    sample = xam_file.stem
    branches = path.add_branch(path.ALIGN_STEP, branch, dict())
    # Determine the final output directory.
    result_dir = path.build(path.CMD_DIR_SEGS,
                            {path.TOP: out_dir,
                             path.SAMPLE: sample,
                             path.CMD: path.ALIGN_STEP,
                             path.BRANCHES: branches})
    if need_write(result_dir, force):
        # Sort and index the XAM file.
        xam_input_dir = tmp_dir.joinpath("input")
        xam_sorted = path.buildpar(path.XAM_SEGS,
                                   {path.TOP: xam_input_dir,
                                    path.SAMPLE: sample,
                                    path.CMD: path.ALIGN_STEP,
                                    path.BRANCHES: branches,
                                    path.REF: fasta.stem,
                                    path.EXT: xam_file.suffix})
        run_sort_xam(xam_file, xam_sorted, n_procs=n_procs)
        run_index_xam(xam_sorted, n_procs=n_procs)
        # Split the XAM file into one file for each reference.
        release_dir = tmp_dir.joinpath("release")
        paired = xam_paired(run_flagstat(xam_file, n_procs=n_procs))
        split_references(xam_sorted,
                         fasta=fasta,
                         paired=paired,
                         phred_arg=format_phred_arg(phred_enc),
                         top=release_dir,
                         branches=branches,
                         n_procs=n_procs,
                         **kwargs)
        release_to_out(out_dir,
                       release_dir,
                       path.build(path.CMD_DIR_SEGS,
                                  {path.TOP: release_dir,
                                   path.SAMPLE: sample,
                                   path.CMD: path.ALIGN_STEP,
                                   path.BRANCHES: branches}))
    return result_dir


@run_func(CMD_SPLITBAM, with_tmp=True, pass_keep_tmp=True)
def run(fasta: str | Path, *,
        # Inputs
        input_path: Iterable[str | Path],
        phred_enc: int,
        # Outputs
        out_dir: str | Path,
        tmp_dir: Path,
        keep_tmp: bool,
        branch: str,
        # Bowtie2
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
        bt2_score_min_e2e: str,
        bt2_score_min_loc: str,
        bt2_i: int,
        bt2_x: int,
        bt2_gbar: int,
        bt2_l: int,
        bt2_s: str,
        bt2_d: int,
        bt2_r: int,
        bt2_dpad: int,
        bt2_orient: str,
        # Samtools
        min_mapq: int,
        min_reads: int,
        sep_strands: bool,
        f1r2_fwd: bool,
        rev_label: str,
        # Parallelization
        max_procs: int,
        force: bool) -> list[Path]:
    """ Trim FASTQ files and align them to reference sequences. """
    # Check for external dependencies.
    require_dependency(SAMTOOLS_CMD, __name__)
    require_dependency(BOWTIE2_CMD, __name__)
    require_dependency(BOWTIE2_BUILD_CMD, __name__)
    # Split each input XAM file.
    return dispatch(split_xam_file,
                    max_procs=max_procs,
                    pass_n_procs=True,
                    args=as_list_of_tuples(map(Path, input_path)),
                    kwargs=dict(fasta=Path(fasta),
                                phred_enc=phred_enc,
                                out_dir=Path(out_dir),
                                tmp_dir=tmp_dir,
                                keep_tmp=keep_tmp,
                                branch=branch,
                                force=force,
                                bt2_local=bt2_local,
                                bt2_discordant=bt2_discordant,
                                bt2_mixed=bt2_mixed,
                                bt2_dovetail=bt2_dovetail,
                                bt2_contain=bt2_contain,
                                bt2_score_min_e2e=bt2_score_min_e2e,
                                bt2_score_min_loc=bt2_score_min_loc,
                                bt2_i=bt2_i,
                                bt2_x=bt2_x,
                                bt2_gbar=bt2_gbar,
                                bt2_l=bt2_l,
                                bt2_s=bt2_s,
                                bt2_d=bt2_d,
                                bt2_r=bt2_r,
                                bt2_dpad=bt2_dpad,
                                bt2_orient=bt2_orient,
                                min_mapq=min_mapq,
                                min_reads=min_reads,
                                sep_strands=sep_strands,
                                f1r2_fwd=f1r2_fwd,
                                rev_label=rev_label))


# Parameters for command line interface
params = [
    # Inputs
    arg_fasta,
    arg_input_path,
    opt_phred_enc,
    # Outputs
    opt_out_dir,
    opt_tmp_pfx,
    opt_force,
    opt_keep_tmp,
    opt_branch,
    # Bowtie2
    opt_bt2_local,
    opt_bt2_discordant,
    opt_bt2_mixed,
    opt_bt2_dovetail,
    opt_bt2_contain,
    opt_bt2_i,
    opt_bt2_x,
    opt_bt2_score_min_loc,
    opt_bt2_score_min_e2e,
    opt_bt2_s,
    opt_bt2_l,
    opt_bt2_d,
    opt_bt2_r,
    opt_bt2_gbar,
    opt_bt2_dpad,
    opt_bt2_orient,
    # Samtools
    opt_min_mapq,
    opt_min_reads,
    opt_sep_strands,
    opt_f1r2_fwd,
    opt_rev_label,
    # Parallelization
    opt_max_procs,
    opt_force,
]


@command(CMD_SPLITBAM, params=params)
def cli(*args, **kwargs):
    """ Trim FASTQ files and align them to reference sequences. """
    return run(*args, **kwargs)
