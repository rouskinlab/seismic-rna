from logging import getLogger
from pathlib import Path

from click import command

from .fq2bam import get_bam_files
from .fqutil import FastqUnit
from ..core import docdef, path
from ..core.cli import (arg_fasta,
                        opt_fastqs, opt_fastqi, opt_fastqp,
                        opt_dmfastqs, opt_dmfastqi, opt_dmfastqp,
                        opt_phred_enc,
                        opt_out_dir, opt_temp_dir,
                        opt_rerun, opt_save_temp,
                        opt_parallel, opt_max_procs,
                        opt_fastqc, opt_qc_extract,
                        opt_cutadapt,
                        opt_cut_a1, opt_cut_g1, opt_cut_a2, opt_cut_g2,
                        opt_cut_o, opt_cut_e, opt_cut_q1, opt_cut_q2,
                        opt_cut_m, opt_cut_indels,
                        opt_cut_discard_trimmed, opt_cut_discard_untrimmed,
                        opt_cut_nextseq,
                        opt_bt2_local, opt_bt2_unal,
                        opt_bt2_discordant, opt_bt2_mixed,
                        opt_bt2_dovetail, opt_bt2_contain,
                        opt_bt2_i, opt_bt2_x,
                        opt_bt2_score_min_loc, opt_bt2_score_min_e2e,
                        opt_bt2_s, opt_bt2_l, opt_bt2_d, opt_bt2_r,
                        opt_bt2_gbar, opt_bt2_dpad, opt_bt2_orient)
from ..core.depend import confirm_dependency
from ..core.parallel import lock_temp_dir
from ..core.shell import (BOWTIE2_CMD, BOWTIE2_BUILD_CMD, CUTADAPT_CMD,
                          FASTQC_CMD, SAMTOOLS_CMD)

logger = getLogger(__name__)

# Parameters for command line interface
params = [
    # Inputs
    arg_fasta,
    opt_fastqp,
    opt_fastqi,
    opt_fastqs,
    opt_dmfastqp,
    opt_dmfastqi,
    opt_dmfastqs,
    opt_phred_enc,
    # Outputs
    opt_out_dir,
    opt_temp_dir,
    opt_rerun,
    opt_save_temp,
    # Parallelization
    opt_parallel,
    opt_max_procs,
    # FASTQC
    opt_fastqc,
    opt_qc_extract,
    # Cutadapt
    opt_cutadapt,
    opt_cut_a1,
    opt_cut_g1,
    opt_cut_a2,
    opt_cut_g2,
    opt_cut_o,
    opt_cut_e,
    opt_cut_q1,
    opt_cut_q2,
    opt_cut_m,
    opt_cut_indels,
    opt_cut_discard_trimmed,
    opt_cut_discard_untrimmed,
    opt_cut_nextseq,
    # Bowtie2
    opt_bt2_local,
    opt_bt2_discordant,
    opt_bt2_mixed,
    opt_bt2_dovetail,
    opt_bt2_contain,
    opt_bt2_unal,
    opt_bt2_i,
    opt_bt2_x,
    opt_bt2_score_min_e2e,
    opt_bt2_score_min_loc,
    opt_bt2_s,
    opt_bt2_l,
    opt_bt2_gbar,
    opt_bt2_d,
    opt_bt2_r,
    opt_bt2_dpad,
    opt_bt2_orient,
]


@command(path.MOD_ALIGN, params=params)
def cli(*args, **kwargs):
    """ Trim FASTQ files with Cutadapt (optional), validate them with
    FASTQC (optional), and align them to a FASTA file with Bowtie2. """
    return run(*args, **kwargs)


@lock_temp_dir
@docdef.auto()
def run(*,
        # Inputs
        fasta: str,
        fastqs: tuple[str, ...],
        fastqi: tuple[str, ...],
        fastqp: tuple[str, ...],
        dmfastqs: tuple[str, ...],
        dmfastqi: tuple[str, ...],
        dmfastqp: tuple[str, ...],
        phred_enc: int,
        # Outputs
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # FASTQC
        fastqc: bool,
        qc_extract: bool,
        # Cutadapt
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: tuple[str, ...],
        cut_a1: tuple[str, ...],
        cut_g2: tuple[str, ...],
        cut_a2: tuple[str, ...],
        cut_o: int,
        cut_e: float,
        cut_indels: bool,
        cut_nextseq: bool,
        cut_discard_trimmed: bool,
        cut_discard_untrimmed: bool,
        cut_m: int,
        # Bowtie2
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
        bt2_unal: bool,
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
        bt2_orient: str) -> list[Path]:
    """
    Run the alignment module.

    Align the reads to the set of reference sequences and output one BAM file
    for each sample aligned to each reference in the directory 'output'.
    Temporary intermediary files are written in the directory 'temp' and then
    deleted after they are no longer needed.
    """

    # Check for external dependencies.
    if fastqc:
        confirm_dependency(FASTQC_CMD, __name__)
    if cut:
        confirm_dependency(CUTADAPT_CMD, __name__)
    confirm_dependency(BOWTIE2_CMD, __name__)
    confirm_dependency(BOWTIE2_BUILD_CMD, __name__)
    confirm_dependency(SAMTOOLS_CMD, __name__)

    # FASTQ files of read sequences may come from up to seven different
    # sources (i.e. each argument beginning with "fq_unit"). This step
    # collects all of them into one list (fq_units) and also bundles
    # together pairs of FASTQ files containing mate 1 and mate 2 reads.
    fq_units = list(FastqUnit.from_paths(fastqs=list(map(Path, fastqs)),
                                         fastqi=list(map(Path, fastqi)),
                                         fastqp=list(map(Path, fastqp)),
                                         dmfastqs=list(map(Path, dmfastqs)),
                                         dmfastqi=list(map(Path, dmfastqi)),
                                         dmfastqp=list(map(Path, dmfastqp)),
                                         phred_enc=phred_enc))

    # Generate and return a BAM file for every FASTQ-reference pair.
    return get_bam_files(fq_units=fq_units,
                         fasta=Path(fasta),
                         out_dir=Path(out_dir),
                         temp_dir=Path(temp_dir),
                         save_temp=save_temp,
                         rerun=rerun,
                         max_procs=max_procs,
                         parallel=parallel,
                         fastqc=fastqc,
                         qc_extract=qc_extract,
                         cut=cut,
                         cut_q1=cut_q1,
                         cut_q2=cut_q2,
                         cut_g1=cut_g1,
                         cut_a1=cut_a1,
                         cut_g2=cut_g2,
                         cut_a2=cut_a2,
                         cut_o=cut_o,
                         cut_e=cut_e,
                         cut_indels=cut_indels,
                         cut_nextseq=cut_nextseq,
                         cut_discard_trimmed=cut_discard_trimmed,
                         cut_discard_untrimmed=cut_discard_untrimmed,
                         cut_m=cut_m,
                         bt2_local=bt2_local,
                         bt2_discordant=bt2_discordant,
                         bt2_mixed=bt2_mixed,
                         bt2_dovetail=bt2_dovetail,
                         bt2_contain=bt2_contain,
                         bt2_unal=bt2_unal,
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
                         bt2_orient=bt2_orient)
