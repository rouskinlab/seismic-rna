from pathlib import Path

from click import command

from .fqunit import FastqUnit
from .write import align_samples
from ..core.arg import (CMD_ALIGN,
                        arg_fasta,
                        opt_fastqz,
                        opt_fastqy,
                        opt_fastqx,
                        opt_dmfastqz,
                        opt_dmfastqy,
                        opt_dmfastqx,
                        opt_phred_enc,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_force,
                        opt_keep_tmp,
                        opt_max_procs,
                        opt_fastp,
                        opt_fastp_5,
                        opt_fastp_3,
                        opt_fastp_w,
                        opt_fastp_m,
                        opt_fastp_poly_g,
                        opt_fastp_poly_g_min_len,
                        opt_fastp_poly_x,
                        opt_fastp_poly_x_min_len,
                        opt_fastp_adapter_trimming,
                        opt_fastp_adapter_1,
                        opt_fastp_adapter_2,
                        opt_fastp_adapter_fasta,
                        opt_fastp_detect_adapter_for_pe,
                        opt_fastp_min_length,
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
                        opt_bt2_un,
                        opt_min_mapq,
                        opt_min_reads,
                        opt_sep_strands,
                        opt_rev_label,
                        opt_f1r2_fwd,
                        optional_path,
                        extra_defaults)
from ..core.extern import (BOWTIE2_CMD,
                           BOWTIE2_BUILD_CMD,
                           FASTP_CMD,
                           SAMTOOLS_CMD,
                           require_dependency)
from ..core.run import run_func


@run_func(CMD_ALIGN,
          with_tmp=True,
          pass_keep_tmp=True,
          extra_defaults=extra_defaults)
def run(fasta: str, *,
        # Inputs
        fastqz: tuple[str, ...],
        fastqy: tuple[str, ...],
        fastqx: tuple[str, ...],
        dmfastqz: tuple[str, ...],
        dmfastqy: tuple[str, ...],
        dmfastqx: tuple[str, ...],
        phred_enc: int,
        # Outputs
        out_dir: str,
        tmp_dir: Path,
        keep_tmp: bool,
        # Fastp
        fastp: bool,
        fastp_5: bool,
        fastp_3: bool,
        fastp_w: int,
        fastp_m: int,
        fastp_poly_g: str,
        fastp_poly_g_min_len: int,
        fastp_poly_x: bool,
        fastp_poly_x_min_len: int,
        fastp_adapter_trimming: bool,
        fastp_adapter_1: str,
        fastp_adapter_2: str,
        fastp_adapter_fasta: str | None,
        fastp_detect_adapter_for_pe: bool,
        fastp_min_length: int,
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
        bt2_un: bool,
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
    if fastp:
        require_dependency(FASTP_CMD, __name__)
    require_dependency(BOWTIE2_CMD, __name__)
    require_dependency(BOWTIE2_BUILD_CMD, __name__)
    require_dependency(SAMTOOLS_CMD, __name__)
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
    return align_samples(
        fq_units=fq_units,
        fasta=Path(fasta),
        out_dir=Path(out_dir),
        tmp_dir=tmp_dir,
        keep_tmp=keep_tmp,
        force=force,
        max_procs=max_procs,
        fastp=fastp,
        fastp_5=fastp_5,
        fastp_3=fastp_3,
        fastp_w=fastp_w,
        fastp_m=fastp_m,
        fastp_poly_g=fastp_poly_g,
        fastp_poly_g_min_len=fastp_poly_g_min_len,
        fastp_poly_x=fastp_poly_x,
        fastp_poly_x_min_len=fastp_poly_x_min_len,
        fastp_adapter_trimming=fastp_adapter_trimming,
        fastp_adapter_1=fastp_adapter_1,
        fastp_adapter_2=fastp_adapter_2,
        fastp_adapter_fasta=optional_path(fastp_adapter_fasta),
        fastp_detect_adapter_for_pe=fastp_detect_adapter_for_pe,
        fastp_min_length=fastp_min_length,
        bt2_local=bt2_local,
        bt2_discordant=bt2_discordant,
        bt2_mixed=bt2_mixed,
        bt2_dovetail=bt2_dovetail,
        bt2_contain=bt2_contain,
        bt2_un=bt2_un,
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
        rev_label=rev_label
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
    # Outputs
    opt_out_dir,
    opt_tmp_pfx,
    opt_keep_tmp,
    # Fastp
    opt_fastp,
    opt_fastp_5,
    opt_fastp_3,
    opt_fastp_w,
    opt_fastp_m,
    opt_fastp_poly_g,
    opt_fastp_poly_g_min_len,
    opt_fastp_poly_x,
    opt_fastp_poly_x_min_len,
    opt_fastp_adapter_trimming,
    opt_fastp_adapter_1,
    opt_fastp_adapter_2,
    opt_fastp_adapter_fasta,
    opt_fastp_detect_adapter_for_pe,
    opt_fastp_min_length,
    # Bowtie2
    opt_bt2_local,
    opt_bt2_discordant,
    opt_bt2_mixed,
    opt_bt2_dovetail,
    opt_bt2_contain,
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
    opt_bt2_un,
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


@command(CMD_ALIGN, params=params)
def cli(*args, **kwargs):
    """ Trim FASTQ files and align them to reference sequences. """
    return run(*args, **kwargs)

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
