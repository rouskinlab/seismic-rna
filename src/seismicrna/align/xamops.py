"""

Alignment XAM Generation Module

========================================================================

Alignment Score Parameters for Bowtie2

Consider this example: Ref = ACGT, Read = AG

Assume that we want to minimize the number of edits needed to convert
the reference into the read sequence. The smallest number of edits is
two, specifically these two deletions (/) from the reference: [A/G/]
which gets a score of (2 * match - 2 * gap_open - 2 * gap_extend).

But there are two alternative alignments, each with 3 edits:
[Ag//] and [A//g] (substitutions marked in lowercase). Each gets the
score (match - substitution - gap_open - 2 * gap_extend).

In order to favor the simpler alignment with two edits,
(2 * match - 2 * gap_open - 2 * gap_extend) must be greater than
(match - substitution - gap_open - 2 * gap_extend). This inequality
simplifies to (substitution > gap_open - match).

Thus, the substitution penalty and match bonus must be relatively large,
and the gap open penalty small. We want to avoid introducing too many
gaps, especially to prevent the introduction of an insertion and a
deletion from scoring better than one substitution.

Consider this example: Ref = ATAT, Read = ACTT

The simplest alignment (the smallest number of mutations) is ActT, which
gets a score of (2 * match - 2 * substitution). Another alignment with
indels is A{C}T/T, where {C} means a C was inserted into the read and
the / denotes an A deleted from the read. This alignment scores
(3 * match - 2 * gap_open - 2 * gap_extend).

Thus, (2 * match - 2 * substitution) must be greater than
(3 * match - 2 * gap_open - 2 * gap_extend), which simplifies to
(2 * gap_open + 2 * gap_extend > match + 2 * substitution).

There are two easy solutions to these inequalities:
- Bowtie v2.5 defaults: 6 > 5 - 2 and 2*5 + 2*3 > 2 + 2*6
- Set every score to 1: 1 > 1 - 1 and 2*1 + 2*1 > 1 + 2*1

"""

import re
from os import linesep
from pathlib import Path
from subprocess import CompletedProcess
from typing import Iterable

from .fqunit import FastqUnit, format_phred_arg
from ..core import path
from ..core.arg import (BOWTIE2_ORIENT,
                        TRIM_POLY_G_AUTO,
                        TRIM_POLY_G_NO,
                        TRIM_POLY_G_YES)
from ..core.extern import (BOWTIE2_CMD,
                           BOWTIE2_BUILD_CMD,
                           ECHO_CMD,
                           FASTP_CMD,
                           ShellCommand,
                           args_to_cmd,
                           cmds_to_pipe,
                           cmds_to_subshell)
from ..core.logs import logger
from ..core.ngs import (collate_xam_cmd,
                        run_flagstat,
                        sort_xam_cmd,
                        view_xam_cmd,
                        xam_to_fastq_cmd,
                        xam_paired,
                        FLAG_PAIRED,
                        FLAG_UNMAP,
                        FLAG_SECONDARY,
                        FLAG_QCFAIL,
                        FLAG_DUPLICATE,
                        FLAG_SUPPLEMENTARY)

# SAM filters
EXCLUDE_FLAGS = (FLAG_UNMAP
                 | FLAG_SECONDARY
                 | FLAG_QCFAIL
                 | FLAG_DUPLICATE
                 | FLAG_SUPPLEMENTARY)

# Bowtie2 parameters
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "1,1"
READ_GAP_PENALTY = "1,1"

# Fastp always outputs FASTQ files using Phred+33 encoding.
FASTP_PHRED_OUT = format_phred_arg(33)


def fastp_cmd(fq_inp: FastqUnit,
              fq_out: FastqUnit | None, *,
              fastp_html: Path,
              fastp_json: Path,
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
              fastp_adapter_fasta: Path | None,
              fastp_detect_adapter_for_pe: bool,
              fastp_min_length: int,
              n_procs: int):
    args = [FASTP_CMD,
            "--thread", n_procs,
            "--dont_eval_duplication",
            "--disable_quality_filtering"]
    # Length filter
    if fastp_min_length > 0:
        args.extend(["--length_required", fastp_min_length])
    else:
        args.append("--disable_length_filtering")
    # Quality trimming
    if fastp_5:
        args.append("--cut_front")
    if fastp_3:
        args.append("--cut_tail")
    if fastp_5 or fastp_3:
        args.extend(["--cut_window_size", fastp_w])
        args.extend(["--cut_mean_quality", fastp_m])
    # Poly(N) tail trimming
    if fastp_poly_g == TRIM_POLY_G_NO:
        args.append("--disable_trim_poly_g")
    else:
        if fastp_poly_g == TRIM_POLY_G_YES:
            args.append("--trim_poly_g")
        elif fastp_poly_g != TRIM_POLY_G_AUTO:
            raise ValueError(f"Invalid fastp_poly_g: {repr(fastp_poly_g)}")
        args.extend(["--poly_g_min_len", fastp_poly_g_min_len])
    if fastp_poly_x:
        args.extend(["--trim_poly_x", "--poly_x_min_len", fastp_poly_x_min_len])
    # Adapter trimming
    if fastp_adapter_trimming:
        if fastp_detect_adapter_for_pe and fq_inp.paired:
            args.append("--detect_adapter_for_pe")
        if fastp_adapter_1:
            args.extend(["--adapter_sequence", fastp_adapter_1])
        if fastp_adapter_2:
            if fq_inp.paired:
                args.extend(["--adapter_sequence_r2", fastp_adapter_2])
            else:
                logger.warning(
                    f"Ignored fastp_adapter_2 ({repr(fastp_adapter_2)}) "
                    f"for {fq_inp} because it has single-end reads"
                )
        if fastp_adapter_fasta:
            args.extend(["--adapter_fasta", fastp_adapter_fasta])
    else:
        args.append("--disable_adapter_trimming")
    # Input files
    if fq_inp.phred_enc == 64:
        args.append("--phred64")
    elif fq_inp.phred_enc != 33:
        raise ValueError("fastp requires a Phred encoding of +33 or +64, "
                         f"but got +{fq_inp.phred_enc}")
    if fq_inp.interleaved:
        args.append("--interleaved_in")
    for flag, value in zip(["-i", "-I"],
                           fq_inp.paths.values(),
                           strict=False):
        args.extend([flag, value])
    # Output files
    if fq_out is not None:
        for flag, value in zip(["-o", "-O"],
                               fq_out.paths.values(),
                               strict=False):
            args.extend([flag, value])
    else:
        # Output to stdout.
        args.append("--stdout")
    # Fastp report files
    args.extend(["--html", fastp_html])
    args.extend(["--json", fastp_json])
    return args_to_cmd(args)


run_fastp = ShellCommand("trimming base calls and filtering reads from",
                         fastp_cmd)


def get_bowtie2_index_paths(prefix: Path):
    """ Return the Bowtie 2 index paths for a FASTA file. """
    suffix = prefix.suffix
    if suffix in path.FASTA_EXTS:
        logger.warning(f"Bowtie 2 index prefix {prefix} has a FASTA extension")
    return [prefix.with_suffix(suffix + ext) for ext in path.BOWTIE2_INDEX_EXTS]


def bowtie2_build_cmd(fasta: Path, prefix: Path, *, n_procs: int = 1):
    """ Build a Bowtie2 index of a FASTA file. """
    # Generate and run the command. Use quiet mode (-q) to suppress the
    # default output, which is quite verbose.
    args = [BOWTIE2_BUILD_CMD, "-q", "--threads", n_procs, fasta, prefix]
    return args_to_cmd(args)


run_bowtie2_build = ShellCommand("building Bowtie 2 index for",
                                 bowtie2_build_cmd)


def _get_from_fq_inp(fq_inp: FastqUnit | None, attr: str):
    if fq_inp is None:
        raise TypeError(f"{attr} must be specified if fq_inp is None")
    return getattr(fq_inp, attr)


def bowtie2_cmd(fq_inp: FastqUnit | None,
                sam_out: Path | None, *,
                paired: bool | None = None,
                phred_arg: str | None = None,
                index_pfx: Path,
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
                fq_unal: Path | None = None,
                n_procs: int):
    if paired is None:
        paired = _get_from_fq_inp(fq_inp, "paired")
    if phred_arg is None:
        phred_arg = _get_from_fq_inp(fq_inp, "phred_arg")
    args = [BOWTIE2_CMD,
            # Resources
            "--threads", n_procs,
            # Alignment setup
            "--local" if bt2_local else "--end-to-end",
            "--gbar", bt2_gbar,
            "--dpad", bt2_dpad,
            "-L", bt2_l,
            "-i", bt2_s,
            "-D", bt2_d,
            "-R", bt2_r,
            # Scoring
            phred_arg,
            "--ignore-quals",
            "--ma", MATCH_BONUS if bt2_local else "0",
            "--mp", MISMATCH_PENALTY,
            "--np", N_PENALTY,
            "--rfg", REF_GAP_PENALTY,
            "--rdg", READ_GAP_PENALTY,
            # Filtering
            "--score-min", (bt2_score_min_loc if bt2_local
                            else bt2_score_min_e2e),
            "-I", bt2_i,
            "-X", bt2_x]
    # Mate pair orientation
    if bt2_orient not in BOWTIE2_ORIENT:
        raise ValueError(f"Invalid value for bt2_orient: {repr(bt2_orient)}")
    # Options for paired-end reads
    args.append(f"--{bt2_orient}")
    if not bt2_discordant:
        args.append("--no-discordant")
    if not bt2_contain:
        args.append("--no-contain")
    if bt2_dovetail:
        args.append("--dovetail")
    if not bt2_mixed:
        args.append("--no-mixed")
    # Filtering
    args.append("--no-unal")
    # Input and output files
    if fq_unal is not None:
        opts_unal = ["--un"]
        if paired:
            opts_unal.append("conc")
        if fq_unal.suffix.endswith(".gz"):
            opts_unal.append("gz")
        opt_unal = "-".join(opts_unal)
        args.extend([opt_unal, fq_unal])
    if sam_out is not None:
        args.extend(["-S", sam_out])
    args.extend(["-x", index_pfx])
    if fq_inp is not None:
        args.extend(fq_inp.bowtie2_inputs)
    else:
        # To provide paired-end reads on standard input, SEISMIC-RNA
        # only supports interleaved mates.
        args.extend(["--interleaved" if paired else "-U", "-"])
    return args_to_cmd(args)


def parse_bowtie2(process: CompletedProcess):
    """ Get the number of reads input and aligned. """
    indentation = 2
    pattern1 = re.compile(r"^(\s*)(\d+) ([ \w>]+); of these:$")
    pattern2 = re.compile(r"^(\s*)(\d+) \(\d+\.\d+%\) ([ \w>]+); of these:$")
    pattern3 = re.compile(r"^(\s*)(\d+) \(\d+\.\d+%\) ([ \w>]+)$")
    patterns = pattern1, pattern2, pattern3

    def parse_match(m: re.Match[str]):
        spaces, n_item, item = m.groups()
        return len(spaces) // indentation, item, int(n_item)

    n_reads = dict()
    names: tuple[str, ...] = tuple()
    lines = iter(process.stderr.split(linesep))
    # Read through the lines until one matches the first pattern.
    try:
        while not (match := pattern1.match(next(lines).rstrip())):
            pass
    except StopIteration:
        return n_reads
    # Read the remaining lines.
    while True:
        if match:
            # If the line matches, then find the name of the alignment
            # and the number of times it occurs.
            level, name, count = parse_match(match)
            # Key for dictionary.
            names = names[: level] + (name,)
            key = ", ".join(names)
            # Check if the key was encountered previously.
            if (prev := n_reads.get(key)) is None:
                # If not, then add the key.
                n_reads[key] = count
            elif prev != count:
                # If so, then confirm the count matches the previous.
                raise ValueError(
                    f"Inconsistent counts for {repr(key)}: {prev} ≠ {count}"
                )
        # Read the next line.
        try:
            line = next(lines).rstrip()
        except StopIteration:
            return n_reads
        # Try to match the line with each pattern, until one matches.
        match = None
        for pattern in patterns:
            if match := pattern.match(line):
                # The pattern matches.
                break


def xamgen_cmd(fq_inp: FastqUnit,
               bam_out: Path, *,
               fastp: bool,
               fastp_dir: Path,
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
               fastp_adapter_fasta: Path | None,
               fastp_detect_adapter_for_pe: bool,
               fastp_min_length: int,
               index_pfx: Path,
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
               fq_unal: Path | None = None,
               min_mapq: int | None = None,
               n_procs: int = 1):
    """ Wrap QC, alignment, and post-processing into one pipeline. """
    cmds = list()
    # Reserve one processor each for view and sort.
    n_procs_bowtie2 = max(n_procs - 2, 1)
    if fastp:
        # Allocate half of the remaining processors to Fastp.
        n_procs_fastp = max(n_procs_bowtie2 // 2, 1)
        n_procs_bowtie2 = max(n_procs_bowtie2 - n_procs_fastp, 1)
        if fq_inp.one_ref:
            fastp_html = f"{fq_inp.ref}__fastp.html"
            fastp_json = f"{fq_inp.ref}__fastp.json"
        else:
            fastp_html = "fastp.html"
            fastp_json = "fastp.json"
        cmds.append(fastp_cmd(
            fq_inp,
            None,
            fastp_html=fastp_dir.joinpath(fastp_html),
            fastp_json=fastp_dir.joinpath(fastp_json),
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
            fastp_adapter_fasta=fastp_adapter_fasta,
            fastp_detect_adapter_for_pe=fastp_detect_adapter_for_pe,
            fastp_min_length=fastp_min_length,
            n_procs=n_procs_fastp,
        ))
        # The input for Bowtie2 comes from what Fastp pipes out.
        bowtie2_fq_inp = None
        phred_arg = FASTP_PHRED_OUT
    else:
        # The input for Bowtie2 comes from the input FASTQ.
        bowtie2_fq_inp = fq_inp
        phred_arg = None
    cmds.append(bowtie2_cmd(
        bowtie2_fq_inp,
        None,
        paired=fq_inp.paired,
        phred_arg=phred_arg,
        index_pfx=index_pfx,
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
        fq_unal=fq_unal,
        n_procs=n_procs_bowtie2,
    ))
    # Filter out any unaligned or otherwise unsuitable reads.
    if fq_inp.paired:
        # Require the paired flag.
        flags_exc = EXCLUDE_FLAGS
        flags_req = FLAG_PAIRED
    else:
        # Exclude the paired flag and require no flags.
        flags_exc = EXCLUDE_FLAGS | FLAG_PAIRED
        flags_req = None
    cmds.append(view_xam_cmd(None,
                             None,
                             min_mapq=min_mapq,
                             flags_req=flags_req,
                             flags_exc=flags_exc,
                             bam=True))
    cmds.append(sort_xam_cmd(None, bam_out))
    return cmds_to_pipe(cmds)


run_xamgen = ShellCommand("aligning, filtering, and sorting by position",
                          xamgen_cmd,
                          parse_bowtie2)


def flags_cmds(xam_inp: Path,
               xam_out: Path | None, *,
               tmp_pfx: Path | None = None,
               flags_req: int | Iterable[int] = (),
               flags_exc: int | Iterable[int] = (),
               n_procs: int = 1):
    """ Filter a XAM file based on flags, then collate the output. """
    if not isinstance(xam_inp, Path):
        raise TypeError(f"Expected xam_inp to be a Path, "
                        f"but got {type(xam_inp).__name__}")
    # Pre-filter the reads for specific flags.
    if isinstance(flags_req, int):
        flags_req = [flags_req]
    elif not isinstance(flags_req, list):
        flags_req = list(flags_req)
    if isinstance(flags_exc, int):
        flags_exc = [flags_exc]
    elif not isinstance(flags_exc, list):
        flags_exc = list(flags_exc)
    num_flags = len(flags_req)
    if len(flags_exc) != num_flags:
        raise ValueError(f"Numbers of flags to require ({num_flags}) and "
                         f"exclude ({len(flags_exc)}) are different")
    if num_flags == 0:
        # View the XAM file with no flags.
        return [view_xam_cmd(xam_inp, xam_out, n_procs=n_procs)]
    n_procs_per_view = max(n_procs // (num_flags + 1), 1)
    multi_flags = num_flags > 1
    view_xam_cmds = [view_xam_cmd(xam_inp,
                                  None if multi_flags else xam_out,
                                  sam=multi_flags,
                                  with_header=(multi_flags and i == 0),
                                  flags_req=flags_req[i],
                                  flags_exc=flags_exc[i],
                                  n_procs=n_procs_per_view)
                     for i in range(num_flags)]
    if not multi_flags:
        # A single flag selection can be returned as-is.
        return view_xam_cmds
    # Reads from multiple flag selections must be collated.
    n_procs_collate = max(n_procs - n_procs_per_view * num_flags, 1)
    return [cmds_to_subshell(view_xam_cmds),
            collate_xam_cmd(None,
                            xam_out,
                            tmp_pfx=tmp_pfx,
                            fast=True,
                            n_procs=n_procs_collate)]


def flags_cmd(*args, **kwargs):
    """ Filter a XAM file based on flags, then collate the output. """
    return cmds_to_pipe(flags_cmds(*args, **kwargs))


run_flags = ShellCommand("filtering and collating", flags_cmd)


def realign_cmd(xam_inp: Path,
                xam_out: Path, *,
                paired: bool | None = None,
                tmp_pfx: Path | None = None,
                flags_req: int | Iterable[int] = (),
                flags_exc: int | Iterable[int] = (),
                min_mapq: int | None = None,
                n_procs: int = 1,
                **kwargs):
    """ Re-align reads that are already in a XAM file. """
    if paired is None:
        paired = xam_paired(run_flagstat(xam_inp, n_procs=n_procs))
    # Reserve processors for the flags commands, conversion to FASTQ,
    # and filtering.
    num_flags = (1 if isinstance(flags_req, int)
                 else len(flags_req := list(flags_req)))
    n_procs_bowtie2 = max(n_procs - max(num_flags, 1) - 2, 1)
    cmds = flags_cmds(xam_inp,
                      None,
                      tmp_pfx=tmp_pfx,
                      flags_req=flags_req,
                      flags_exc=flags_exc)
    # Convert the reads to FASTQ format.
    cmds.append(xam_to_fastq_cmd(None, None))
    # Re-align the reads from the BAM file using Bowtie2.
    cmds.append(bowtie2_cmd(None,
                            None,
                            paired=paired,
                            n_procs=n_procs_bowtie2,
                            **kwargs))
    # Filter low-quality alignments.
    cmds.append(view_xam_cmd(None, xam_out, min_mapq=min_mapq))
    return cmds_to_pipe(cmds)


run_realign = ShellCommand("filtering, collating, realigning, and formatting",
                           realign_cmd)


def export_cmd(xam_in: Path,
               xam_out: Path, *,
               ref: str,
               header: str,
               ref_file: Path | None = None,
               n_procs: int = 1):
    """ Select and export one reference to its own XAM file. """
    n_procs_per_view = max((n_procs - 1) // 2, 1)
    # Pipe the header line containing only this one reference.
    echo_step = args_to_cmd([ECHO_CMD, header])
    # Select only the reads that aligned to the reference; ignore the
    # original header so that the output XAM file contains only the
    # one reference in the XAM file.
    ref_step = view_xam_cmd(xam_in,
                            None,
                            sam=True,
                            with_header=False,
                            ref=ref,
                            n_procs=n_procs_per_view)
    # Merge the one header line and the reads for the reference.
    merge_step = cmds_to_subshell([echo_step, ref_step])
    # Export the reads into a XAM file.
    export_step = view_xam_cmd(None,
                               xam_out,
                               refs_file=ref_file,
                               n_procs=n_procs_per_view)
    return cmds_to_pipe([merge_step, export_step])


run_export = ShellCommand("selecting reference and exporting",
                          export_cmd)

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
