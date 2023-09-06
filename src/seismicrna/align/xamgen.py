"""

Alignment SAM/BAM/CRAM Generation Module

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
from logging import getLogger
from os import linesep
from pathlib import Path
from subprocess import CompletedProcess

from .fqunit import FastqUnit
from ..core import path
from ..core.cli import BOWTIE2_ORIENT
from ..core.shell import (BOWTIE2_CMD, BOWTIE2_BUILD_CMD, ECHO_CMD,
                          args_to_cmd, cmds_to_pipe, cmds_to_subshell,
                          ShellCommand)
from ..core.xam import (sort_xam_cmd, view_xam_cmd,
                        FLAG_PAIRED, FLAG_UNMAP, FLAG_SECONDARY, FLAG_QCFAIL,
                        FLAG_DUPLICATE, FLAG_SUPPLEMENTARY)

logger = getLogger(__name__)

# Bowtie2 parameters
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "1,1"
READ_GAP_PENALTY = "1,1"


def get_bowtie2_index_paths(prefix: Path):
    """ Return the Bowtie 2 index paths for a FASTA file. """
    suffix = prefix.suffix
    if suffix in path.FASTA_EXTS:
        logger.warning(f"Bowtie 2 index prefix {prefix} has a FASTA extension")
    return [prefix.with_suffix(suffix + ext) for ext in path.BOWTIE2_INDEX_EXTS]


def bowtie2_build_cmd(fasta: Path, prefix: Path, *, n_procs: int = 1):
    """ Build a Bowtie2 index of a FASTA file. """
    # Generate and run the command. Use quiet mode because otherwise,
    # Bowtie2-Build produces extremely verbose output.
    args = [BOWTIE2_BUILD_CMD, "-q", "--threads", n_procs, fasta, prefix]
    return args_to_cmd(args)


run_bowtie2_build = ShellCommand("building Bowtie 2 index for",
                                 bowtie2_build_cmd)


def bowtie2_cmd(fq_inp: FastqUnit,
                sam_out: Path | None, *,
                index_pfx: Path,
                n_procs: int,
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
                bt2_orient: str):
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
            fq_inp.phred_arg,
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
    if not bt2_unal:
        args.append("--no-unal")
    # Mate pair orientation
    if bt2_orient not in BOWTIE2_ORIENT:
        logger.warning(f"Invalid mate orientation for Bowtie2: '{bt2_orient}'. "
                       f"Setting to '{BOWTIE2_ORIENT[0]}'")
        bt2_orient = BOWTIE2_ORIENT[0]
    args.append(f"--{bt2_orient}")
    if not bt2_discordant:
        args.append("--no-discordant")
    if not bt2_contain:
        args.append("--no-contain")
    if bt2_dovetail:
        args.append("--dovetail")
    if not bt2_mixed:
        args.append("--no-mixed")
    # Formatting
    args.append("--xeq")
    # Input and output files
    if sam_out is not None:
        args.extend(["-S", sam_out])
    args.extend(["-x", index_pfx])
    args.extend(fq_inp.bowtie2_inputs)
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
    lines = iter(process.stderr.decode().split(linesep))
    # Read through the lines until one matches the first pattern.
    while not (match := pattern1.match(line := next(lines, "").rstrip())):
        if not line:
            # Prevent an infinite loop if lines becomes exhausted.
            return n_reads
    # Read the rest of the lines until the iterator is exhausted.
    while line:
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
                    f"Inconsistent counts for '{key}': {prev} ≠ {count}")
        # Read the next line, defaulting to an empty string.
        line = next(lines, "").rstrip()
        # Try to match the line with each pattern, until one matches.
        match = None
        for pattern in patterns:
            if match := pattern.match(line):
                # The pattern matches.
                break
    return n_reads


run_bowtie2 = ShellCommand("aligning", bowtie2_cmd, parse_bowtie2)


def xamgen_cmd(fq_inp: FastqUnit,
               bam_out: Path | None, *,
               min_mapq: int,
               n_procs: int,
               **kwargs):
    """ Wrap alignment and post-processing into one pipeline. """
    bowtie2_step = bowtie2_cmd(fq_inp, None, n_procs=n_procs, **kwargs)
    # Filter out any unaligned or otherwise unsuitable reads.
    flags_exc = (FLAG_UNMAP | FLAG_SECONDARY | FLAG_QCFAIL
                 | FLAG_DUPLICATE | FLAG_SUPPLEMENTARY)
    # Ensure all output reads have the right pairing status.
    if fq_inp.paired:
        # Require the paired flag.
        flags_req = FLAG_PAIRED
    else:
        # Exclude the paired flag and require no flags.
        flags_exc |= FLAG_PAIRED
        flags_req = None
    view_xam_step = view_xam_cmd(None, None,
                                 min_mapq=min_mapq,
                                 flags_req=flags_req,
                                 flags_exc=flags_exc,
                                 bam=True,
                                 n_procs=n_procs)
    sort_xam_step = sort_xam_cmd(None, bam_out, n_procs=n_procs)
    return cmds_to_pipe([bowtie2_step, view_xam_step, sort_xam_step])


run_xamgen = ShellCommand("aligning, filtering, and sorting by position",
                          xamgen_cmd,
                          parse_bowtie2)


def export_cmd(xam_in: Path | None,
               cram_out: Path | None, *,
               ref: str,
               header: str,
               ref_file: Path,
               n_procs: int):
    """ Wrap sorting and cramming into one pipeline. """
    # Pipe the header line.
    echo_step = args_to_cmd([ECHO_CMD, header])
    # Select only the reads that aligned to the reference, and ignore
    # the original header.
    ref_step = view_xam_cmd(xam_in, None,
                            sam=True, with_header=False,
                            ref=ref, n_procs=n_procs)
    # Merge the one header line and the reads for the reference.
    merge_step = cmds_to_subshell([echo_step, ref_step])
    # Sort reads by name so that mates are adjacent.
    sort_step = sort_xam_cmd(None, None, name=True, n_procs=n_procs)
    # Cram the reads into a CRAM file.
    cram_step = view_xam_cmd(None, cram_out, refs_file=ref_file,
                             n_procs=n_procs)
    return cmds_to_pipe([merge_step, sort_step, cram_step])


run_export = ShellCommand("selecting reference, sorting by name, and cramming",
                          export_cmd)
