import re
from logging import getLogger
from pathlib import Path
from subprocess import CompletedProcess

from . import path
from .shell import PipelineStep, ParsedPipelineStep, make_cmd, SAMTOOLS_CMD

logger = getLogger(__name__)

# SAM file format specifications
SAM_HEADER = '@'
SAM_DELIM = '\t'
SAM_ALIGN_SCORE = 'AS:i:'
SAM_EXTRA_SCORE = 'XS:i:'
FLAG_PAIRED = 2 ** 0
FLAG_PROPER = 2 ** 1
FLAG_UNMAP = 2 ** 2
FLAG_MUNMAP = 2 ** 3
FLAG_REVERSE = 2 ** 4
FLAG_MREVERSE = 2 ** 5
FLAG_FIRST = 2 ** 6
FLAG_SECOND = 2 ** 7
FLAG_SECONDARY = 2 ** 8
FLAG_QCFAIL = 2 ** 9
FLAG_DUPLICATE = 2 ** 10
FLAG_SUPPLEMENTARY = 2 ** 11
MAX_FLAG = sum([FLAG_PAIRED, FLAG_PROPER,
                FLAG_UNMAP, FLAG_MUNMAP,
                FLAG_REVERSE, FLAG_MREVERSE,
                FLAG_FIRST, FLAG_SECOND,
                FLAG_SECONDARY, FLAG_QCFAIL,
                FLAG_DUPLICATE, FLAG_SUPPLEMENTARY])


def get_bam_index(bam: Path):
    """ Get the path to an index of a BAM file. """
    return bam.with_suffix(path.BAI_EXT)


def get_args_index_bam(bam: Path, _: None = None, *, n_procs: int = 1):
    """ Build an index of a BAM file using `samtools index`. """
    return [SAMTOOLS_CMD, "index", "-@", n_procs - 1, bam]


index_bam = PipelineStep(get_args_index_bam, make_cmd, "indexing")


def get_args_sort_xam(xam_inp: Path | None,
                      xam_out: Path | None, *,
                      name: bool = False,
                      n_procs: int = 1):
    """ Sort a SAM or BAM file using `samtools sort`. """
    args = [SAMTOOLS_CMD, "sort", "-@", n_procs - 1]
    if name:
        # Sort by name instead of coordinate.
        args.append("-n")
    if xam_out:
        args.extend(["-o", xam_out])
    else:
        # To increase speed, do not compress on stdout.
        args.append("-u")
    if xam_inp:
        args.append(xam_inp)
    return args


sort_xam = PipelineStep(get_args_sort_xam, make_cmd, "sorting")


def get_args_view_xam(xam_inp: Path | None,
                      xam_out: Path | None, *,
                      bam: bool = False,
                      cram: bool = False,
                      header: bool = False,
                      min_mapq: int | None = None,
                      flags_req: int | None = None,
                      flags_exc: int | None = None,
                      ref: str | None = None,
                      end5: int | None = None,
                      end3: int | None = None,
                      refs_file: Path | None = None,
                      n_procs: int = 1):
    """ Convert between SAM and BAM formats, extract reads aligning to a
    specific reference/section, and filter by flag and mapping quality
    using `samtools view`. """
    args = [SAMTOOLS_CMD, "view", "-@", n_procs - 1]
    # Read filters
    if min_mapq is not None:
        # Require minimum mapping quality.
        args.extend(["--min-mq", min_mapq])
    if flags_req is not None:
        # Require these flags.
        args.extend(["-f", flags_req])
    if flags_exc is not None:
        # Exclude these flags.
        args.extend(["-F", flags_exc])
    # Output format
    if cram:
        args.append("--cram")
        if bam:
            logger.warning("Both BAM and CRAM flags were set: using CRAM")
        if not refs_file:
            logger.warning("Missing reference file for CRAM output")
    elif bam:
        args.append("--bam")
    else:
        # The header argument affects only output in SAM format.
        args.append("--with-header" if header else "--no-header")
    # Reference file
    if refs_file:
        args.extend(["--reference", refs_file])
    # Input and output files
    if xam_out:
        args.extend(["-o", xam_out])
    else:
        # To increase speed, do not compress on stdout.
        args.append("-u")
    if xam_inp:
        args.append(xam_inp)
    # Reference and section specification
    if ref is not None:
        if end5 is not None and end3 is not None:
            # View only reads aligning to a section of this reference.
            args.append(f"{ref}:{end5}-{end3}")
        else:
            # View only reads aligning to this reference.
            args.append(ref)
            if end5 is not None:
                logger.warning(f"Got end5 = {end5} but not end3")
            if end3 is not None:
                logger.warning(f"Got end3 = {end3} but not end5")
    elif end5 is not None or end5 is not None:
        logger.warning(f"Options end5 and end3 require a reference name")
    return args


view_xam = PipelineStep(get_args_view_xam, make_cmd, "viewing")


def get_args_flagstat(xam_inp: Path | None,
                      _: None = None, *,
                      n_procs: int = 1):
    """ Compute the statistics with `samtools flagstat`. """
    args = [SAMTOOLS_CMD, "flagstat", "-@", n_procs - 1]
    if xam_inp:
        args.append(xam_inp)
    return args


def parse_stdout_flagstat(process: CompletedProcess):
    """ Convert the output into a dict with one entry per line. """
    stdout = process.stdout.decode()
    stats_pattern = "([0-9]+) [+] ([0-9]+) ([A-Za-z0-9 ]+)"
    return {stat.strip(): (int(n1), int(n2))
            for n1, n2, stat in map(re.Match.groups,
                                    re.finditer(stats_pattern, stdout))}


run_flagstat = ParsedPipelineStep(get_args_flagstat,
                                  make_cmd,
                                  parse_stdout_flagstat,
                                  "computing flagstats")


def count_single_paired(flagstats: dict):
    """ Count the records in a SAM/BAM file given an output dict from
    `get_flagstats()`. """
    mapped, _ = flagstats["mapped"]
    # Count paired-end reads with both mates in the file.
    paired_self_mate, _ = flagstats["with itself and mate mapped"]
    paired_two, extra = divmod(paired_self_mate, 2)
    if extra:
        raise ValueError(f"Reads with themselves and their mates mapped must "
                         f"be even, but got {paired_self_mate}")
    # Count paired-end reads with only one mate in the file.
    paired_one, _ = flagstats["singletons"]
    # Count single-end reads.
    singles = mapped - (paired_one + paired_self_mate)
    logger.debug(f"Single-end: {singles}\n"
                 f"Paired-end, one mate: {paired_one}\n"
                 f"Paired-end, two mates: {paired_two}")
    return paired_two, paired_one, singles


def count_total_records(flagstats: dict):
    """ Count the unique records in a SAM/BAM file. """
    return sum(count_single_paired(flagstats))
