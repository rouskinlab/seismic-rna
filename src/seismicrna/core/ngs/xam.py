import re
from logging import getLogger
from pathlib import Path
from subprocess import CompletedProcess

from ..extern import SAMTOOLS_CMD, args_to_cmd, ShellCommand

logger = getLogger(__name__)

# SAM file format specifications
SAM_HEADER = "@"
SAM_DELIM = "\t"
SAM_NOREF = "*"
SAM_SEQLINE = "@SQ"
SAM_SEQNAME = "SN:"
SAM_SEQLEN = "LN:"
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
MAX_FLAG = sum([FLAG_PAIRED,
                FLAG_PROPER,
                FLAG_UNMAP,
                FLAG_MUNMAP,
                FLAG_REVERSE,
                FLAG_MREVERSE,
                FLAG_FIRST,
                FLAG_SECOND,
                FLAG_SECONDARY,
                FLAG_QCFAIL,
                FLAG_DUPLICATE,
                FLAG_SUPPLEMENTARY])


def calc_extra_threads(n_procs: int):
    """ Calculate the number of extra threads to use (option -@). """
    try:
        if not isinstance(n_procs, int):
            raise TypeError(
                f"n_procs must be int, but got {type(n_procs).__name__}"
            )
        if n_procs < 1:
            raise ValueError(f"n_procs must be ≥ 1, but got {n_procs}")
        return n_procs - 1
    except (TypeError, ValueError) as error:
        logger.warning(error)
        return 0


def index_xam_cmd(bam: Path, *, n_procs: int = 1):
    """ Build an index of a XAM file using `samtools index`. """
    return args_to_cmd([SAMTOOLS_CMD, "index",
                        "-@", calc_extra_threads(n_procs),
                        bam])


run_index_xam = ShellCommand("indexing alignment map",
                             index_xam_cmd,
                             opath=False)


def sort_xam_cmd(xam_inp: Path | None,
                 xam_out: Path | None, *,
                 tmp_pfx: Path | None = None,
                 name: bool = False,
                 n_procs: int = 1):
    """ Sort a SAM or BAM file using `samtools sort`. """
    args = [SAMTOOLS_CMD, "sort",
            "-@", calc_extra_threads(n_procs)]
    if name:
        # Sort by name instead of coordinate.
        args.append("-n")
    # Write temporary files with this prefix.
    if tmp_pfx is None and xam_out is not None:
        tmp_pfx = xam_out.with_suffix("")
    if tmp_pfx is not None:
        args.extend(["-T", tmp_pfx])
    if xam_out:
        args.extend(["-o", xam_out])
    else:
        # To increase speed, do not compress on stdout.
        args.append("-u")
    if xam_inp:
        args.append(xam_inp)
    return args_to_cmd(args)


run_sort_xam = ShellCommand("sorting alignment map", sort_xam_cmd)


def collate_xam_cmd(xam_inp: Path | None,
                    xam_out: Path | None, *,
                    tmp_pfx: Path | None = None,
                    fast: bool = False,
                    n_procs: int = 1):
    """ Collate a SAM or BAM file using `samtools collate`. """
    args = [SAMTOOLS_CMD, "collate",
            "-@", calc_extra_threads(n_procs)]
    if fast:
        # Use fast mode (outputs primary alignments only).
        args.append("-f")
    # Write temporary files with this prefix.
    if tmp_pfx is None and xam_out is not None:
        tmp_pfx = xam_out.with_suffix("")
    if tmp_pfx is not None:
        args.extend(["-T", tmp_pfx])
    if xam_out:
        args.extend(["-o", xam_out])
    else:
        # Output to stdout.
        args.append("-O")
        # To increase speed, do not compress on stdout.
        args.append("-u")
    # samtools collate requires an input argument; - means to read from
    # standard input.
    args.append(xam_inp if xam_inp else "-")
    return args_to_cmd(args)


def view_xam_cmd(xam_inp: Path | None,
                 xam_out: Path | None, *,
                 sam: bool = False,
                 bam: bool = False,
                 cram: bool = False,
                 with_header: bool = False,
                 only_header: bool = False,
                 min_mapq: int = 0,
                 flags_req: int = 0,
                 flags_exc: int = 0,
                 ref: str | None = None,
                 end5: int | None = None,
                 end3: int | None = None,
                 refs_file: Path | None = None,
                 n_procs: int = 1):
    """ Convert between SAM and BAM formats, extract reads aligning to a
    specific reference/section, and filter by flag and mapping quality
    using `samtools view`. """
    args = [SAMTOOLS_CMD, "view",
            "-@", calc_extra_threads(n_procs)]
    # Read filters
    if min_mapq:
        # Require minimum mapping quality.
        args.extend(["-q", min_mapq])
    if flags_req:
        # Require these flags.
        args.extend(["-f", flags_req])
    if flags_exc:
        # Exclude these flags.
        args.extend(["-F", flags_exc])
    # Header
    if with_header:
        args.append("-h")
    if only_header:
        args.append("-H")
    # Reference file
    if refs_file:
        args.extend(["-T", refs_file])
    # Output format
    if cram:
        args.append("-C")
        if sam:
            logger.warning("Both SAM and CRAM flags were set: using CRAM")
        if bam:
            logger.warning("Both BAM and CRAM flags were set: using CRAM")
        if not refs_file:
            logger.warning("Missing reference file for CRAM output")
    elif bam:
        args.append("-b")
        if sam:
            logger.warning("Both BAM and SAM flags were set: using BAM")
    # Input and output files
    if xam_out:
        args.extend(["-o", xam_out])
    elif not sam:
        # To increase speed, do not compress binary standard output.
        args.append("-u")
    if xam_inp:
        args.append(xam_inp)
    # Reference and section specification
    if ref:
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
    return args_to_cmd(args)


def flagstat_cmd(xam_inp: Path | None, *, n_procs: int = 1):
    """ Compute the statistics with `samtools flagstat`. """
    args = [SAMTOOLS_CMD, "flagstat",
            "-@", calc_extra_threads(n_procs)]
    if xam_inp:
        args.append(xam_inp)
    return args_to_cmd(args)


def parse_flagstat(process: CompletedProcess):
    """ Convert the output into a dict with one entry per line. """
    stats_pattern = "([0-9]+) [+] ([0-9]+) ([A-Za-z0-9 ]+)"
    return {stat.strip(): (int(n1), int(n2))
            for n1, n2, stat in map(re.Match.groups,
                                    re.finditer(stats_pattern, process.stdout))}


run_flagstat = ShellCommand("computing flagstats",
                            flagstat_cmd,
                            parse_flagstat,
                            opath=False)


def count_single_paired(flagstats: dict):
    """ Count the records in a SAM/BAM file given an output dict from
    `get_flagstats()`. """
    mapped, _ = flagstats["primary mapped"]
    # Count properly paired reads.
    paired_end_reads_proper, _ = flagstats["properly paired"]
    paired_end_pairs_proper, extra = divmod(paired_end_reads_proper, 2)
    if extra:
        raise ValueError("Number of properly paired reads must be even, "
                         f"but got {paired_end_reads_proper}")
    # Count improper paired-end reads (either only one mate exists or
    # both mates exist but did not pair properly with each other).
    paired_end_reads_both_mates, _ = flagstats["with itself and mate mapped"]
    paired_end_reads_both_mates_improper = (paired_end_reads_both_mates
                                            - paired_end_reads_proper)
    if paired_end_reads_both_mates_improper < 0:
        raise ValueError(
            "Number of paired-end reads with both mates mapped and paired "
            "properly with each other must be no more than the number of "
            "paired-end reads with both mates mapped (though not necessarily "
            f"with each other), but got {paired_end_reads_proper} "
            f"> {paired_end_reads_both_mates}, which indicates a bug"
        )
    paired_end_reads_one_mate, _ = flagstats["singletons"]
    paired_end_reads_improper = (paired_end_reads_both_mates_improper
                                 + paired_end_reads_one_mate)
    # Count single-end reads.
    paired_end_reads = paired_end_reads_proper + paired_end_reads_improper
    single_end_reads = mapped - paired_end_reads
    if single_end_reads < 0:
        raise ValueError(
            "Number of paired-end reads must be no more than the number of "
            f"mapped reads in total, but got {paired_end_reads} > {mapped}, "
            "which indicates a bug"
        )
    logger.debug(f"Proper pairs: {paired_end_pairs_proper}\n"
                 f"Other paired-end reads: {paired_end_reads_improper}\n"
                 f"Single-end reads: {single_end_reads}")
    return paired_end_pairs_proper, paired_end_reads_improper, single_end_reads


def count_total_reads(flagstats: dict):
    """ Count the total records in a SAM/BAM file. """
    return sum(count_single_paired(flagstats))


def xam_paired(flagstats: dict):
    """ Determine if the reads are single-end or paired-end. """
    # Determine if there are any paired-end and single-end reads.
    paired_two, paired_one, singles = count_single_paired(flagstats)
    paired = paired_two > 0 or paired_one > 0
    single = singles > 0
    # SEISMIC-RNA currently cannot handle both types in one sample.
    if paired and single:
        raise ValueError(f"Got both single-end (n = {singles}) and paired-end "
                         f"(n = {paired_one + paired_two}) reads")
    # The pairing status cannot be determined if there are no reads.
    if not paired and not single:
        logger.warning("Got 0 reads: neither single-end nor paired-end")
    # Return True if definitely paired-end, False if single-end or if
    # there are no reads.
    return paired


def idxstats_cmd(xam_inp: Path):
    """ Count the number of reads aligning to each reference. """
    return args_to_cmd([SAMTOOLS_CMD, "idxstats", xam_inp])


def parse_idxstats(process: CompletedProcess):
    """ Map each reference to the number of reads aligning to it. """
    counts = dict()
    for line in process.stdout.splitlines():
        if stripped_line := line.rstrip():
            ref, length, mapped, unmapped = stripped_line.split(SAM_DELIM)
            if ref != SAM_NOREF:
                counts[ref] = int(mapped)
    # Sort the references in from most to least abundant.
    return {ref: counts[ref]
            for ref in sorted(counts, key=counts.__getitem__, reverse=True)}


run_idxstats = ShellCommand("counting reads for each reference",
                            idxstats_cmd,
                            parse_idxstats,
                            opath=False)


def ref_header_cmd(xam_inp: Path, *, n_procs: int):
    """ Get the header line for each reference. """
    return view_xam_cmd(xam_inp,
                        None,
                        sam=True,
                        only_header=True,
                        n_procs=n_procs)


def parse_ref_header(process: CompletedProcess):
    """ Map each reference to its header line. """
    for line in process.stdout.splitlines():
        if line.startswith(SAM_SEQLINE):
            # Find the field that has the reference name.
            for field in line.split(SAM_DELIM):
                if field.startswith(SAM_SEQNAME):
                    # Yield the reference name and the full line.
                    ref = field[len(SAM_SEQNAME):]
                    yield ref, line.rstrip()
                    break
            else:
                # The sequence name was not found.
                logger.warning(f"Failed to find sequence name in line:\n{line}")


run_ref_header = ShellCommand("getting header line for each reference",
                              ref_header_cmd,
                              parse_ref_header,
                              opath=False)


def xam_to_fastq_cmd(xam_inp: Path | None,
                     fq_out: Path | None, *,
                     flags_req: int | None = None,
                     flags_exc: int | None = None,
                     label_12: bool = False,
                     n_procs: int = 1):
    """ Convert XAM format to FASTQ format, and filter by flags. """
    args = [SAMTOOLS_CMD, "fastq",
            "-@", calc_extra_threads(n_procs)]
    if flags_req is not None:
        # Require these flags.
        args.extend(["-f", flags_req])
    if flags_exc is not None:
        # Exclude these flags.
        args.extend(["-F", flags_exc])
    if not label_12:
        # Do not add /1 or /2 labels to names of paired-end reads.
        args.append("-n")
    args.append(xam_inp if xam_inp is not None else "-")
    if fq_out is not None:
        args.extend([">", fq_out])
    return args_to_cmd(args)

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
