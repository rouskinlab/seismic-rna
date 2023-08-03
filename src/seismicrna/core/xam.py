from logging import getLogger
from pathlib import Path
from math import isnan, nan
import re
from typing import BinaryIO

from . import path
from .shell import run_cmd, SAMTOOLS_CMD

logger = getLogger(__name__)

# SAM file format specifications
SAM_HEADER = b"@"
SAM_DELIMITER = b"\t"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
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


def index_bam(bam: Path, n_procs: int = 1):
    """ Build an index of a BAM file using `samtools index`. """
    logger.info(f"Began building BAM index of {bam}")
    index = bam.with_suffix(path.BAI_EXT)
    cmd = [SAMTOOLS_CMD, "index", "-@", n_procs - 1, bam]
    # Build the BAM index.
    run_cmd(cmd, check_created=[index])
    logger.info(f"Ended building BAM index of {bam}: "
                f"{bam.with_suffix(path.BAI_EXT)}")
    return index


def sort_xam(xam_inp: Path, xam_out: Path, *,
             name: bool = False, n_procs: int = 1):
    """ Sort a SAM or BAM file using `samtools sort`. """
    logger.info(f"Began sorting {xam_inp}")
    cmd = [SAMTOOLS_CMD, "sort", "-@", n_procs - 1]
    if name:
        # Sort by name instead of coordinate.
        cmd.append("-n")
    cmd.extend(["-o", xam_out, xam_inp])
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created directory: {xam_out.parent}")
    # Sort the SAM/BAM file.
    run_cmd(cmd, check_created=[xam_out])
    logger.info(f"Ended sorting {xam_inp} to {xam_out}")


def view_xam(xam_inp: Path,
             xam_out: Path, *,
             ref: str | None = None,
             end5: int | None = None,
             end3: int | None = None,
             flags_req: int | None = None,
             flags_exc: int | None = None,
             n_procs: int = 1):
    """ Convert between SAM and BAM formats, or extract reads aligning
    to a specific reference/section using `samtools view`. """
    logger.info(f"Began viewing {xam_inp}")
    cmd = [SAMTOOLS_CMD, "view", "-@", n_procs - 1, "-h"]
    if flags_req is not None:
        # Require these flags.
        cmd.extend(["-f", flags_req])
    if flags_exc is not None:
        # Exclude these flags.
        cmd.extend(["-F", flags_exc])
    if xam_out.suffix == path.BAM_EXT:
        # Write a binary (BAM) file.
        cmd.append("-b")
    # Input and output files
    cmd.extend(["-o", xam_out, xam_inp])
    # Reference and section specification
    if ref is not None:
        if end5 is not None and end3 is not None:
            # View only reads aligning to a section of this reference.
            cmd.append(f"{ref}:{end5}-{end3}")
        else:
            # View only reads aligning to this reference.
            cmd.append(ref)
            if end5 is not None:
                logger.warning(f"Got end5 = {end5} but not end3")
            if end3 is not None:
                logger.warning(f"Got end3 = {end3} but not end5")
    elif end5 is not None or end5 is not None:
        logger.warning(f"Options end5 and end3 require a reference name")
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created directory: {xam_out.parent}")
    # View the SAM/BAM file.
    run_cmd(cmd, check_created=[xam_out])
    logger.info(f"Ended viewing {xam_inp} as {xam_out}")


def get_flagstat(xam_inp: Path, n_procs: int = 1):
    """ Return a dict of the output from `samtools flagstat`. """
    logger.info(f"Began computing flag statistics on {xam_inp}")
    # Use samtools flagstat to compute the flag statistics.
    cmd = [SAMTOOLS_CMD, "flagstat", "-@", n_procs - 1, xam_inp]
    output = run_cmd(cmd, check_is_before=[xam_inp]).stdout
    # Convert the output into a dict with one entry per line.
    stat_pattern = b"([0-9]+) [+] ([0-9]+) ([A-Za-z0-9 ]+)"
    stats = {stat.strip(): (int(n1), int(n2))
             for n1, n2, stat in map(re.Match.groups,
                                     re.finditer(stat_pattern, output))}
    logger.debug(f"Flag statistics for {xam_inp}:\n{stats}")
    logger.info(f"Ended computing flag statistics on {xam_inp}")
    return stats


def count_xam(xam_inp: Path, n_procs: int = 1):
    """ Count the unique records in a SAM/BAM file. """
    logger.info(f"Began counting records in {xam_inp}")
    # Compute the flag statistics.
    flagstat = get_flagstat(xam_inp, n_procs)
    mapped, _ = flagstat[b"mapped"]
    # Count paired-end reads with both mates in the file.
    paired_two, _ = flagstat[b"with itself and mate mapped"]
    # Count paired-end reads with only one mate in the file.
    paired_one, _ = flagstat[b"singletons"]
    # Count single-end reads.
    singles = mapped - (paired_one + paired_two)
    # Count unique reads (two mates of one pair count as one).
    unique = singles + paired_one + paired_two // 2
    logger.debug(f"Number of records in {xam_inp}:\n"
                 f"Single-end: {singles}\n"
                 f"Paired-end, one mate: {paired_one}\n"
                 f"Paired-end, two mates: {paired_two // 2}")
    logger.info(f"Ended counting {unique} records in {xam_inp}")
    return unique


def dedup_sam(sam_inp: Path, sam_out: Path):
    """ Remove SAM reads that map equally to multiple locations. """

    logger.info(f"Began deduplicating {sam_inp}")

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"([-0-9]+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"([-0-9]+)")

    min_fields = 11
    max_flag = 4095  # 2^12 - 1

    def get_score(line: bytes, ptn: re.Pattern[bytes]):
        """ Get the alignment score from a line in a SAM file. """
        return float(match.groups()[0]) if (match := ptn.search(line)) else nan

    def is_best_align(line: bytes):
        """ Return whether the line contains the best alignment of the
        read it contains. """
        try:
            if isnan(score_a := get_score(line, pattern_a)):
                logger.warning(f"Missing alignment score for {line.decode()}")
                return True
            # Compare using "not >=" instead of "<" because, if score x
            # is missing (NaN), the comparison will be False, but this
            # function should return False only if score x exists and
            # is greater than or equal to score_a.
            return not get_score(line, pattern_x) >= score_a
        except Exception as error:
            raise ValueError(f"Failed to determine if line {line} in {sam_inp} "
                             f"is the best alignment: {error}")

    def read_is_paired(line: bytes):
        info = line.split()
        if len(info) < min_fields:
            raise ValueError(f"Invalid line in {sam_inp}:\n{line.decode()}")
        flag = int(info[1])
        if flag < 0 or flag > max_flag:
            raise ValueError(f"Invalid flag in {sam_inp}: {flag}")
        return bool(flag & FLAG_PAIRED)

    def write_summary_single(written: int, skipped: int, errors: int):
        total = written + skipped + errors
        try:
            fw = 100 * written / total
            fs = 100 * written / total
            fe = 100 * errors / total
        except ZeroDivisionError:
            fw, fs, fe = float("nan"), float("nan"), float("nan")
        return (f"\n\nSummary of deduplicating {sam_inp}\n"
                f"Total reads: {total:>12}\n"
                f"Uniquely-mapped (written): {written:>12} ({fw:>6.2f}%)\n"
                f"Multiply-mapped (skipped): {skipped:>12} ({fs:>6.2f}%)\n"
                f"Reads w/ errors (skipped): {errors:>12} ({fe:>6.2f}%)\n\n")

    def write_summary_paired(pwrit: int, pskip: int, perro: int,
                             mwrit: int, mskip: int, merro: int):
        mates = mwrit + mskip + merro
        pairs = pwrit + pskip + perro
        pairs1 = 2 * pairs - mates
        pairs2 = pairs - pairs1
        try:
            f1 = 100 * pairs1 / pairs
            f2 = 100 * pairs2 / pairs
            fw = 100 * pwrit / pairs
            fs = 100 * pskip / pairs
            fe = 100 * perro / pairs
        except ZeroDivisionError:
            f1, f2, fw, fs, fe = nan, nan, nan, nan, nan
        return (f"\n\nSummary of deduplicating {sam_inp}\n"
                f"Total mates: {mates:>12}\n"
                f"Total pairs: {pairs:>12}\n"
                f"Pairs w/ both mates mapped: {pairs2:>12} ({f2:>6.2f}%)\n"
                f"Pairs w/ one mate unmapped: {pairs1:>12} ({f1:>6.2f}%)\n"
                f"Uniquely-mapped (written):  {pwrit:>12} ({fw:>6.2f}%)\n"
                f"Multiply-mapped (skipped):  {pskip:>12} ({fs:>6.2f}%)\n"
                f"Pairs w/ errors (skipped):  {perro:>12} ({fe:>6.2f}%)\n\n")

    def iter_single(sam: BinaryIO, line: bytes):
        """ For each read, yield the best-scoring alignment, excluding
        reads that aligned equally well to multiple locations. """
        n_copy = 0
        n_skip = 0
        n_erro = 0
        while line:
            try:
                if is_best_align(line):
                    n_copy += 1
                    yield line
                else:
                    n_skip += 1
            except Exception as error:
                n_erro += 1
                logger.error(error)
            line = sam.readline()
        logger.info(write_summary_single(n_copy, n_skip, n_erro))

    def iter_paired(sam: BinaryIO, line1: bytes):
        """ For each pair of reads, yield the pair of alignments for
        which both the forward alignment and the reverse alignment in
        the pair scored best among all alignments for the forward and
        reverse reads, respectively. Exclude pairs for which the forward
        and/or reverse read aligned equally well to multiple locations,
        or for which the best alignments for the forward and reverse
        reads individually are not part of the same alignment pair. """
        n_pairs_written = 0
        n_pairs_skipped = 0
        n_pairs_errors = 0
        n_mates_written = 0
        n_mates_skipped = 0
        n_mates_errors = 0
        while line1:
            if SAM_DELIMITER not in line1:
                # If the line is blank, skip it.
                line1 = sam.readline()
                continue
            # Check if there is at least one more line in the file.
            if (line2 := sam.readline()) and SAM_DELIMITER in line2:
                # Check if the reads on lines 1 and 2 are mates.
                if (line1.split(SAM_DELIMITER, 1)[0]
                        == line2.split(SAM_DELIMITER, 1)[0]):
                    # If they are mates of each other, check if the
                    # lines represent the best alignment for both.
                    try:
                        if is_best_align(line1) and is_best_align(line2):
                            # If so, then yield both mates.
                            n_pairs_written += 1
                            n_mates_written += 2
                            yield line1
                            yield line2
                        else:
                            # Otherwise, skip both mates.
                            n_pairs_skipped += 1
                            n_mates_skipped += 2
                    except Exception as error:
                        n_pairs_errors += 1
                        n_mates_errors += 2
                        logger.error(error)
                    # Read the next line from the file.
                    line1 = sam.readline()
                    continue
            # If line2 does not have the mate of line1, then check only
            # line1 now and consider line2 on the next iteration.
            try:
                if is_best_align(line1):
                    # Yield only line1.
                    n_pairs_written += 1
                    n_mates_written += 1
                    yield line1
                else:
                    # Skip line1.
                    n_pairs_skipped += 1
                    n_mates_skipped += 1
            except Exception as error:
                n_pairs_errors += 1
                n_mates_errors += 1
                logger.error(error)
            # Since line2 was not used yet, consider it on the next
            # iteration instead of reading another line from the file.
            line1 = line2
        logger.info(write_summary_paired(n_pairs_written,
                                         n_pairs_skipped,
                                         n_pairs_errors,
                                         n_mates_written,
                                         n_mates_skipped,
                                         n_mates_errors))

    # Make the output directory.
    sam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Creating directory: {sam_out.parent}")

    # Deduplicate the alignments.
    with (open(sam_inp, "rb") as sami, open(sam_out, "xb") as samo):
        # Copy the entire header from the input to the output SAM file.
        n_header = 0
        while (line_ := sami.readline()).startswith(SAM_HEADER):
            samo.write(line_)
            n_header += 1
        logger.debug(
            f"Copied {n_header} header lines from {sam_inp} to {sam_out}")
        # Copy only the reads that mapped best to one location from the
        # input to the output SAM file.
        if line_:
            is_paired = read_is_paired(line_)
            # Determine whether the reads in the SAM file are paired.
            if is_paired:
                endedness = "paired"
                lines = iter_paired(sami, line_)
            else:
                endedness = "single"
                lines = iter_single(sami, line_)
            logger.debug(f"Deduplicating {endedness}-end SAM file: {sam_inp}")
            # Write the selected lines to the output SAM file.
            for line_ in lines:
                samo.write(line_)
        else:
            logger.critical(f"SAM file {sam_inp} contained no reads")

    if not sam_out.is_file():
        raise FileNotFoundError(f"Failed to deduplicate {sam_inp} to {sam_out}")
    logger.info(f"Ended deduplicating {sam_inp} to {sam_out}")
