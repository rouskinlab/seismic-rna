from collections import defaultdict
from datetime import datetime
from pathlib import Path
from shutil import rmtree
from typing import Iterable

from .fqunit import FastqUnit
from .report import AlignRefReport, AlignSampleReport
from .xamops import (FASTP_PHRED_OUT,
                     run_bowtie2_build,
                     get_bowtie2_index_paths,
                     run_export,
                     run_flags,
                     run_realign,
                     run_xamgen)
from ..core import path
from ..core.arg import CMD_ALIGN
from ..core.logs import logger
from ..core.ngs import (FLAG_PAIRED,
                        FLAG_PROPER,
                        FLAG_FIRST,
                        FLAG_SECOND,
                        FLAG_REVERSE,
                        count_single_paired,
                        count_total_reads,
                        run_flagstat,
                        run_ref_header,
                        run_index_xam,
                        run_idxstats,
                        xam_paired)
from ..core.seq import DNA, get_fasta_seq, parse_fasta, write_fasta
from ..core.task import dispatch


def format_ref_reverse(ref: str, rev_label: str):
    """ Name the reverse strand of a reference. """
    if not rev_label:
        raise ValueError("rev_label must have a value, "
                         f"but got {repr(rev_label)}")
    return f"{ref}{rev_label}"


def write_tmp_ref_files(tmp_dir: Path,
                        refset_path: Path,
                        refs: set[str],
                        n_procs: int = 1):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, tuple[Path, Path]] = dict()
    if refs:
        logger.routine("Began writing temporary FASTA files of references "
                       f"{refs} to {tmp_dir}")
        # Parse the FASTA only if there are any references to write.
        for record in parse_fasta(refset_path, DNA):
            ref, _ = record
            if ref in refs:
                logger.detail(f"Writing FASTA file of reference {repr(ref)}")
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.build(*path.FASTA_STAGE_SEGS,
                                      top=tmp_dir,
                                      stage=path.STAGE_ALIGN_INDEX_DEMULT,
                                      ref=ref,
                                      ext=refset_path.suffix)
                # Create the parent directory.
                ref_path.parent.mkdir(parents=True, exist_ok=True)
                try:
                    # Write the temporary FASTA file.
                    write_fasta(ref_path, [record])
                    # Build a Bowtie2 index of the temporary FASTA file.
                    index_prefix = ref_path.with_suffix("")
                    run_bowtie2_build(ref_path, index_prefix, n_procs=n_procs)
                except Exception as error:
                    # If anything goes wrong when writing and indexing
                    # this reference, log an error message but continue
                    # with the other references.
                    logger.error(error)
                else:
                    # Record the temporary FASTA and index prefix.
                    ref_paths[ref] = ref_path, index_prefix
            else:
                logger.detail(f"Skipped unused reference {repr(ref)}")
    missing = sorted(refs - set(ref_paths.keys()))
    if missing:
        # If any references in refs do not have sequences, then log an
        # error but continue with the references that have sequences.
        logger.error("".join([f"Missing sequences in {refset_path}: ",
                              ", ".join(missing)]))
    return ref_paths


def calc_flags_sep_strands(f1r2_fwd: bool, paired: bool, bt2_mixed: bool):
    """ Calculate flags for separating strands. """
    if paired:
        flags_req_fwd = [FLAG_FIRST | FLAG_PAIRED, FLAG_SECOND | FLAG_PAIRED]
        flags_exc_fwd = [FLAG_SECOND, FLAG_FIRST]
        flags_req_rev = [FLAG_FIRST | FLAG_PAIRED, FLAG_SECOND | FLAG_PAIRED]
        flags_exc_rev = [FLAG_SECOND, FLAG_FIRST]
        if bt2_mixed:
            # When using mixed mode, some paired-end reads may not have
            # a mate; thus, the reads in the BAM file are not guaranteed
            # to alternate between read 1 and read 2.
            # When realigning paired-end reads to the reverse strand,
            # the interleaved BAM file is converted to an interleaved
            # FASTQ file; the presence of any paired-end read without a
            # mate (singleton) will thus throw the interleaving out of
            # register, pairing up reads that are not actually mates.
            # To circumvent this problem, only properly paired reads are
            # allowed to be realigned.
            # It is possible to realign singleton paired-end reads, but
            # only by considering them to be single-ended, aligning them
            # separately from the properly paired reads, and merging the
            # results of the two alignments into one BAM file.
            # That is doable, but much more complicated and only worth
            # implementing if there is a clear need for this feature.
            logger.warning("When mixed mode (--bt2-mixed) is combined with "
                           "strand separation (--sep-strands), any paired-end "
                           "read that failed to align properly with its mate "
                           "will be dropped, not aligned to the reverse strand")
            flags_req_fwd[0] |= FLAG_PROPER
            flags_req_fwd[1] |= FLAG_PROPER
            flags_req_rev[0] |= FLAG_PROPER
            flags_req_rev[1] |= FLAG_PROPER
    else:
        flags_req_fwd = [0]
        flags_exc_fwd = [FLAG_PAIRED]
        flags_req_rev = [0]
        flags_exc_rev = [FLAG_PAIRED]
    if f1r2_fwd:
        # Align forward read 1s and reverse read 2s to the + strand.
        # Align reverse read 1s and forward read 2s to the - strand.
        if paired:
            flags_req_fwd[1] |= FLAG_REVERSE
            flags_exc_rev[1] |= FLAG_REVERSE
        flags_exc_fwd[0] |= FLAG_REVERSE
        flags_req_rev[0] |= FLAG_REVERSE
    else:
        # Align reverse read 1s and forward read 2s to the + strand.
        # Align forward read 1s and reverse read 2s to the - strand.
        if paired:
            flags_req_rev[1] |= FLAG_REVERSE
            flags_exc_fwd[1] |= FLAG_REVERSE
        flags_exc_rev[0] |= FLAG_REVERSE
        flags_req_fwd[0] |= FLAG_REVERSE
    flags = (flags_req_fwd, flags_exc_fwd), (flags_req_rev, flags_exc_rev)
    logger.detail(
        f"Calculated SAM flags for separating strands with parameters "
        f"f1r2_fwd={f1r2_fwd}, paired={paired}, bt2_mixed={bt2_mixed}: {flags}"
    )
    return flags


def separate_strands(xam_file: Path,
                     fasta: Path, *,
                     paired: bool | None,
                     bt2_mixed: bool,
                     rev_label: str,
                     f1r2_fwd: bool,
                     min_mapq: int,
                     keep_tmp: bool,
                     n_procs: int = 1,
                     **kwargs):
    """ Separate a XAM file into two XAM files of reads that aligned to
    the forward and reverse strands, respectively. """
    logger.routine(
        f"Began separating forward and reverse strands in {xam_file}"
    )
    if paired is None:
        paired = xam_paired(run_flagstat(xam_file, n_procs=n_procs))
    out_dir = xam_file.parent
    ref = path.parse(xam_file, *path.XAM_SEGS)[path.REF]
    ref_rev = format_ref_reverse(ref, rev_label)
    # Calculate the flags.
    ((flags_req_fwd, flags_exc_fwd),
     (flags_req_rev, flags_exc_rev)) = calc_flags_sep_strands(f1r2_fwd,
                                                              paired,
                                                              bt2_mixed)
    # Make a temporary directory for all splitting strand operations.
    tmp_dir = out_dir.joinpath(ref)
    tmp_dir.mkdir(parents=False, exist_ok=False)
    logger.detail(f"Created temporary directory {tmp_dir} for aligning "
                  f"reverse-strand reads in {xam_file} to {repr(ref_rev)}")
    try:
        # Write the reverse-strand reference sequence to a FASTA file.
        index_dir = tmp_dir.joinpath("index")
        index_dir.mkdir()
        fasta_rev = index_dir.joinpath(ref_rev).with_suffix(path.FASTA_EXTS[0])
        refseq = get_fasta_seq(fasta, DNA, ref)
        write_fasta(fasta_rev, [(ref_rev, refseq.rc)])
        # Index the reverse-strand reference.
        index_rev = fasta_rev.with_suffix("")
        run_bowtie2_build(fasta_rev, index_rev, n_procs=n_procs)
        # Re-align the reads that had aligned to the reverse strand of
        # the forward-strand reference to the reverse-strand reference.
        bam_rev = out_dir.joinpath(ref_rev).with_suffix(path.BAM_EXT)
        realign_dir = tmp_dir.joinpath("realign")
        realign_dir.mkdir()
        run_realign(xam_file,
                    bam_rev,
                    tmp_pfx=realign_dir.joinpath(ref_rev),
                    index_pfx=index_rev,
                    paired=paired,
                    flags_req=flags_req_rev,
                    flags_exc=flags_exc_rev,
                    min_mapq=min_mapq,
                    bt2_mixed=bt2_mixed,
                    n_procs=n_procs,
                    **kwargs)
        if not keep_tmp:
            rmtree(index_dir)
        # Extract the reads that had aligned to the forward strand.
        bam_fwd = realign_dir.joinpath(ref).with_suffix(path.BAM_EXT)
        run_flags(xam_file,
                  bam_fwd,
                  tmp_pfx=realign_dir.joinpath(ref),
                  flags_req=flags_req_fwd,
                  flags_exc=flags_exc_fwd,
                  n_procs=n_procs)
        # Renaming overwrites the original BAM file of both strands.
        bam_fwd.rename(xam_file)
        logger.detail(f"Overwrote {xam_file} with only the forward-stand reads "
                      f"from {bam_fwd}")
        logger.routine(
            f"Ended separating forward and reverse strands in {xam_file}"
        )
        return bam_rev
    finally:
        # Make sure to delete tmp_dir if keep_tmp is False because if it
        # exists and this function is run again, then tmp_dir.mkdir()
        # will fail.
        if not keep_tmp:
            rmtree(tmp_dir)
            logger.detail(f"Deleted temporary directory {tmp_dir}")


def extract_reference(ref: str,
                      header: str, *,
                      xam_whole: Path,
                      fasta: Path,
                      sample: str,
                      top: Path,
                      min_reads: int,
                      sep_strands: bool,
                      n_procs: int = 1,
                      **kwargs):
    """ Extract one reference from a XAM file. """
    logger.routine(f"Began extracting reference {repr(ref)} from {xam_whole}")
    if min_reads < 0:
        min_reads = 0
        logger.warning(f"min_reads must be ≥ 0, but got {min_reads}: set to 0")
    # Export the reads that align to the given reference.
    xam_ref = path.build(*path.XAM_SEGS,
                         top=top,
                         sample=sample,
                         cmd=CMD_ALIGN,
                         ref=ref,
                         ext=path.BAM_EXT)
    run_export(xam_whole,
               xam_ref,
               ref=ref,
               header=header,
               n_procs=n_procs)
    xam_files = [xam_ref]
    if sep_strands:
        try:
            # Split the XAM file into forward and reverse strands.
            xam_rev = separate_strands(xam_ref,
                                       fasta,
                                       n_procs=n_procs,
                                       **kwargs)
            xam_files.append(xam_rev)
        except Exception:
            # Delete the XAM file containing both strands because its
            # name is the same as the file of only forward-strand reads.
            # If not deleted, it would remain in the output directory
            # and could be given to a future step that expects only
            # forward-strand reads, causing unintended behavior.
            xam_ref.unlink()
            logger.warning(f"Deleted {xam_ref} because separating it into "
                           "forward and reverse strands failed")
            raise
    # Count the reads in each XAM file; delete files with too few.
    nums_reads = dict()
    for xam in xam_files:
        try:
            ref = path.parse(xam, *path.XAM_SEGS)[path.REF]
            if ref in nums_reads:
                raise ValueError(f"Duplicate reference: {repr(ref)}")
            num_reads = count_total_reads(run_flagstat(xam, n_procs=n_procs))
            logger.detail(f"{xam} has {num_reads} read(s)")
            if num_reads < min_reads:
                logger.warning(
                    f"Skipped sample {repr(sample)} reference {repr(ref)}: "
                    f"{num_reads} < {min_reads} read(s)"
                )
                xam.unlink()
                logger.detail(f"Deleted {xam}")
        except Exception as error:
            logger.error(error)
            xam.unlink()
            logger.detail(f"Deleted {xam}")
        else:
            nums_reads[ref] = num_reads
    logger.routine(f"Ended extracting reference {repr(ref)} from {xam_whole}")
    return nums_reads


def split_references(xam_whole: Path, *,
                     fasta: Path,
                     paired: bool | None,
                     phred_arg: str,
                     top: Path,
                     keep_tmp: bool,
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
                     min_mapq: int,
                     min_reads: int,
                     sep_strands: bool,
                     f1r2_fwd: bool,
                     rev_label: str,
                     n_procs: int = 1):
    """ Split a XAM file into one file per reference. """
    logger.routine(f"Began splitting {xam_whole} by reference")
    sample = path.parse(xam_whole, *path.XAM_SEGS)[path.SAMP]
    # Guess how many reads mapped to each reference by the index stats.
    refs_counts = run_idxstats(xam_whole)
    # Guess which references received enough reads.
    guess_refs = {ref for ref, count in refs_counts.items()
                  if count >= min_reads}
    logger.detail(
        f"Guessed that there are ≥ {min_reads} read(s) in each of the "
        f"{len(guess_refs)} reference(s) {sorted(guess_refs)}"
    )
    # Cache the header for each reference that received enough reads.
    ref_headers = {ref: header
                   for ref, header in run_ref_header(xam_whole, n_procs=n_procs)
                   if ref in guess_refs}
    logger.detail(f"Cached SAM headers for the {len(guess_refs)} reference(s) "
                  f"that were guessed to have ≥ {min_reads} read(s)")
    # Split the whole XAM file into one XAM file for each reference that
    # was guessed to have received enough reads.
    nums_reads = dispatch(extract_reference,
                          max_procs=n_procs,
                          pass_n_procs=True,
                          args=list(ref_headers.items()),
                          kwargs=dict(xam_whole=xam_whole,
                                      fasta=fasta,
                                      sample=sample,
                                      paired=paired,
                                      phred_arg=phred_arg,
                                      top=top,
                                      keep_tmp=keep_tmp,
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
                                      sep_strands=sep_strands,
                                      f1r2_fwd=f1r2_fwd,
                                      rev_label=rev_label,
                                      min_reads=min_reads))
    # Collect the number of reads for each reference into one dict.
    refs_counts = dict()
    for num_reads in nums_reads:
        for ref, count in num_reads.items():
            logger.detail(f"Reference {repr(ref)} got {count} read(s)")
            if ref in refs_counts:
                logger.error(f"Reference {repr(ref)} is duplicated")
                xam_ref = path.build(*path.XAM_SEGS,
                                     top=top,
                                     sample=sample,
                                     cmd=CMD_ALIGN,
                                     ref=ref,
                                     ext=path.BAM_EXT)
                try:
                    xam_ref.unlink()
                    logger.detail(f"Deleted {xam_ref}")
                except OSError:
                    pass
            else:
                refs_counts[ref] = count
    # Sort the references in decreasing order of the number of aligned
    # reads, and then alphabetically in case of a tie.
    counts_refs = defaultdict(list)
    for ref, count in refs_counts.items():
        counts_refs[count].append(ref)
    refs_counts = dict()
    for count in sorted(counts_refs, reverse=True):
        for ref in sorted(counts_refs[count]):
            refs_counts[ref] = count
    logger.routine(f"Ended splitting {xam_whole} by reference")
    return refs_counts


def fq_pipeline(fq_inp: FastqUnit,
                fasta: Path,
                bowtie2_index: Path, *,
                out_dir: Path,
                tmp_dir: Path,
                keep_tmp: bool,
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
                fastp_adapter_fasta: Path | None,
                fastp_detect_adapter_for_pe: bool,
                fastp_min_length: int,
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
                min_mapq: int,
                min_reads: int,
                sep_strands: bool,
                f1r2_fwd: bool,
                rev_label: str,
                n_procs: int = 1) -> list[Path]:
    """ Run all stages of the alignment pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    began = datetime.now()
    logger.routine(f"Began processing {fq_inp} through the alignment pipeline")
    # Get attributes of the sample and references.
    sample = fq_inp.sample
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    # Create the output directory: this is necessary for FASTQ files of
    # unaligned reads to have a place to be written.
    align_dir = path.builddir(*path.CMD_DIR_SEGS,
                              top=out_dir,
                              sample=sample,
                              cmd=path.CMD_ALIGN_DIR)
    # Optionally trim the reads with Fastp, and then align them to the
    # reference sequence with Bowtie2.
    xam_whole = path.buildpar(*path.XAM_STAGE_SEGS,
                              top=tmp_dir,
                              sample=sample,
                              cmd=path.CMD_ALIGN_DIR,
                              stage=path.STAGE_ALIGN_MAP,
                              ref=refset,
                              ext=path.BAM_EXT)
    reads_align = run_xamgen(
        fq_inp,
        xam_whole,
        fastp=fastp,
        fastp_dir=align_dir,
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
        index_pfx=bowtie2_index,
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
        fq_unal=(path.build(path.SampSeg,
                            path.CmdSeg,
                            path.DmFastqSeg,
                            top=out_dir,
                            sample=sample,
                            cmd=CMD_ALIGN,
                            ref=(f"{fq_inp.ref}__unaligned"
                                 if fq_inp.ref is not None
                                 else "unaligned"),
                            ext=path.FQ_EXTS[0])
                 if bt2_un
                 else None),
        n_procs=n_procs
    )
    # The number of reads after trimming is defined as the number fed to
    # Bowtie 2, regardless of whether the reads were actually trimmed.
    reads_trim = reads_align.pop("reads", None)
    if reads_trim is None:
        raise RuntimeError("Failed to parse number of reads input to Bowtie2 "
                           f"(perhaps Bowtie2 failed): got {reads_align}")
    logger.detail(f"Determined Bowtie 2 received {reads_trim} reads")
    if fastp:
        # If the reads were trimmed, then the initial number must be
        # found by counting the reads in the input FASTQ.
        reads_init = fq_inp.n_reads
    else:
        # Otherwise, the initial number of reads equals the number fed
        # to Bowtie 2, so we can save time by using that number.
        reads_init = reads_trim
    logger.detail(f"Determined {fq_inp} contained {reads_init} reads")
    # Index the whole XAM file to enable exporting only reads aligning
    # to each reference and to speed counting reads.
    run_index_xam(xam_whole, n_procs=n_procs)
    # Count the reads after filtering.
    flagstats = run_flagstat(xam_whole)
    paired_two, paired_one, singles = count_single_paired(flagstats)
    if fq_inp.paired:
        if singles:
            raise RuntimeError(f"{xam_whole} has {singles} single-end reads")
        reads_filter = {"paired-end, both mates mapped": paired_two,
                        "paired-end, one mate unmapped": paired_one}
    else:
        if n_paired := paired_two + paired_one:
            raise RuntimeError(f"{xam_whole} has {n_paired} paired-end reads")
        reads_filter = {"single-end": singles}
    logger.detail(
        f"Determined {xam_whole} contained {reads_filter} reads after filtering"
    )
    # Split the whole XAM file into one XAM file for each reference.
    reads_refs = split_references(xam_whole,
                                  fasta=fasta,
                                  paired=fq_inp.paired,
                                  phred_arg=(FASTP_PHRED_OUT if fastp
                                             else fq_inp.phred_arg),
                                  top=out_dir,
                                  keep_tmp=keep_tmp,
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
                                  sep_strands=sep_strands,
                                  f1r2_fwd=f1r2_fwd,
                                  rev_label=rev_label,
                                  min_reads=min_reads,
                                  n_procs=n_procs)
    if not keep_tmp:
        # Delete the BAM file of all references.
        xam_whole.unlink(missing_ok=True)
        logger.detail(f"Deleted {xam_whole}")
    ended = datetime.now()
    # Write a report to summarize the alignment.
    if fq_inp.ref is not None:
        # Use the demultiplexed version of the AlignReport.
        report_type = AlignRefReport
    else:
        # Use the non-demultiplexed version of the AlignReport.
        report_type = AlignSampleReport
    report = report_type(sample=sample,
                         ref=fq_inp.ref,
                         paired_end=fq_inp.paired,
                         phred_enc=fq_inp.phred_enc,
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
                         fastp_adapter_fasta=fastp_adapter_fasta,
                         fastp_detect_adapter_for_pe=fastp_detect_adapter_for_pe,
                         fastp_min_length=fastp_min_length,
                         bt2_local=bt2_local,
                         bt2_discordant=bt2_discordant,
                         bt2_mixed=bt2_mixed,
                         bt2_dovetail=bt2_dovetail,
                         bt2_contain=bt2_contain,
                         bt2_score_min=(bt2_score_min_loc if bt2_local
                                        else bt2_score_min_e2e),
                         bt2_i=bt2_i,
                         bt2_x=bt2_x,
                         bt2_gbar=bt2_gbar,
                         bt2_l=bt2_l,
                         bt2_s=bt2_s,
                         bt2_d=bt2_d,
                         bt2_r=bt2_r,
                         bt2_dpad=bt2_dpad,
                         bt2_orient=bt2_orient,
                         bt2_un=bt2_un,
                         min_mapq=min_mapq,
                         sep_strands=sep_strands,
                         f1r2_fwd=f1r2_fwd,
                         rev_label=rev_label,
                         min_reads=min_reads,
                         align_reads_init=reads_init,
                         reads_trim=reads_trim,
                         reads_align=reads_align,
                         reads_filter=reads_filter,
                         reads_refs=reads_refs,
                         began=began,
                         ended=ended)
    report_saved = report.save(out_dir, force=True)
    logger.routine(f"Ended processing {fq_inp} through the alignment pipeline")
    return report_saved.parent


def fqs_pipeline(fq_units: list[FastqUnit],
                 main_fasta: Path, *,
                 max_procs: int,
                 out_dir: Path,
                 tmp_dir: Path,
                 keep_tmp: bool,
                 **kwargs) -> list[Path]:
    """ Run all stages of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
    logger.routine(f"Began running the alignment pipeline")
    # Validate the maximum number of processes.
    if max_procs < 1:
        logger.warning("max_procs must be ≥ 1: setting to 1")
        max_procs = 1
    # Get the name of the reference for every demultiplexed FASTQ.
    tmp_refs = set(filter(None, (fq_unit.ref for fq_unit in fq_units)))
    if tmp_refs:
        logger.detail(f"Found {len(tmp_refs)} references among demultiplexed "
                      f"FASTQ files: {sorted(tmp_refs)}")
    # Write a temporary FASTA file and Bowtie2 index for each
    # demultiplexed FASTQ.
    tmp_fasta_paths = write_tmp_ref_files(tmp_dir,
                                          main_fasta,
                                          tmp_refs,
                                          max_procs)
    # Check if the main FASTA file already has a Bowtie2 index.
    main_index = main_fasta.with_suffix("")
    if all(index.is_file() for index in get_bowtie2_index_paths(main_index)):
        logger.detail(f"A Bowtie 2 index exists for {main_fasta}")
    else:
        logger.detail(f"A Bowtie 2 index does not exist for {main_fasta}")
        main_index = None
    # Make the arguments for each alignment task.
    iter_args: list[tuple[FastqUnit, Path, Path]] = list()
    # One alignment task will be created for each FASTQ unit.
    for fq_unit in fq_units:
        logger.detail(f"Preparing to align {fq_unit}")
        if fq_unit.ref is not None:
            logger.detail(f"{fq_unit} contains reads from 1 reference, "
                          f"{repr(fq_unit.ref)}")
            # If the FASTQ came from demultiplexing (so contains
            # reads from only one reference), then align to the
            # temporary FASTA file containing only that reference.
            try:
                tmp_fasta, tmp_index = tmp_fasta_paths[fq_unit.ref]
            except KeyError:
                # If the FASTA with that reference does not exist,
                # then log an error and skip this FASTQ.
                logger.error(f"Skipped {fq_unit} because {main_fasta} "
                             f"does not contain reference {repr(fq_unit.ref)}")
                continue
            # Add these arguments to the lists of arguments that will be
            # passed to fq_pipeline.
            iter_args.append((fq_unit, tmp_fasta, tmp_index))
            logger.detail(f"Planning to align {fq_unit} to {tmp_fasta}")
        else:
            logger.detail(f"{fq_unit} may contain reads from ≥ 1 reference")
            # If the FASTQ may contain reads from ≥ 1 references,
            # then align to the FASTA file with all references.
            if main_index is None:
                # The FASTA of all the references does not already
                # have a Bowtie2 index, so build a temporary index.
                # Determine the name of the set of references.
                refset = path.parse(main_fasta, path.FastaSeg)[path.REF]
                # Determine the path of the temporary Bowtie 2 index
                # of the main FASTA file.
                main_index = path.build(*path.FASTA_INDEX_DIR_STAGE_SEGS,
                                        top=tmp_dir,
                                        stage=path.STAGE_ALIGN_INDEX,
                                        ref=refset)
                # Make its parent directory if it does not exist.
                main_index.parent.mkdir(parents=True, exist_ok=True)
                logger.detail(f"Created directory {main_index.parent} "
                              f"for Bowtie 2 index of {main_fasta}")
                # Build the Bowtie2 index.
                try:
                    run_bowtie2_build(main_fasta,
                                      main_index,
                                      n_procs=max_procs)
                    # Create a symbolic link to the reference file in
                    # the same directory as the new index.
                    fasta_link = main_index.with_suffix(main_fasta.suffix)
                    fasta_link.symlink_to(main_fasta)
                    logger.detail("Created a temporary symbolic link "
                                  f"{fasta_link} pointing to {main_fasta}")
                    # Add the FASTA link and the Bowtie 2 index to the
                    # set of files to delete after alignment finishes.
                    # Being deleted is the only purpose of fasta_link.
                    tmp_fasta_paths[refset] = fasta_link, main_index
                except Exception as error:
                    logger.error(error)
                    # Reset main_index to None and skip this FASTQ unit.
                    main_index = None
                    continue
            # Add these arguments to the lists of arguments that
            # will be passed to fq_pipeline. Note that main_index
            # could be a pre-built index in the same directory as
            # main_fasta or a temporary index that is deleted when
            # alignment finishes; but only in the latter case is it
            # added to tmp_fasta_paths.
            iter_args.append((fq_unit, main_fasta, main_index))
            logger.detail(f"Planning to align {fq_unit} to {main_fasta}")
    # Generate alignment map (XAM) files.
    xam_dirs = dispatch(fq_pipeline,
                        max_procs,
                        args=iter_args,
                        kwargs=dict(out_dir=out_dir,
                                    tmp_dir=tmp_dir,
                                    keep_tmp=keep_tmp,
                                    **kwargs))
    # Return the final alignment map (XAM) directories.
    logger.routine(f"Ended running the alignment pipeline")
    return xam_dirs


def figure_alignments(fq_units: list[FastqUnit], refs: set[str]):
    """ Every expected alignment of a sample to a reference. """
    # Map each combination of a sample and reference to a FASTQ unit.
    alignments: dict[tuple[str, str], FastqUnit] = dict()
    # Keep track of any duplicate sample-reference pairs.
    duplicates: set[tuple[str, str]] = set()
    for fq_unit in fq_units:
        # Determine which references the FASTQ reads could come from.
        if fq_unit.ref is None:
            # The FASTQ contains reads from potentially all references.
            fq_refs = refs
        else:
            # The FASTQ contains reads from only one reference.
            # Confirm that the reference actually exists.
            if fq_unit.ref not in refs:
                logger.error(f"No reference {repr(fq_unit.ref)} for {fq_unit}")
                continue
            fq_refs = {fq_unit.ref}
        # Add each sample-reference pair to the expected alignments.
        for ref in fq_refs:
            sample_ref = fq_unit.sample, ref
            if sample_ref in duplicates:
                # Skip the sample-reference pair if it is a duplicate.
                continue
            try:
                # Test if the sample-reference pair is already in the
                # dict of alignments. If so, then remove it.
                alignments.pop(sample_ref)
            except KeyError:
                # If not, then add the FASTQ to the dict of alignments,
                # keyed by its sample-reference pair.
                alignments[sample_ref] = fq_unit
            else:
                # If so, then flag it as a duplicate.
                logger.warning(f"Duplicate sample and reference: {sample_ref}")
                duplicates.add(sample_ref)
    # Return a duplicate-free dict of alignments.
    return alignments


def check_fqs_xams(alignments: dict[tuple[str, str], FastqUnit],
                   out_dir: Path):
    """ Return every FASTQ unit on which alignment must be run and every
    expected XAM file that already exists. """
    fqs_missing: dict[tuple[str, str], FastqUnit] = dict()
    xams_extant: list[Path] = list()
    for (sample, ref), fq_unit in alignments.items():
        # Determine the path of the XAM file expected to result from the
        # alignment of the sample to the reference.
        for ext in path.XAM_EXTS:
            xam_expect = path.build(*path.XAM_SEGS,
                                    top=out_dir,
                                    cmd=CMD_ALIGN,
                                    sample=sample,
                                    ref=ref,
                                    ext=ext)
            if xam_expect.is_file():
                # If the XAM file already exists, then add it to the
                # dict of XAM files that have already been aligned.
                xams_extant.append(xam_expect)
                break
        else:
            # If at least one XAM file for a FASTQ unit does not exist,
            # then align the FASTQ.
            fqs_missing[sample, ref] = fq_unit
    return fqs_missing, xams_extant


def merge_nondemult_fqs(fq_units: Iterable[FastqUnit]):
    """ For every FASTQ that is not demultiplexed, merge all the keys
    that map to the FASTQ into one key: (sample, None). Merging ensures
    that every non-demultiplexed FASTQ is aligned only once to the whole
    set of references, not once for every reference in the set. This
    function is essentially the inverse of `figure_alignments`. """
    merged: dict[tuple[str, str | None], FastqUnit] = dict()
    for fq_unit in fq_units:
        merged[fq_unit.sample, fq_unit.ref] = fq_unit
    return list(merged.values())


def list_fqs_xams(fq_units: list[FastqUnit],
                  refs: set[str],
                  out_dir: Path):
    """ List every FASTQ to align and every extant XAM file. """
    # Determine all possible alignments of a sample and reference.
    alignments = figure_alignments(fq_units, refs)
    # Determine which alignments need to be or have already been run.
    fqs_missing, xams_extant = check_fqs_xams(alignments, out_dir)
    # Merge entries for each non-demultiplexed FASTQ.
    fqs_merged = merge_nondemult_fqs(fqs_missing.values())
    return fqs_merged, set(xams_extant)


def align_samples(fq_units: list[FastqUnit],
                  fasta: Path, *,
                  out_dir: Path,
                  force: bool,
                  **kwargs) -> list[Path]:
    """ Run the alignment pipeline and return a tuple of all XAM files
    from the pipeline. """
    if not fq_units:
        logger.warning("No FASTQ files or pairs of FASTQ files were given")
        return list()
    if force:
        # force all alignments.
        fqs_to_align = fq_units
        xams_extant = set()
    else:
        # Get the names of all reference sequences.
        refs = set(parse_fasta(fasta, None))
        # Run only the alignments whose outputs do not yet exist.
        fqs_to_align, xams_extant = list_fqs_xams(fq_units, refs, out_dir)
    if fqs_to_align:
        # Align all FASTQs that need to be aligned.
        xam_dirs_new = set(fqs_pipeline(fqs_to_align,
                                        fasta,
                                        out_dir=out_dir,
                                        **kwargs))
    else:
        logger.warning("All given FASTQ files have already been aligned")
        xam_dirs_new = set()
    # Merge the existing and new XAM paths into a tuple of strings.
    return list({xam.parent for xam in xams_extant} | xam_dirs_new)

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
