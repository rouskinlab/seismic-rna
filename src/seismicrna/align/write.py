from datetime import datetime
from itertools import chain
from logging import getLogger
from pathlib import Path
from shutil import copyfile
from typing import Iterable

from .fqops import FastqUnit, run_fastqc, run_cutadapt
from .report import AlignRefReport, AlignSampleReport
from .xamops import (run_bowtie2_build,
                     get_bowtie2_index_paths,
                     run_export,
                     run_xamgen)
from ..core import path
from ..core.arg import CMD_ALIGN, CMD_QC
from ..core.ngs import (count_single_paired,
                        count_total_reads,
                        run_flagstat,
                        run_ref_header,
                        run_index_xam,
                        run_index_fasta,
                        run_idxstats)
from ..core.parallel import dispatch
from ..core.seq import DNA, parse_fasta, write_fasta

logger = getLogger(__name__)


def write_temp_ref_files(temp_dir: Path,
                         refset_path: Path,
                         refs: set[str],
                         n_procs: int = 1):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, tuple[Path, Path]] = dict()
    if refs:
        # Parse the FASTA only if there are any references to write.
        for record in parse_fasta(refset_path, DNA):
            ref, _ = record
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.build(*path.FASTA_STEP_SEGS, top=temp_dir,
                                      step=path.STEP_ALIGN_INDEX_DEMULT,
                                      ref=ref, ext=refset_path.suffix)
                # Create the parent directory.
                ref_path.parent.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Created directory: {ref_path.parent}")
                try:
                    # Write the temporary FASTA file.
                    write_fasta(ref_path, [record])
                    # Build a Bowtie2 index of the temporary FASTA file.
                    index_prefix = ref_path.with_suffix("")
                    run_bowtie2_build(ref_path, index_prefix, n_procs=n_procs)
                except Exception as error:
                    logger.critical(
                        f"Failed to generate reference {ref_path}: {error}")
                else:
                    # Record the temporary FASTA and index prefix.
                    ref_paths[ref] = ref_path, index_prefix
    if missing := sorted(refs - set(ref_paths.keys())):
        logger.critical(f"Missing references in {refset_path}: "
                        + ", ".join(missing))
    return ref_paths


def fq_pipeline(fq_inp: FastqUnit,
                fasta: Path,
                bowtie2_index: Path, *,
                out_dir: Path,
                temp_dir: Path,
                keep_temp: bool,
                fastqc: bool,
                qc_extract: bool,
                cut: bool,
                cut_q1: int,
                cut_q2: int,
                cut_g1: str,
                cut_a1: str,
                cut_g2: str,
                cut_a2: str,
                cut_o: int,
                cut_e: float,
                cut_indels: bool,
                cut_nextseq: bool,
                cut_discard_trimmed: bool,
                cut_discard_untrimmed: bool,
                cut_m: int,
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
                cram: bool,
                n_procs: int = 1) -> list[Path]:
    """ Run all steps of the alignment pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    began = datetime.now()
    # Get attributes of the sample and references.
    sample = fq_inp.sample
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    # Determine the output directory for the finished XAM files.
    xams_out_dir = path.builddir(path.SampSeg,
                                 path.CmdSeg,
                                 top=out_dir,
                                 sample=sample,
                                 cmd=CMD_ALIGN)
    # Determine the path for FASTQC output files.
    if fq_inp.ref is not None:
        # If the input FASTQ files are demultiplexed, then include the
        # name of the reference in the path of the FASTQC output files.
        fqc_segs = path.FASTQC_DEMULT_SEGS
        fqc_vals = {path.TOP: out_dir,
                    path.SAMP: sample,
                    path.CMD: CMD_QC,
                    path.REF: fq_inp.ref}
        # Use the demultiplexed version of the AlignReport.
        report_type = AlignRefReport
    else:
        # Otherwise, use only the name of the sample, command, and step.
        fqc_segs = path.FASTQC_SEGS
        fqc_vals = {path.TOP: out_dir, path.SAMP: sample, path.CMD: CMD_QC}
        # Use the non-demultiplexed version of the AlignReport.
        report_type = AlignSampleReport
    if fastqc:
        # Run FASTQC on the input FASTQ files.
        fqc_out = path.build(*fqc_segs, **fqc_vals, step=path.STEP_QC_INIT)
        try:
            run_fastqc(fq_inp, fqc_out, extract=qc_extract, n_procs=n_procs)
        except Exception as error:
            logger.error(f"Failed to run FASTQC on {fq_inp}: {error}")
    if cut:
        # Trim adapters and low-quality bases with Cutadapt.
        fq_cut = fq_inp.to_new(path.StepSeg,
                               top=temp_dir,
                               sample=sample,
                               step=path.STEP_ALIGN_TRIM)
        run_cutadapt(fq_inp, fq_cut,
                     n_procs=n_procs,
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
                     cut_m=cut_m)
        if fastqc:
            # Run FASTQC after trimming with Cutadapt.
            fqc_out = path.build(*fqc_segs, **fqc_vals, step=path.STEP_QC_TRIM)
            try:
                run_fastqc(fq_cut, fqc_out, extract=qc_extract, n_procs=n_procs)
            except Exception as error:
                logger.error(f"Failed to run FASTQC on {fq_inp}: {error}")
    else:
        fq_cut = None
    # Align the FASTQ to the reference sequence using Bowtie2.
    xam_whole = path.build(*path.XAM_STEP_SEGS,
                           top=temp_dir,
                           sample=sample,
                           cmd=CMD_ALIGN,
                           step=path.STEP_ALIGN_MAP,
                           ref=refset,
                           ext=path.BAM_EXT)
    reads_align = run_xamgen(fq_inp if fq_cut is None else fq_cut,
                             xam_whole,
                             index_pfx=bowtie2_index,
                             n_procs=n_procs,
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
                                      else None))
    if not keep_temp and fq_cut is not None:
        # Delete the trimmed FASTQ file (do not delete the input FASTQ
        # file even if trimming was not performed).
        for fq_file in fq_cut.paths.values():
            fq_file.unlink(missing_ok=True)
    # The number of reads after trimming is defined as the number fed to
    # Bowtie 2, regardless of whether the reads were actually trimmed.
    reads_trim = reads_align.pop("reads", None)
    if reads_trim is None:
        raise RuntimeError("Failed to parse number of reads input to Bowtie2 "
                           f"(perhaps Bowtie2 failed): got {reads_align}")
    # If the reads were trimmed, then the initial number must be found
    # by counting the reads in the input FASTQ. Otherwise, the initial
    # number equals the number after trimming.
    reads_init = fq_inp.n_reads if cut else reads_trim
    # Index the whole XAM file to enable exporting only reads aligning
    # to each reference and to speed counting reads.
    run_index_xam(xam_whole, n_procs=n_procs)
    # Count the reads after filtering.
    flagstats = run_flagstat(xam_whole)
    paired_two, paired_one, singles = count_single_paired(flagstats)
    if fq_inp.paired:
        if singles:
            raise ValueError(f"{xam_whole} got {singles} single-end reads")
        reads_filter = {"paired-end, both mates mapped": paired_two,
                        "paired-end, one mate unmapped": paired_one}
    else:
        if n_paired := paired_two + paired_one:
            raise ValueError(f"{xam_whole} got {n_paired} paired-end reads")
        reads_filter = {"single-end": singles}
    # Guess how many reads mapped to each reference by the index stats.
    reads_refs = run_idxstats(xam_whole)
    # Guess which references received enough reads.
    sufficient_refs = {ref for ref, reads_ref in reads_refs.items()
                       if reads_ref >= min_reads}
    # Cache the sequence for each reference that received enough reads.
    ref_seqs = {ref: seq for ref, seq in parse_fasta(fasta, DNA)
                if ref in sufficient_refs}
    # Cache the header for each reference that received enough reads.
    ref_headers = {ref: header
                   for ref, header in run_ref_header(xam_whole, n_procs=n_procs)
                   if ref in sufficient_refs}
    # Determine the path of the output FASTA and its index.
    refs_file = xams_out_dir.joinpath(fasta.name)
    refs_file_index = refs_file.with_suffix(f"{refs_file.suffix}{path.FAI_EXT}")
    # In case the FASTA file and its index already exist, delete them.
    refs_file.unlink(missing_ok=True)
    refs_file_index.unlink(missing_ok=True)
    # Split the whole XAM file into one XAM file for each reference that
    # was guessed to have received enough reads.
    xams_out: list[Path] = list()
    # Iterate over list(sufficient_refs) instead of just sufficient_refs
    # so that items can be deleted from sufficient_refs without altering
    # the object being iterated over.
    for ref in list(sufficient_refs):
        try:
            # Export the reads that align to the given reference.
            xam_ref = path.build(*path.XAM_SEGS,
                                 top=out_dir,
                                 sample=sample,
                                 cmd=CMD_ALIGN,
                                 ref=ref,
                                 ext=(path.CRAM_EXT if cram else path.BAM_EXT))
            if xam_ref.parent != xams_out_dir:
                raise path.PathValueError(f"{xam_ref} is not in {xams_out_dir}")
            exp_kwargs = dict(ref=ref, header=ref_headers[ref], n_procs=n_procs)
            if cram:
                # Write the one reference sequence to a temporary FASTA.
                # Do NOT use overwrite=True because if refs_file is a
                # link, then the file it is linked to will be erased.
                write_fasta(refs_file, [(ref, ref_seqs[ref])])
                exp_kwargs.update(dict(ref_file=refs_file))
            run_export(xam_whole, xam_ref, **exp_kwargs)
        except Exception as error:
            logger.error(f"Failed to output {ref}: {error}")
        else:
            # Count the number of reads accurately and update the guess.
            reads_refs[ref] = count_total_reads(run_flagstat(xam_ref))
            logger.debug(f"File {xam_ref} received {reads_refs[ref]} reads")
            if reads_refs[ref] >= min_reads:
                xams_out.append(xam_ref)
            else:
                # Delete the XAM file if it did not get enough reads.
                xam_ref.unlink()
                # Delete the reference from sufficient_refs.
                sufficient_refs.remove(ref)
                logger.debug(f"Deleted {xam_ref} because it received "
                             f"{reads_refs[ref]} reads (< {min_reads})")
        finally:
            # Remove the temporary FASTA files.
            refs_file.unlink()
            refs_file_index.unlink()
    # The whole XAM file is no longer needed.
    if not keep_temp:
        xam_whole.unlink()
    if cram:
        try:
            # Make a hard link to the original FASTA file.
            refs_file.hardlink_to(fasta)
        except OSError as error:
            logger.warning(f"Copying {refs_file} to {fasta} because hard "
                           f"linking failed: {error}")
            copyfile(fasta, refs_file)
        # Index the new hard-linked FASTA file.
        run_index_fasta(refs_file)
    ended = datetime.now()
    if insufficient_refs := set(reads_refs) - sufficient_refs:
        logger.warning(f"Skipped references with fewer than {min_reads} reads: "
                       f"{sorted(insufficient_refs)}")
    # Write a report to summarize the alignment.
    report = report_type(sample=sample,
                         ref=fq_inp.ref,
                         paired_end=fq_inp.paired,
                         phred_enc=fq_inp.phred_enc,
                         fastqc=fastqc,
                         cut=cut,
                         cut_q1=cut_q1,
                         cut_q2=cut_q2,
                         cut_g1=list(cut_g1),
                         cut_a1=list(cut_a1),
                         cut_g2=list(cut_g2),
                         cut_a2=list(cut_a2),
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
                         reads_init=reads_init,
                         reads_trim=reads_trim,
                         reads_align=reads_align,
                         reads_filter=reads_filter,
                         reads_refs=reads_refs,
                         began=began,
                         ended=ended)
    report.save(out_dir, overwrite=True)
    # Return a list of name-sorted XAM files, each of which contains a
    # set of reads that all align to the same reference.
    return xams_out


def fqs_pipeline(fq_units: list[FastqUnit],
                 main_fasta: Path, *,
                 max_procs: int,
                 parallel: bool,
                 out_dir: Path,
                 temp_dir: Path,
                 keep_temp: bool,
                 **kwargs) -> list[Path]:
    """ Run all steps of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
    # Validate the maximum number of processes.
    if max_procs < 1:
        logger.warning("Maximum CPUs must be ≥ 1: setting to 1")
        max_procs = 1
    # Get the name of the reference for every demultiplexed FASTQ.
    temp_refs = set(filter(None, (fq_unit.ref for fq_unit in fq_units)))
    # Write a temporary FASTA file and Bowtie2 index for each
    # demultiplexed FASTQ.
    temp_fasta_paths = write_temp_ref_files(temp_dir, main_fasta,
                                            temp_refs, max_procs)
    # Check if the main FASTA file already has a Bowtie2 index.
    main_index = main_fasta.with_suffix("")
    if not all(index.is_file() for index in get_bowtie2_index_paths(main_index)):
        # Bowtie2 index does not already exist.
        main_index = None
    # Make the arguments for each alignment task.
    iter_args: list[tuple[FastqUnit, Path, Path]] = list()
    # One alignment task will be created for each FASTQ unit.
    for fq_unit in fq_units:
        if fq_unit.ref is not None:
            # If the FASTQ came from demultiplexing (so contains
            # reads from only one reference), then align to the
            # temporary FASTA file containing only that reference.
            try:
                temp_fasta, temp_index = temp_fasta_paths[fq_unit.ref]
            except KeyError:
                # If the FASTA with that reference does not exist,
                # then log an error and skip this FASTQ.
                logger.error(
                    f"Skipped {fq_unit} because reference '{fq_unit.ref}' "
                    f"was not found in FASTA file {main_fasta}")
                continue
            # Add these arguments to the lists of arguments that
            # will be passed to fq_pipeline.
            iter_args.append((fq_unit, temp_fasta, temp_index))
            logger.debug(f"Added task: align {fq_unit} to {temp_index}")
        else:
            # If the FASTQ may contain reads from ≥ 1 references,
            # then align to the FASTA file with all references.
            if main_index is None:
                # The FASTA of all the references does not already
                # have a Bowtie2 index, so build a temporary index.
                # Determine the name of the set of references.
                refset = path.parse(main_fasta, path.FastaSeg)[path.REF]
                # Determine the path of the temporary Bowtie2 index
                # of the main FASTA file.
                main_index = path.build(*path.FASTA_INDEX_DIR_STEP_SEGS,
                                        top=temp_dir,
                                        step=path.STEP_ALIGN_INDEX,
                                        ref=refset)
                # Make its parent directory if it does not exist.
                main_index.parent.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Created directory: {main_index.parent}")
                # Build the Bowtie2 index.
                try:
                    run_bowtie2_build(main_fasta, main_index,
                                      n_procs=max_procs)
                    # Create a symbolic link to the reference file in
                    # the same directory as the new index.
                    fasta_link = main_index.with_suffix(main_fasta.suffix)
                    fasta_link.symlink_to(main_fasta)
                    # Add the FASTA link and the Bowtie2 index to the
                    # set of files to delete after alignment finishes.
                    # Being deleted is the only purpose of fasta_link.
                    temp_fasta_paths[refset] = fasta_link, main_index
                except Exception as error:
                    logger.critical(
                        f"Failed to index {main_fasta} with Bowtie2: {error}")
                    # Reset main_index to None and skip this FASTQ unit.
                    main_index = None
                    continue
            # Add these arguments to the lists of arguments that
            # will be passed to fq_pipeline. Note that main_index
            # could be a pre-built index in the same directory as
            # main_fasta or a temporary index that is deleted when
            # alignment finishes; but only in the latter case is it
            # added to temp_fasta_paths.
            iter_args.append((fq_unit, main_fasta, main_index))
            logger.debug(f"Added task: align {fq_unit} to {main_index}")
    # Generate alignment map (XAM) files.
    xam_lists = dispatch(fq_pipeline,
                         max_procs,
                         parallel,
                         hybrid=True,
                         args=iter_args,
                         kwargs=dict(out_dir=out_dir,
                                     temp_dir=temp_dir,
                                     keep_temp=keep_temp,
                                     **kwargs))
    xams = list(chain(*xam_lists))
    if not keep_temp:
        # Delete the temporary files.
        for ref_file, index_prefix in temp_fasta_paths.values():
            # Reference file
            ref_file.unlink(missing_ok=True)
            logger.debug(f"Deleted temporary reference file: {ref_file}")
            # Index files
            for index_file in get_bowtie2_index_paths(index_prefix):
                index_file.unlink(missing_ok=True)
                logger.debug(f"Deleted temporary index file: {index_file}")
    # Return the final alignment map (XAM) files.
    return xams


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
                logger.error(f"No reference '{fq_unit.ref}' for {fq_unit}")
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
        xams_new = set(fqs_pipeline(fqs_to_align,
                                    fasta,
                                    out_dir=out_dir,
                                    **kwargs))
    else:
        logger.warning("All given FASTQ files have already been aligned")
        xams_new = set()
    # Merge the existing and new XAM paths into a tuple of strings.
    return list(xams_extant | xams_new)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
