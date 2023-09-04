from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Iterable

from .fqops import FastqUnit, run_fastqc, run_cutadapt
from .report import AlignReport
from .xamgen import (build_bowtie2_index, get_bowtie2_index_paths,
                     run_export, run_xamgen)
from ..core import path
from ..core.cmd import CMD_ALIGN, CMD_QC
from ..core.parallel import dispatch
from ..core.seq import parse_fasta, write_fasta
from ..core.xam import (count_single_paired, count_total_records,
                        get_bam_index, index_bam, run_flagstat)

logger = getLogger(__name__)

OUT_EXT = path.CRAM_EXT


def write_temp_ref_files(temp_dir: Path,
                         refset_path: Path,
                         refs: set[str],
                         n_procs: int = 1):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, tuple[Path, Path]] = dict()
    if refs:
        # Parse the FASTA only if there are any references to write.
        for record in parse_fasta(refset_path):
            ref = record[0]
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.build(*path.FASTA_STEP_SEGS,
                                      top=temp_dir, step=path.STEPS_ALIGN[0],
                                      ref=ref, ext=refset_path.suffix)
                # Create the parent directory.
                ref_path.parent.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Created directory: {ref_path.parent}")
                try:
                    # Write the temporary FASTA file.
                    write_fasta(ref_path, [record])
                    # Build a Bowtie2 index of the temporary FASTA file.
                    index_prefix = ref_path.with_suffix("")
                    build_bowtie2_index(ref_path, index_prefix, n_procs=n_procs)
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
                save_temp: bool,
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
                bt2_orient: str,
                min_mapq: int,
                min_reads: int,
                n_procs: int = 1) -> list[Path]:
    """ Run all steps of the alignment pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    # Get attributes of the sample and references.
    sample = fq_inp.sample
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    refs = list(ref for ref, _ in parse_fasta(fasta))
    # Determine the path for FASTQC output files.
    if fq_inp.ref:
        # If the input FASTQ files are demultiplexed, then include the
        # name of the reference in the path of the FASTQC output files.
        fqc_segs = path.FASTQC_DEMULT_SEGS
        fqc_vals = {path.TOP: out_dir, path.SAMP: sample, path.CMD: CMD_QC,
                    path.REF: fq_inp.ref}
    else:
        # Otherwise, use only the name of the sample, command, and step.
        fqc_segs = path.FASTQC_SEGS
        fqc_vals = {path.TOP: out_dir, path.SAMP: sample, path.CMD: CMD_QC}
    if fastqc:
        # Run FASTQC on the input FASTQ files.
        fqc_out = path.build(*fqc_segs, **fqc_vals, step=path.STEPS_QC[0])
        try:
            run_fastqc(fq_inp, fqc_out, extract=qc_extract, n_procs=n_procs)
        except Exception as error:
            logger.error(f"Failed to run FASTQC on {fq_inp}: {error}")
    if cut:
        # Trim adapters and low-quality bases with Cutadapt.
        fq_cut = fq_inp.to_new(path.StepSeg,
                               top=temp_dir,
                               step=path.STEPS_ALIGN[1])
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
            fqc_out = path.build(*fqc_segs, **fqc_vals, step=path.STEPS_QC[1])
            try:
                run_fastqc(fq_cut, fqc_out, extract=qc_extract, n_procs=n_procs)
            except Exception as error:
                logger.error(f"Failed to run FASTQC on {fq_inp}: {error}")
    else:
        fq_cut = None
    # Align the FASTQ to the reference sequence using Bowtie2.
    bam_whole = path.build(*path.XAM_STEP_SEGS,
                           top=temp_dir, sample=sample,
                           cmd=CMD_ALIGN, step=path.STEPS_ALIGN[4],
                           ref=refset, ext=path.BAM_EXT)
    reads_align = run_xamgen(fq_inp if fq_cut is None else fq_cut,
                             bam_whole,
                             index_pfx=bowtie2_index,
                             n_procs=n_procs,
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
                             bt2_orient=bt2_orient,
                             min_mapq=min_mapq)
    if not save_temp and fq_cut is not None:
        # Delete the trimmed FASTQ file (do not delete the input FASTQ
        # file even if trimming was not performed).
        for fq_file in fq_cut.paths.values():
            fq_file.unlink(missing_ok=True)
    # The number of reads after trimming is defined as the number fed to
    # Bowtie 2, regardless of whether the reads were actually trimmed.
    reads_trim = reads_align.pop("reads")
    # If the reads were trimmed, then the initial number must be found
    # by counting the reads in the input FASTQ. Otherwise, the initial
    # number equals the number after trimming.
    reads_init = fq_inp.n_reads if cut else reads_trim
    # Count the reads after filtering.
    flagstats = run_flagstat(bam_whole, None)
    paired_two, paired_one, singles = count_single_paired(flagstats)
    if fq_inp.paired:
        if singles:
            raise ValueError(f"{bam_whole} got {singles} single-end reads")
        reads_filter = {"paired-end, both mates mapped": paired_two,
                        "paired-end, one mate unmapped": paired_one}
    else:
        if n_paired := paired_two + paired_one:
            raise ValueError(f"{bam_whole} got {n_paired} paired-end reads")
        reads_filter = {"single-end": singles}
    # Index the sorted BAM file in order to split it by reference.
    bam_index = get_bam_index(bam_whole)
    index_bam(bam_whole, bam_index, n_procs=n_procs)
    # Link the file of reference sequences to the output directory.
    xams_out_dir = path.builddir(path.SampSeg, path.CmdSeg,
                                 top=out_dir, sample=sample, cmd=CMD_ALIGN)
    refs_file = xams_out_dir.joinpath(fasta.name)
    if not refs_file.is_file():
        refs_file.hardlink_to(fasta)
    # Split the indexed BAM file into one file per reference.
    xams_out: list[Path] = list()
    reads_refs = dict()
    for ref in refs:
        try:
            xam_ref = path.build(*path.XAM_SEGS,
                                 top=out_dir, sample=sample, cmd=CMD_ALIGN,
                                 ref=ref, ext=OUT_EXT)
            if xam_ref.parent != xams_out_dir:
                raise path.PathValueError(f"{xam_ref} is not in {xams_out_dir}")
            # Export the reads that align to the given reference.
            run_export(bam_whole, xam_ref, refs_file=refs_file, n_procs=n_procs)
            # Count the reads in the name-sorted file.
            n_reads = count_total_records(run_flagstat(xam_ref, None,
                                                       n_procs=n_procs))
            reads_refs[ref] = n_reads
            if n_reads >= min_reads:
                logger.debug(f"{ref} got {n_reads} reads")
                xams_out.append(xam_ref)
            else:
                logger.warning(f"{ref} got {n_reads} reads, which is less than "
                               f"the minimum ({min_reads}): deleting {xam_ref}")
                xam_ref.unlink(missing_ok=True)
        except Exception as error:
            logger.error(f"Failed to output {ref}: {error}")
    # The pre-split BAM file and its index are no longer needed.
    if not save_temp:
        bam_whole.unlink(missing_ok=True)
        bam_index.unlink(missing_ok=True)
    # Write a report to summarize the alignment.
    report = AlignReport(
        out_dir=out_dir,
        sample=sample,
        demultiplexed=fq_inp.ref is not None,
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
        bt2_unal=bt2_unal,
        bt2_score_min=bt2_score_min_loc if bt2_local else bt2_score_min_e2e,
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
        reads_init=reads_init,
        reads_trim=reads_trim,
        reads_align=reads_align,
        reads_filter=reads_filter,
        reads_refs=reads_refs,
    )
    report.save()
    # Return a list of name-sorted BAM files, each of which contains a
    # set of reads that all align to the same reference.
    return xams_out


def fqs_pipeline(fq_units: list[FastqUnit],
                 main_fasta: Path, *,
                 max_procs: int,
                 parallel: bool,
                 out_dir: Path,
                 temp_dir: Path,
                 save_temp: bool,
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
                logger.critical(
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
                                        step=path.STEPS_ALIGN[0],
                                        ref=refset)
                # Make its parent directory if it does not exist.
                main_index.parent.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Created directory: {main_index.parent}")
                # Build the Bowtie2 index.
                try:
                    build_bowtie2_index(main_fasta, main_index,
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
    # Generate binary alignment map (BAM) files.
    bam_lists = dispatch(fq_pipeline, max_procs, parallel, hybrid=True,
                         args=iter_args, kwargs=dict(out_dir=out_dir,
                                                     temp_dir=temp_dir,
                                                     save_temp=save_temp,
                                                     **kwargs))
    bams = list(chain(*bam_lists))
    if not save_temp:
        # Delete the temporary files.
        for ref_file, index_prefix in temp_fasta_paths.values():
            # Reference file
            ref_file.unlink(missing_ok=True)
            logger.debug(f"Deleted temporary reference file: {ref_file}")
            # Index files
            for index_file in get_bowtie2_index_paths(index_prefix):
                index_file.unlink(missing_ok=True)
                logger.debug(f"Deleted temporary index file: {index_file}")
    # Return the final binary alignment map (BAM) files.
    return bams


def figure_alignments(fq_units: list[FastqUnit], refs: set[str]):
    """ Return a `dict` of every expected alignment of a sample to a
    reference sequence. Check for and remove duplicates. """
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
            fq_refs = [fq_unit.ref]
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


def check_fqs_bams(alignments: dict[tuple[str, str], FastqUnit],
                   out_dir: Path):
    """ Return every FASTQ unit on which alignment must be run and every
    expected BAM file that already exists. """
    alignments_missing: dict[tuple[str, str], FastqUnit] = dict()
    bams_existing: list[Path] = list()
    for (sample, ref), fq_unit in alignments.items():
        # Determine the path of the BAM file expected to result from the
        # alignment of the sample to the reference.
        bam_expect = path.build(*path.XAM_SEGS,
                                top=out_dir, cmd=CMD_ALIGN,
                                sample=sample, ref=ref, ext=OUT_EXT)
        if bam_expect.is_file():
            # If the BAM file already exists, then add it to the dict of
            # BAM files that have already been aligned.
            bams_existing.append(bam_expect)
        else:
            # If at least one BAM file for a FASTQ unit does not exist,
            # then align the FASTQ.
            alignments_missing[sample, ref] = fq_unit
    return alignments_missing, bams_existing


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


def list_fqs_bams(fq_units: list[FastqUnit],
                  refs: set[str],
                  out_dir: Path):
    """ List every FASTQ that needs to be aligned and every expected BAM
    file that already exists. """
    # Determine all possible alignments of a sample and reference.
    alignments = figure_alignments(fq_units, refs)
    # Determine which alignments need to be / have already been run.
    alignments_missing, bams_existing = check_fqs_bams(alignments, out_dir)
    # Merge entries for each non-demultiplexed FASTQ.
    fqs_to_align = merge_nondemult_fqs(alignments_missing.values())
    return fqs_to_align, set(bams_existing)


def get_bam_files(fq_units: list[FastqUnit],
                  fasta: Path, *,
                  out_dir: Path,
                  rerun: bool,
                  **kwargs) -> list[Path]:
    """ Run the alignment pipeline and return a tuple of all BAM files
    from the pipeline. """
    if not fq_units:
        logger.warning("No FASTQ files or pairs of FASTQ files were given")
        return list()
    if rerun:
        # Rerun all alignments.
        fqs_to_align = fq_units
        bams = set()
    else:
        # Get the names of all reference sequences.
        refs = {ref for ref, _ in parse_fasta(fasta)}
        # Run only the alignments whose outputs do not yet exist.
        fqs_to_align, bams = list_fqs_bams(fq_units, refs, out_dir)
    if fqs_to_align:
        # Align all FASTQs that need to be aligned.
        bams_new = set(fqs_pipeline(fqs_to_align, fasta,
                                    out_dir=out_dir, **kwargs))
    else:
        logger.warning("All given FASTQ files have already been aligned")
        bams_new = set()
    # Merge the existing and new BAM paths into a tuple of strings.
    return list(bams | bams_new)
