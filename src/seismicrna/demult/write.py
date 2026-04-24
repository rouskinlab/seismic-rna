from typing import Iterable
import numpy as np
import gzip
import shutil
from collections import defaultdict
from itertools import zip_longest
from datetime import datetime
from pathlib import Path
from ..align.fqunit import FastqUnit, DuplicateSampleReferenceError, fastq_gz
from ..core.logs import logger
from ..core import path
from ..core.seq import parse_fasta
from ..core.extern.shell import args_to_cmd, run_cmd, SEQKIT_CMD
from ..core.task import dispatch, as_list_of_tuples
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.path import symlink_if_needed, mkdir_if_needed, FQ_ALL_EXTS
from ..core.seq import DNA
from .barcode import RefBarcodes

PART_STR = ".part_"
ALLOWED_SUFFIXES = sorted(
    (ext for ext in FQ_ALL_EXTS if ext.startswith('.')),
    key=len,
    reverse=True
)

def check_demult_fqs(demult_fqs: dict[tuple[str, str], FastqUnit],
                     out_dir: Path,
                     branches: dict[str, str]):
    """ Return every FASTQ unit on which demultiplexing must be run and every
    expected demultiplexed file that already exists. """
    fqs_missing: dict[tuple[str, str], FastqUnit] = dict()
    demult_fqs_extant: list[Path] = list()
    for (sample, ref), fq_unit in demult_fqs.items():
        # Determine the path of the demultiplexed file expected to result from the
        # demultiplexing of the sample by reference.
        if fq_unit.paired:
            fqs_found = 0
            for exts, segs in zip((path.FQ1_EXTS, path.FQ2_EXTS),
                                 (path.DMFASTQ1_SEGS, path.DMFASTQ2_SEGS)):
                for ext in exts:
                    demult_expect = path.build(segs,
                                                {path.TOP: out_dir,
                                                path.STEP: path.DEMULT_STEP,
                                                path.BRANCHES: branches,
                                                path.SAMPLE: sample,
                                                path.REF: ref,
                                                path.EXT: ext})
                    if demult_expect.is_file():
                        fqs_found += 1
                        # If the demultiplexed file already exists, then add it to the
                        # dict of demultiplexed files that have already been generated.
                        demult_fqs_extant.append(demult_expect)
                        break
            if fqs_found < 2:
                # If at least one demultiplexed file for a FASTQ unit does not exist,
                # then demultiplex the FASTQ.
                fqs_missing[sample, ref] = fq_unit
        else:
            for ext in path.FQ_EXTS:
                demult_expect = path.build(path.DmFastqSeg,
                                           {path.TOP: out_dir,
                                            path.STEP: path.DEMULT_STEP,
                                            path.BRANCHES: branches,
                                            path.SAMPLE: sample,
                                            path.REF: ref,
                                            path.EXT: ext})
                if demult_expect.is_file():
                    # If the demultiplexed file already exists, then add it to the
                    # dict of demultiplexed files that have already been generated.
                    demult_fqs_extant.append(demult_expect)
                    break
            else:
                # If at least one demultiplexed file for a FASTQ unit does not exist,
                # then demultiplex the FASTQ.
                fqs_missing[sample, ref] = fq_unit
    return fqs_missing, demult_fqs_extant


def merge_fqs(fq_units: Iterable[FastqUnit]):
    """ For every FASTQ that is not demultiplexed, merge all the keys
    that map to the FASTQ into one key: (sample, None). Merging ensures
    that every non-demultiplexed FASTQ is aligned only once to the whole
    set of references, not once for every reference in the set. This
    function is essentially the inverse of `figure_alignments`. """
    merged: dict[tuple[str, str | None], FastqUnit] = dict()
    for fq_unit in fq_units:
        merged[fq_unit.sample, fq_unit.ref] = fq_unit
    return list(merged.values())


def list_demulted_fqs(demult_fqs: dict[tuple[str, str], FastqUnit],
                      out_dir: Path,
                      branches: dict[str, str]):
    """ List every FASTQ to demultiplex and every extant demultiplexed file. """
    # Determine which FASTQs need to be or have already been run.
    fqs_missing, demult_fqs_extant = check_demult_fqs(demult_fqs, out_dir, branches)
    fqs_merged = merge_fqs(fqs_missing.values())
    return fqs_merged, set(demult_fqs_extant)


def list_demult(fq_units: list[FastqUnit], refs: set[str]):
    """ List every expected demultiplexed FASTQ from a multiplexed FASTQ. """
    logger.routine("Began listing demultiplexed FASTQs")
    # Map each combination of a sample and reference to a FASTQ unit.
    demult_fqs: dict[tuple[str, str], FastqUnit] = dict()
    for fq_unit in fq_units:
        fq_refs = refs
        # Add each sample-reference pair to the expected demultiplexed fastqs.
        sample = fq_unit.sample
        for ref in fq_refs:
            logger.detail(
                f"Adding demultiplexed FASTQ {repr(ref)} for sample {repr(sample)}"
            )
            sample_ref = sample, ref
            if sample_ref in demult_fqs:
                raise DuplicateSampleReferenceError(sample_ref)
            demult_fqs[sample_ref] = fq_unit
    logger.routine("Ended listing demultiplexed FASTQs")
    return demult_fqs


def demult_samples(fq_units: list[FastqUnit],
                   fasta: Path, *,
                   refs_meta: Path,
                   barcode_start: int,
                   barcode_end: int,
                   mismatch_tolerance: int,
                   index_tolerance: int,
                   allow_n: bool,
                   read_pos: int | None,
                   barcode: tuple[tuple[str, DNA, int]],
                   out_dir: Path,
                   branch: str,
                   force: bool,
                   **kwargs) -> list[Path]:
    """ Run the demult pipeline and return a tuple of all fastq files
    from the pipeline. """
    if not fq_units:
        logger.detail("No FASTQ files or pairs of files were given to demult")
        return list()
    # Even though the ancestors argument of path.add_branch() is empty,
    # use it instead of just "branches = {path.DEMULT_STEP: branch}" in
    # order to validate the branch name.
    branches = path.add_branch(path.ALIGN_STEP, branch, dict())
    # List the names of all reference sequences.
    refs = set(parse_fasta(fasta, None))

    if read_pos is None:
        read_pos = barcode_start

    if barcode_start and barcode_end:
        coords = [(ref, barcode_start, barcode_end, read_pos) for ref in refs]
    elif barcode_start or barcode_end:
        if barcode_start:
            raise ValueError("Missing --barcode-end")
        if barcode_end:
            raise ValueError("Missing --barcode-start")
    else:
        coords = []

    if refs_meta:
        barcodes = RefBarcodes(ref_seqs=parse_fasta(fasta, DNA),
                               refs_meta_file=refs_meta,
                               coords=coords,
                               bcs=barcode,
                               mismatches=mismatch_tolerance,
                               index_tolerance=index_tolerance,
                               allow_n=allow_n)
    else:
        barcodes = RefBarcodes(ref_seqs=parse_fasta(fasta, DNA),
                               coords=coords,
                               bcs=barcode,
                               mismatches=mismatch_tolerance,
                               index_tolerance=index_tolerance,
                               allow_n=allow_n)

    # List all demultiplexed FASTQs and check for duplicates.
    demult_fqs = list_demult(fq_units, refs | barcodes.uniq_names)

    if force:
        # force all alignments.
        fqs_to_demult = fq_units
        demult_fqs_extant = set()
    else:
        # Demultiplex only the fastqs whose outputs do not yet exist.
        fqs_to_demult, demult_fqs_extant = list_demulted_fqs(demult_fqs, out_dir, branches)

    if fqs_to_demult:
        # Demultiplex all FASTQs that need to be demultiplexed.
        fq_dirs_new = set(demult_fqs_pipeline(fqs_to_demult,
                                              fasta=fasta,
                                              barcodes=barcodes,
                                              mismatch_tolerance=mismatch_tolerance,
                                              index_tolerance=index_tolerance,
                                              out_dir=out_dir,
                                              branches=branches,
                                              **kwargs))
    else:
        logger.warning("All given FASTQ files have already been demultiplexed: "
                       "use --force to overwrite")
        fq_dirs_new = set()

    # Merge the existing and new FASTQ paths into a tuple of strings.
    return list({fq.parent for fq in demult_fqs_extant} | fq_dirs_new)


def to_range(pos: int, tolerance: int):
    if tolerance > 0:
        range = f"{pos-tolerance},{pos+tolerance}"
    else:
        range = str(pos)
    return range


def demult_fqs_pipeline(fq_units: list[FastqUnit],
                        fasta: Path,
                        barcodes: RefBarcodes,
                        index_tolerance: int,
                        mismatch_tolerance: int, *,
                        num_cpus: int,
                        tmp_dir: Path,
                        **kwargs) -> list[Path]:
    """ Run all stages of demultiplexing for one or more FASTQ files or pairs
    of mated FASTQ files. """
    logger.routine("Began running the demultiplexing pipeline")
    # Validate the maximum number of processes.
    if num_cpus < 1:
        logger.warning("--num_cpus must be ≥ 1: setting to 1")
        num_cpus = 1
    # Make the arguments for each demultiplexing task.
    iter_args: list[tuple[FastqUnit, str, str]] = list()

    for fq_unit in fq_units:
            iter_args.append((fq_unit, barcodes))
            # Add these arguments to the lists of arguments that
            # will be passed to fq_pipeline.
    # Generate demultiplexed FASTQ files.
    fq_dirs = dispatch(demult_fq_pipeline,
                       num_cpus=num_cpus,
                       pass_num_cpus=True,
                       as_list=True,
                       ordered=False,
                       raise_on_error=False,
                       args=iter_args,
                       kwargs=dict(tmp_dir=tmp_dir, force_serial=True, **kwargs))
    logger.routine("Ended running the demultiplexing pipeline")
    # Return the final demultiplexed FASTQ directories.
    return fq_dirs


def demult_fq_pipeline(fq_inp: FastqUnit,
                       barcodes: RefBarcodes,
                       *,
                       out_dir: Path,
                       tmp_dir: Path,
                       num_cpus: int = 1,
                       **kwargs):
    """ Run all stages of the demult pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    began = datetime.now()
    release_dir, working_dir = get_release_working_dirs(tmp_dir)
    logger.routine(f"Began processing {fq_inp} through the demult pipeline")
    # Get attributes of the sample and references.

    num_split = num_cpus-1
    split_fqs = split_fq(fq_inp, working_dir, num_split=num_split)

    args = [(fq_inp, fqs, barcodes) for fqs in split_fqs]
    demult_parts = dispatch(process_fq_part,
                            num_cpus=num_cpus,
                            pass_num_cpus=False,
                            as_list=True,
                            ordered=True,
                            raise_on_error=False,
                            args=args,
                            kwargs=dict(release_dir=release_dir,
                                        **kwargs))

    demult_parts = np.array(demult_parts)

    assembled_parts = list()
    n_mates = demult_parts.shape[1]
    n_parts = demult_parts.shape[2]

    for mate in range(n_mates):
        for part in range(n_parts):
            parts = [fq_mates[mate][part] for fq_mates in demult_parts]
            assembled_parts.append(parts)

    arg_parts = [arg_part for arg_part in assembled_parts if any(fq.exists() for fq in arg_part)]
    dispatch(merge_parts,
             num_cpus=num_cpus,
             pass_num_cpus=False,
             as_list=True,
             ordered=False,
             raise_on_error=False,
             args=as_list_of_tuples(arg_parts))
    out_path = release_to_out(out_dir=out_dir, release_dir=release_dir, initial_path=release_dir/Path(fq_inp.sample)/"demult")

    return out_path


def split_fq(fq_inp: FastqUnit,
             working_dir: Path,
             num_split: int):
    split_dir = working_dir/fq_inp.sample/"split"
    for fq in fq_inp.paths.values():
        fq = Path(fq)
        num_split = max(num_split, 1)
        if num_split == 1:
            split_path = (split_dir / strip_all_fq_suffixes(fq.name)).with_suffix(".part_001"+"".join(fq.suffixes))
            mkdir_if_needed(split_path.parent)
            # Rather than "splitting" into 1 file, symlink the original file.
            symlink_if_needed(split_path, fq)
        else:
            args = [
                    SEQKIT_CMD,
                    "split2",
                    "-p", str(num_split),
                    "-j", str(num_split),
                    str(fq),
                    "-O", split_dir
                ]
            cmd = args_to_cmd(args)
            run_cmd(cmd)
    return get_split_paths(split_dir, fq_inp, num_split)


def strip_all_fq_suffixes(path: str | Path):
    name = Path(path).name
    for ext in ALLOWED_SUFFIXES:
        if name.endswith(ext):
            name = name[:-len(ext)]
            # Aditionally strip the .part_ suffix from split fastqs
            suffix = Path(name).suffix
            if suffix.startswith(PART_STR):
                name = name[:-len(suffix)]
            return name


def rename_fq_part(fq_path: Path) -> Path:
    suffixes = fq_path.suffixes
    base_name = strip_all_fq_suffixes(fq_path)
    part = next((s for s in suffixes if s.startswith(PART_STR)), None)
    if part is None:
        return fq_path
    final_suffixes = [s for s in suffixes if s != part]
    new_name = f"{part[len(PART_STR):]}_{base_name}" + "".join(final_suffixes)
    new_fq = fq_path.with_name(new_name)
    fq_path.rename(new_fq)
    return new_fq


def get_split_paths(split_dir, fq_inp, num_parts):
    fq_parts = []
    for inp_path in fq_inp.paths.values():
        p = Path(inp_path)
        suffix = get_fq_suffix(p)
        base_name = strip_all_fq_suffixes(p)
        fq_parts.append((base_name, suffix))
    split_paths = [
        tuple(split_dir / f"{stem}.part_{i:03d}{ext}" for stem, ext in fq_parts)
        for i in range(1, num_parts + 1)
    ]
    return [
        tuple(rename_fq_part(fq) for fq in fq_tuple)
        for fq_tuple in split_paths
    ]


def get_part(fq_path: Path):
    return fq_path.name.split("_")[0]


def remove_suffixes(path: Path):
    while path.suffix:
        path = path.with_suffix("")
    return path


def get_fq_suffix(path: Path):
    name = path.name
    for ext in sorted(ALLOWED_SUFFIXES, key=len, reverse=True):
        if name.endswith(ext):
            return ext
    return ""

def process_fq_part(fq_inp: FastqUnit,
                    fqs: tuple[Path, ...],
                    barcodes: RefBarcodes,
                    *,
                    release_dir: Path,
                    branches: dict[str, str],
                    **kwargs):
    # Define segment patterns and corresponding mappings.
    fq_seg_patterns = [path.Fastq1Seg, path.Fastq2Seg, path.FastqSeg]
    fq_seg_map = {
        path.Fastq1Seg: path.DMFASTQ1_SEGS,
        path.Fastq2Seg: path.DMFASTQ2_SEGS,
        path.FastqSeg:  path.DMFASTQ_SEGS,
    }

    # Build FastqUnit using the provided FASTQ paths.
    fq_unit_args = dict(zip(fq_inp.paths.keys(), fqs))
    fq_unit_part = FastqUnit(**fq_unit_args,
                             phred_enc=fq_inp.phred_enc,
                             one_ref=fq_inp.one_ref)

    # Prepare output paths mapping for each barcode.
    out_fq_unit_args = defaultdict(dict)
    for fq_type, fq_path in fq_unit_part.paths.items():
        # Get the complete extension (e.g. '.fastq.gz')
        fq_extension = "".join(fq_path.suffixes)
        seg_pattern = next((pattern for pattern in fq_seg_patterns
                            if path.path_matches(fq_path.name, [pattern])), None)
        if not seg_pattern:
            logger.error(f"{fq_path} does not match any known FASTQ pattern.")
            continue

        dm_fastq_segs = fq_seg_map[seg_pattern]
        # Process only the first three barcode items.
        for name in barcodes.uniq_names:
            tmp_fq_path = path.buildpar(
                dm_fastq_segs,
                {
                    path.TOP: release_dir,
                    path.SAMPLE: fq_inp.sample,
                    path.STEP: path.DEMULT_STEP,
                    path.BRANCHES: branches,
                    path.REF: f"{get_part(fq_path)}_{name}",
                    path.EXT: dm_fastq_segs[-1].exts[0],
                }
            )
            # Ensure the output file retains the original suffix (e.g., gzipped)
            out_fq_path = remove_suffixes(tmp_fq_path).with_suffix(fq_extension)
            out_fq_unit_args[name][fq_type] = out_fq_path

    # Build output FastqUnit objects per barcode.
    out_fqs = {
            name: FastqUnit(**paths,
                            phred_enc=fq_inp.phred_enc,
                            one_ref=fq_inp.one_ref)
            for name, paths in out_fq_unit_args.items()
    }

    demult_fqs = demult_ahocorasick(fq_unit_part, out_fqs, barcodes)

    return demult_fqs


def get_open_func(fq_path: Path):
    return gzip.open if fastq_gz(fq_path) else open


def check_matches(matches: Iterable[tuple[tuple[int, str, set]]], barcodes):
    samples = set()
    for match in matches:
        for read in match:
            if read is not None:
                end_index, barcode_id = read
                name = barcodes.name_map[barcode_id]
                valid_ends = barcodes.valid_positions[barcode_id]
                if end_index in valid_ends:
                    samples.add(name)
    if not samples:
        return None

    if len(samples) > 1:
        logger.warning(f"Read matched more than one sample {samples}")
        # raise ValueError(f"Read matched more than one sample {samples}")
        return None
    else:
        return samples.pop()


def demult_ahocorasick(fq_unit: FastqUnit,
                       out_fqs: dict[int, FastqUnit],
                       barcodes: RefBarcodes,
                       buffer_limit: int = 1000):
    read_buffer = defaultdict(list)

    # Choose the appropriate open function based on file suffix.
    open_funcs = dict()

    for fq in fq_unit.paths.values():
        if not fq.exists(): # TODO Make sure this works
            return (tuple(),)
        open_func = gzip.open if fastq_gz(fq) else open
        open_funcs[fq] = open_func

    automaton = barcodes.automaton
    rc_automaton = barcodes.rc_automaton
    # Process FASTQ records using FastqUnit.iter_records()
    count = 0
    total = 0
    for segs, recs in fq_unit.iter_records(segments=[barcodes.read_pos_range, barcodes.rc_read_pos_range]):
        seqs = [" ".join(seg) for seg in segs]
        bcs = [(match, match_rc) for match, match_rc in zip_longest(automaton.iter(seqs[0]), rc_automaton.iter(seqs[1]))]
        name = check_matches(bcs, barcodes)
        if name:
            for fq_lines, fq_path in zip(recs, out_fqs[name].paths.values()):
                read_buffer[fq_path].append(fq_lines)
                # Flush if the buffer for this barcode is large.
                if len(read_buffer[fq_path]) >= buffer_limit:
                    with open_func(fq_path, "at") as out_file:
                        out_file.write("\n".join("\n".join(record) for record in read_buffer[fq_path]) + "\n")
                    read_buffer[fq_path].clear()
            count += 1
        else:
            pass
        total += 1
    # Write remaining buffered records.
    for fq_path, records in read_buffer.items():
        if records:
            with open_func(fq_path, "at") as out_file:
                out_file.write("\n".join("\n".join(record) for record in records) + "\n")

    logger.detail(f"Identified {count} out of {total} reads ({100*(count/total):.3f}%)")

    return tuple((tuple(out_fqs[name].paths.values()) for name in barcodes.uniq_names))


def merge_parts(parts = list[Path]):
    parts = [part for part in parts if part.exists()]
    if len(parts) == 0:
        logger.warning("No files to merge")
        return None
    base_name = parts[0].name.split('_', 1)[1]
    if not all(
        ('_' in part.name and part.name.split('_', 1)[1] == base_name)
        for part in parts
    ):
        logger.error("To merge parts, all files must have the same base name.")
    sorted_parts = sorted(parts, key=lambda f: int(f.name.split('_', 1)[0]))

    merged_path = parts[0].with_name(base_name)
    with merged_path.open("wb") as merged:
        for part in sorted_parts:
            with part.open("rb") as part_handle:
                shutil.copyfileobj(part_handle, merged)
            part.unlink()

    return merged_path
