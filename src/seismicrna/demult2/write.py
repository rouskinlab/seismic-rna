from typing import Iterable
import numpy as np
import gzip
import shutil

from datetime import datetime
from pathlib import Path
from ..align.fqunit import FastqUnit, DuplicateSampleReferenceError
from ..core.logs import logger
from ..core import path
from ..core.seq import parse_fasta
from ..core.extern.shell import args_to_cmd, run_cmd, SEQKIT_CMD, UMI_TOOLS_CMD, UMI_TOOLS_EXTRACT_CMD
from ..core.task import dispatch, as_list_of_tuples
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.seq import DNA
from .barcode import RefBarcodes

PART_STR = ".part_"

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
                                 (path.DMFASTQ1_SEGS2, path.DMFASTQ2_SEGS2)):
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
                               bcs=barcode)

    # List all demultiplexed FASTQs and check for duplicates.
    demult_fqs = list_demult(fq_units, refs | barcodes.names)

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
                        max_procs: int,
                        tmp_dir: Path,
                        **kwargs) -> list[Path]:
    """ Run all stages of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
    logger.routine("Began running the demultiplexing pipeline")
    # Validate the maximum number of processes.
    if max_procs < 1:
        logger.warning("max_procs must be ≥ 1: setting to 1")
        max_procs = 1
    # Make the arguments for each demultiplexing task.
    iter_args: list[tuple[FastqUnit, str, str]] = list()
    bcs = barcodes.barcodes
    rc_bcs = barcodes.rc_barcodes
    bcs = list(map(str, bcs))
    rc_bcs = list(map(str, rc_bcs))

    sub_patterns = list()

    bc_len = len(bcs[0])

    for read_pos, bc_list in barcodes.by_pos.items():
        bc_str = "|".join(map(str, bc_list))
        sub_patterns.append(f"(.{{{to_range(read_pos, index_tolerance)}}}(?P<umi_1>{bc_str}){{s<={mismatch_tolerance}}})")
        # sub_patterns.append(f"(.{{{to_range(read_pos, index_tolerance)}}}(?P<umi_1>.{{{bc_len}}}){{s<={mismatch_tolerance}}})")
    for rc_read_pos, rc_bc_list in barcodes.rc_by_pos.items():
        rc_bc_str = "|".join(map(str, rc_bc_list))
        sub_patterns.append(f"(.{{{to_range(rc_read_pos, index_tolerance)}}}(?P<umi_1>{rc_bc_str}){{s<={mismatch_tolerance}}})")
        # sub_patterns.append(f"(.{{{to_range(rc_read_pos, index_tolerance)}}}(?P<umi_1>.{{{bc_len}}}){{s<={mismatch_tolerance}}})")

    bc_pattern = "|".join(sub_patterns)
    bc_pattern2 = bc_pattern # TODO: Implement different barcodes for each read
    # One demultiplexing task will be created for each FASTQ unit.
    for fq_unit in fq_units:
            iter_args.append((fq_unit, bc_pattern, bc_pattern2, barcodes))
            # Add these arguments to the lists of arguments that
            # will be passed to fq_pipeline.
    # Generate demultiplexed FASTQ files.
    fq_dirs = dispatch(demult_fq_pipeline,
                       max_procs,
                       args=iter_args,
                       kwargs=dict(tmp_dir=tmp_dir, **kwargs))
    logger.routine("Ended running the demultiplexing pipeline")
    # Return the final demultiplexed FASTQ directories.
    return fq_dirs


def demult_fq_pipeline(fq_inp: FastqUnit,
                       bc_pattern: str,
                       bc_pattern2: str,
                       barcodes: RefBarcodes,
                       *,
                       out_dir: Path,
                       tmp_dir: Path,
                    #    batch_size: int,
                       n_procs: int = 1,
                       **kwargs):
    """ Run all stages of the demult pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    began = datetime.now()
    release_dir, working_dir = get_release_working_dirs(tmp_dir)
    logger.routine(f"Began processing {fq_inp} through the demult pipeline")
    # Get attributes of the sample and references.

    # fq1, fq2 = fq_inp.paths.values()
    # run_ahocorasick(fq_path=fq1, fq2_path=fq2, barcodes=barcodes)
    # return "Done"
    split_fqs = split_fq(fq_inp, working_dir, batch_size=10000, n_procs=n_procs)

    barcodes = barcodes.as_dict
    args = [(fq_inp, fqs, bc_pattern, bc_pattern2, barcodes) for fqs in split_fqs]
    demult_parts = dispatch(process_fq_part,
                            n_procs,
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

    dispatch(merge_parts,
             n_procs,
             args=as_list_of_tuples(assembled_parts),
             pass_n_procs=False)

    out_path = release_to_out(out_dir=out_dir, release_dir=release_dir, initial_path=release_dir/Path(fq_inp.sample)/"demult")

    return out_path


def split_fq(fq_inp: FastqUnit,
             working_dir: Path,
             batch_size: int,
             n_procs: int):
    for fq in fq_inp.paths.values():
        split_dir = working_dir/fq_inp.sample/"split"
        args = [
                SEQKIT_CMD,
                "split2",
                "-p", str(n_procs), #str(min(n_procs, fq_inp.n_reads/batch_size)),
                "-j", str(n_procs),
                str(fq),
                "-O", split_dir
            ]
        cmd = args_to_cmd(args)
        run_cmd(cmd)
    return get_split_paths(split_dir, fq_inp, n_procs)


def rename_fq_part(fq_path: Path) -> Path:
    suffixes = fq_path.suffixes
    base_name = fq_path.name
    for s in reversed(suffixes):
        if base_name.endswith(s):
            base_name = base_name[:-len(s)]
    part = next((s for s in suffixes if s.startswith(PART_STR)), None)
    if part is None:
        return fq_path
    final_suffixes = [s for s in suffixes if s != part]
    new_name = f"{part[len(PART_STR):]}_{base_name}" + "".join(final_suffixes)
    new_fq = fq_path.with_name(new_name)
    fq_path.rename(new_fq)
    return new_fq


def get_split_paths(split_dir, fq_inp, num_parts):
    fq_parts = [
        (Path(p).name.rsplit(".", maxsplit=2)[0], "".join(Path(p).suffixes)) # Can FASTQs processed by SEISMIC have periods in their names?
        for p in fq_inp.paths.values()
    ]
    split_paths = [
        tuple(split_dir / f"{stem}.part_{i:03d}{ext}" for stem, ext in fq_parts)
        for i in range(1, num_parts + 1)
    ]
    return [
        tuple(rename_fq_part(fq) for fq in fq_tuple)
        for fq_tuple in split_paths
    ]


def run_ahocorasick(fq_path: Path,
                    fq2_path: Path,
                    barcodes: RefBarcodes):

    if fq_path.suffix == ".gz":
        open_func = gzip.open
    else:
        open_func = open

    with open_func(fq_path, "rt") as fq1:
        with open_func(fq2_path, "rt") as fq2:
            n = 0
            while True:
                # print(n)
                fq1_lines = [fq1.readline() for _ in range(4)]
                fq2_lines = [fq2.readline() for _ in range(4)]
                if not fq1_lines[0]:
                    break
                header = fq1_lines[0].strip()
                name = header.split()[0]
                seq = fq1_lines[1] + fq2_lines[1]
                for end_index, name in barcodes.automaton.iter(seq):
                    pass
                    # print((end_index, name))
                    # logger.error((start_index, end_index, (insert_order, original_value)))
                n+=1



def extract_barcodes(fq_paths: tuple[Path, ...],
                     bc_pattern: str,
                     bc_pattern2: str):
    """
    Run umi_tools extract on the given FASTQ(s).
    """
    ext_fq_paths = tuple([fq.parent.parent/"extracted"/fq.name for fq in fq_paths]) # TODO: Improve path logic
    path.mkdir_if_needed(ext_fq_paths[0].parent)
    args = [
        UMI_TOOLS_CMD,
        UMI_TOOLS_EXTRACT_CMD,
        "--extract-method=regex",
        "--either-read",
        f"--bc-pattern='{bc_pattern}'",
        "-I", fq_paths[0],
        "-S", ext_fq_paths[0]
        ]

    if len(fq_paths) > 1:
        assert(len(fq_paths) == 2)
        args.extend([
        f"--bc-pattern2='{bc_pattern2}'",
        "--read2-in", fq_paths[1],
        "--read2-out", ext_fq_paths[1]
    ])
    logger.action(f"Extracting barcodes from {fq_paths}")
    args = list(map(str, args))
    cmd = " ".join(args)
    cmd = args

    run_cmd(cmd, shell=False)

    # Clean up input paths to save disk space.
    for fq in fq_paths:
        fq.unlink()
    return ext_fq_paths

def get_part(fq_path: Path):
    return fq_path.name.split("_")[0]


def remove_suffixes(path: Path):
    while path.suffix:
        path = path.with_suffix("")
    return path


def process_fq_part(fq_inp: FastqUnit,
                    fqs: tuple[Path, ...],
                    bc_pattern: str,
                    bc_pattern2: str,
                    barcodes: dict[str, DNA],
                    *,
                    release_dir: Path,
                    branches: dict[str, str],
                    n_procs: int,
                    **kwargs):

    ext_fqs = extract_barcodes(fqs,
                               bc_pattern,
                               bc_pattern2)

    fq_seg_patterns = [path.Fastq1Seg,
                       path.Fastq2Seg,
                       path.FastqSeg]
    fq_seg_map = {path.Fastq1Seg: path.DMFASTQ1_SEGS2,
                  path.Fastq2Seg: path.DMFASTQ2_SEGS2,
                  path.FastqSeg: path.DMFASTQ_SEGS2}

    args = list()
    for fq_path in ext_fqs:
        fq_extension = "".join(fq_path.suffixes)
        seg_pattern = next((pattern for pattern in fq_seg_patterns if path.path_matches(fq_path.name, [pattern])), None)
        if not seg_pattern:
            logger.error(f"{fq_path} does not match any known FASTQ pattern.")
        out_paths = dict()
        for ref, barcode in barcodes.items():
            dm_fastq_segs = fq_seg_map[seg_pattern]
            tmp_fq_path = path.buildpar(dm_fastq_segs,
                                        {path.TOP: release_dir,
                                         path.SAMPLE:fq_inp.sample,
                                         path.STEP: path.DEMULT_STEP,
                                         path.BRANCHES: branches,
                                         path.REF: f"{get_part(fq_path)}_{ref}",
                                         path.EXT: dm_fastq_segs[-1].exts[0]})
            # Ensure suffix is the same as input suffix (in case gzipped).
            out_fq_path = remove_suffixes(tmp_fq_path).with_suffix(fq_extension)
            out_paths[barcode] = out_fq_path
        args.append((fq_path, out_paths))
    demult_fqs = dispatch(demult_fq,
                          n_procs,
                          args=args,
                          pass_n_procs=False,
                          kwargs=kwargs)
    return demult_fqs

def demult_fq(fq_path: Path,
              out_paths: dict[str, DNA],
              **kwargs):

    if fq_path.suffix == ".gz":
        open_func = gzip.open
    else:
        open_func = open

    out_streams = dict()
    for barcode in out_paths:
        out_streams[barcode] = open_func(out_paths[barcode], "wt")
    with open_func(fq_path, "rt") as fq:
        while True:
            fq_lines = [fq.readline() for _ in range(4)]
            if not fq_lines[0]:
                break
            header = fq_lines[0].strip()
            name = header.split()[0]
            if "_" in name:
                read_bc = DNA(name.rsplit("_", 1)[1])
            else:
                read_bc = "unknown"
            if read_bc in out_streams:
                out_streams[read_bc].writelines(fq_lines)
            elif read_bc.rc in out_streams:
                out_streams[read_bc.rc].writelines(fq_lines)
            else:
                pass
                # logger.warning(f"Unknown barcode encountered {read_bc}")
    for handle in out_streams.values():
        handle.close()
    return list(out_paths.values())

def merge_parts(parts = list[Path]):

    base_name = parts[0].name.split('_', 1)[1]
    if not all(
        ('_' in part.name and part.name.split('_', 1)[1] == base_name)
        for part in parts
    ):
        raise ValueError("To merge parts, all files must have the same base name.")
    sorted_parts = sorted(parts, key=lambda f: int(f.name.split('_', 1)[0]))

    merged_path = parts[0].with_name(base_name)

    with merged_path.open("wb") as merged:
        for part in sorted_parts:
            with part.open("rb") as part_handle:
                shutil.copyfileobj(part_handle, merged)
            part.unlink()

    return merged_path