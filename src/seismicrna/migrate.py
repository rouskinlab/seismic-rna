import gzip
import os
import pickle
import shutil
from hashlib import md5
from pathlib import Path
from typing import Iterable

import brotli
from click import command

from .cluster.report import ClusterReport
from .core import path
from .core.arg import (CMD_MIGRATE,
                       arg_input_path,
                       opt_out_dir,
                       opt_num_cpus)
from .core.io import BadChecksumError, save_brickle
from .core.report import Report, NumBatchesF, ChecksumsF, RefseqChecksumF
from .core.run import run_func
from .core.task import as_list_of_tuples, dispatch
from .mask.report import MaskReport
from .relate.io import RelateBatchIO
from .relate.report import RelateReport

JSON_INDENT = " " * 4
JSON_DELIM = f",\n{JSON_INDENT}"


class FindAndReplaceError(ValueError):
    pass


def find_and_replace(file: Path, find: str, replace: str, count: int = 1):
    if file.suffix.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file, "rt") as f:
        text = f.read()
    find_count = text.count(find)
    if 0 <= count != find_count:
        raise FindAndReplaceError(
            f"{file} contains {repr(find)} {find_count} time(s)"
        )
    text = text.replace(find, replace)
    with open_func(file, "wt") as f:
        f.write(text)
    return find_count


def update_brickle(file: str | Path, md5_checksum: str):
    # Load the item from the file.
    with open(file, "rb") as f:
        data = f.read()
    md5_digest = md5(data).hexdigest()
    if md5_digest != md5_checksum:
        raise BadChecksumError(
            f"Expected MD5 digest of {file} to be {md5_checksum}, "
            f"but got {md5_digest}"
        )
    item = pickle.loads(brotli.decompress(data))
    # Edit the item's attributes.
    setattr(item, "branches", dict())
    if isinstance(item, RelateBatchIO):
        # Rename self._end5s, self._end3s, and self._muts
        for end in [5, 3]:
            setattr(item, f"seg_end{end}s", getattr(item, f"_end{end}s"))
            delattr(item, f"_end{end}s")
        setattr(item, "muts", getattr(item, "_muts"))
        delattr(item, "_muts")
    # Save the updated item back to the same file.
    return save_brickle(item, file, brotli_level=10, force=True)


def update_brickles(report: Report, top: str | Path):
    report_file = report.get_path(top)
    num_batches = report.get_field(NumBatchesF)
    checksums = report.get_field(ChecksumsF)
    if isinstance(report, RelateReport):
        if "readnames" in checksums:
            for batch in range(num_batches):
                batch_file = report_file.parent.joinpath(
                    f"readnames-batch-{batch}.brickle"
                )
                md5_checksum = checksums["readnames"][batch]
                checksums["readnames"][batch] = update_brickle(batch_file,
                                                               md5_checksum)
        report_type = "relate"
    elif isinstance(report, MaskReport):
        report_type = "mask"
    elif isinstance(report, ClusterReport):
        report_type = "cluster"
    else:
        raise TypeError(report)
    for batch in range(num_batches):
        batch_file = report_file.parent.joinpath(
            f"{report_type}-batch-{batch}.brickle"
        )
        md5_checksum = checksums[report_type][batch]
        checksums[report_type][batch] = update_brickle(batch_file,
                                                       md5_checksum)
    return checksums


def migrate_align_dir(align_dir: Path):
    for file in align_dir.iterdir():
        if file.name.endswith("align-report.json"):
            find_and_replace(file,
                             '"Branches": [],',
                             '"Branches": {},')


def migrate_relate_ref_dir(ref_dir: Path):
    report_file = ref_dir.joinpath("relate-report.json")
    find_and_replace(report_file,
                     '"Branches": [],',
                     '"Branches": {},')
    try:
        find_and_replace(report_file,
                         "MD5 checksums of batches",
                         "Checksums of batches (SHA-512)")
    except FindAndReplaceError:
        # The report is pooled.
        return
    find_and_replace(report_file,
                     "MD5 checksum of reference sequence",
                     "Checksum of reference sequence (SHA-512)")
    report = RelateReport.load(report_file)
    top, _ = report.parse_path(report_file)
    sha512_checksum = update_brickle(report_file.with_name("refseq.brickle"),
                                     report.get_field(RefseqChecksumF))
    setattr(report, RefseqChecksumF.key, sha512_checksum)
    sha512_checksums = update_brickles(report, top)
    setattr(report, ChecksumsF.key, sha512_checksums)
    report.save(top, force=True)


def migrate_mask_reg_dir(reg_dir: Path):
    report_file = reg_dir.joinpath("mask-report.json")
    find_and_replace(report_file,
                     '"Branches": [],',
                     '"Branches": {},')
    try:
        find_and_replace(report_file,
                         "MD5 checksums of batches",
                         "Checksums of batches (SHA-512)")
    except FindAndReplaceError:
        # The report is joined.
        return
    report = MaskReport.load(report_file)
    top, _ = report.parse_path(report_file)
    sha512_checksums = update_brickles(report, top)
    setattr(report, ChecksumsF.key, sha512_checksums)
    report.save(top, force=True)


def migrate_mask_ref_dir(ref_dir: Path, num_cpus: int):
    return dispatch(migrate_mask_reg_dir,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    ordered=False,
                    raise_on_error=False,
                    as_list=True,
                    args=as_list_of_tuples(ref_dir.iterdir()))


def migrate_cluster_reg_dir(reg_dir: Path):
    report_file = reg_dir.joinpath("cluster-report.json")
    find_and_replace(report_file,
                     '"Branches": [],',
                     '"Branches": {},')
    try:
        find_and_replace(report_file,
                         "MD5 checksums of batches",
                         "Checksums of batches (SHA-512)")
    except FindAndReplaceError:
        # The report is joined.
        return
    report = ClusterReport.load(report_file)
    top, _ = report.parse_path(report_file)
    sha512_checksums = update_brickles(report, top)
    setattr(report, ChecksumsF.key, sha512_checksums)
    report.save(top, force=True)


def migrate_cluster_ref_dir(ref_dir: Path, num_cpus: int):
    return dispatch(migrate_cluster_reg_dir,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    ordered=False,
                    raise_on_error=False,
                    as_list=True,
                    args=as_list_of_tuples(ref_dir.iterdir()))


def migrate_fold_reg_dir(reg_dir: Path):
    for file in reg_dir.iterdir():
        if file.name.endswith("fold-report.json"):
            find_and_replace(file,
                             '"Branches": [],',
                             '"Branches": {},')


def migrate_fold_ref_dir(ref_dir: Path, num_cpus: int):
    return dispatch(migrate_fold_reg_dir,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    ordered=False,
                    raise_on_error=False,
                    as_list=True,
                    args=as_list_of_tuples(ref_dir.iterdir()))


def migrate_sample_dir(sample_dir: Path, num_cpus: int):
    align_dir = sample_dir.joinpath("align")
    if align_dir.is_dir():
        migrate_align_dir(align_dir)
    relate_dir = sample_dir.joinpath("relate")
    if relate_dir.is_dir():
        dispatch(migrate_relate_ref_dir,
                 num_cpus=num_cpus,
                 pass_num_cpus=False,
                 ordered=False,
                 raise_on_error=False,
                 as_list=True,
                 args=as_list_of_tuples(relate_dir.iterdir()))
    mask_dir = sample_dir.joinpath("mask")
    if mask_dir.is_dir():
        dispatch(migrate_mask_ref_dir,
                 num_cpus=num_cpus,
                 pass_num_cpus=True,
                 ordered=False,
                 raise_on_error=False,
                 as_list=True,
                 args=as_list_of_tuples(mask_dir.iterdir()))
    cluster_dir = sample_dir.joinpath("cluster")
    if cluster_dir.is_dir():
        dispatch(migrate_cluster_ref_dir,
                 num_cpus=num_cpus,
                 pass_num_cpus=True,
                 ordered=False,
                 raise_on_error=False,
                 as_list=True,
                 args=as_list_of_tuples(cluster_dir.iterdir()))
    fold_dir = sample_dir.joinpath("fold")
    if fold_dir.is_dir():
        dispatch(migrate_fold_ref_dir,
                 num_cpus=num_cpus,
                 pass_num_cpus=True,
                 ordered=False,
                 raise_on_error=False,
                 as_list=True,
                 args=as_list_of_tuples(fold_dir.iterdir()))


def migrate_out_dir(out_dir: Path, num_cpus: int):
    dispatch(migrate_sample_dir,
             num_cpus=num_cpus,
             pass_num_cpus=True,
             ordered=False,
             raise_on_error=False,
             as_list=True,
             args=as_list_of_tuples(out_dir.iterdir()))


@run_func(CMD_MIGRATE)
def run(input_path: Iterable[str | Path], *,
        out_dir: str | Path,
        num_cpus: int):
    """ Migrate output directories from v0.23 to v0.24 """
    input_path = list(input_path)
    if len(input_path) != 1:
        raise ValueError("seismic migrate can process 1 directory at a time, "
                         f"but got {len(input_path)}")
    input_dir = path.sanitize(input_path[0], strict=True)
    out_dir = path.sanitize(out_dir, strict=False)
    if out_dir.exists():
        if os.path.samefile(input_dir, out_dir):
            message = ("For safety, seismic migrate refuses to overwrite "
                       "existing directories, so the output directory must "
                       "be different from the input directory, but got "
                       f"input_dir={input_dir}, out_dir={out_dir}. Please "
                       "specify a different output directory that does not "
                       "yet exist with -o")
        else:
            message = ("For safety, seismic migrate refuses to overwrite "
                       "existing directories, so the output directory must "
                       f"not exist, but got out_dir={out_dir}, which exists. "
                       "Please specify a different output directory that "
                       "does not yet exist with -o")
        raise FileExistsError(message)
    try:
        shutil.copytree(input_dir, out_dir)
        migrate_out_dir(out_dir, num_cpus=num_cpus)
    except Exception as error:
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        raise error


params = [
    arg_input_path,
    opt_out_dir,
    opt_num_cpus
]


@command(CMD_MIGRATE, params=params)
def cli(*args, **kwargs):
    """ Migrate output directories from v0.23 to v0.24 """
    run(*args, **kwargs)
