import gzip
import os
import shutil
from pathlib import Path
from typing import Iterable

from click import command

from .core import path
from .core.arg import (CMD_MIGRATE,
                       arg_input_path,
                       opt_out_dir,
                       opt_num_cpus)
from .core.run import run_func
from .core.task import as_list_of_tuples, dispatch


class FindAndReplaceError(ValueError):
    pass


def find_and_replace(file: Path, find: str, replace: str, count: int = 1):
    """ Replace occurrences of a string in a plain-text or gzipped file.

    Parameters
    ----------
    file: Path
        Path to the file to modify in place.
    find: str
        Substring to search for.
    replace: str
        Replacement string.
    count: int
        Expected number of occurrences of `find`; if the actual count
        differs a FindAndReplaceError is raised.  Pass -1 to skip the
        check.
    """
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


def _rewrite_report_content(report_file: Path, old_step: str, new_step: str):
    """ Update JSON fields inside a report file that reference the old
    step name: the checksum dict key (e.g., ``"relate":`` -> ``"idmut":``)
    and, for filter reports, the n_mask_iter -> n_filter_iter rename. """
    try:
        find_and_replace(report_file,
                         f'"{old_step}":',
                         f'"{new_step}":',
                         count=-1)
    except FindAndReplaceError:
        pass
    if old_step == "mask":
        try:
            find_and_replace(report_file,
                             '"n_mask_iter":',
                             '"n_filter_iter":',
                             count=-1)
        except FindAndReplaceError:
            pass


def _migrate_step_layout(sample_dir: Path, old_step: str, new_step: str):
    """ Rename ``sample_dir/<old_step>/`` -> ``sample_dir/<new_step>/``,
    including recursively renaming any files whose names start with
    ``<old_step>-`` to ``<new_step>-`` (e.g., for relate -> idmut, an
    ``<old_step>-report.json`` becomes ``<new_step>-report.json`` and
    each ``<old_step>-batch-N.brickle`` becomes ``<new_step>-batch-N.brickle``),
    and updating step-name references inside the renamed report JSON
    files. Idempotent: skips if the new directory already exists. """
    old_step_dir = sample_dir.joinpath(old_step)
    new_step_dir = sample_dir.joinpath(new_step)
    if new_step_dir.exists():
        return
    if not old_step_dir.is_dir():
        return
    # First rename files inside, then the directory itself.
    for entry in list(old_step_dir.rglob("*")):
        if entry.is_file() and entry.name.startswith(f"{old_step}-"):
            new_name = entry.name.replace(f"{old_step}-",
                                          f"{new_step}-",
                                          1)
            entry.rename(entry.with_name(new_name))
    old_step_dir.rename(new_step_dir)
    # Update step-name references in renamed report JSON files.
    for report_file in new_step_dir.rglob(f"{new_step}-report.json"):
        _rewrite_report_content(report_file, old_step, new_step)


def migrate_sample_dir(sample_dir: Path, num_cpus: int):
    _migrate_step_layout(sample_dir, "relate", "idmut")
    _migrate_step_layout(sample_dir, "mask", "filter")


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
    """ Migrate output directories from v0.24 to v0.25 """
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
    """ Migrate output directories from v0.24 to v0.25 """
    run(*args, **kwargs)
