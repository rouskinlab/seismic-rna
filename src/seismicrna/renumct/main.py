from collections import defaultdict
from pathlib import Path
from typing import Iterable

from click import command

from ..core import path
from ..core.arg import (CMD_RENUMCT,
                        opt_ct_pos_5,
                        opt_inplace,
                        opt_out_dir,
                        opt_force,
                        opt_num_cpus)
from ..core.logs import logger
from ..core.rna import renumber_ct as renumber_ct
from ..core.run import run_func
from ..core.task import dispatch


@run_func(CMD_RENUMCT)
def run(*,
        ct_pos_5: Iterable[tuple[str, int]],
        inplace: bool,
        out_dir: str | Path,
        force: bool,
        num_cpus: int):
    """ Renumber connectivity table (CT) files given a 5' position. """
    # For each start position, find all files to renumber.
    start_files = {start: list(path.find_files(Path(files),
                                               [path.ConnectTableSeg]))
                   for files, start in ct_pos_5}
    # Check for files with multiple start positions.
    file_starts = defaultdict(set)
    for start, files in start_files.items():
        for file in files:
            file_starts[file].add(start)
    multi_starts = {file: sorted(starts)
                    for file, starts in file_starts.items()
                    if len(starts) > 1}
    if multi_starts:
        logger.warning(f"Got multiple start positions for {multi_starts}; "
                       f"using the largest start position for each file")
    # Use the largest start position for each file.
    file_start = {file: max(starts) for file, starts in file_starts.items()}
    # Determine the output path of each renumbered file.
    if inplace:
        # If modifying in-place, then the output and input paths match.
        out_path = file_start
    else:
        # Otherwise, determine the output paths using out_dir.
        out_dir = path.sanitize(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = path.transpaths(out_dir, *file_start, strict=True)
    file_out = dict(zip(file_start, out_path, strict=True))
    # Generate the positional arguments for renumber_ct.
    args = [(file, file_out[file], file_start[file]) for file in file_start]
    # Renumber the files; if modifying in-place, force must be True.
    return dispatch(renumber_ct,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=args,
                    kwargs=dict(force=force or inplace))


params = [
    opt_ct_pos_5,
    opt_inplace,
    opt_out_dir,
    opt_force,
    opt_num_cpus,
]


@command(CMD_RENUMCT, params=params)
def cli(*args, **kwargs):
    """ Renumber connectivity table (CT) files given a 5' position. """
    return run(*args, **kwargs)
