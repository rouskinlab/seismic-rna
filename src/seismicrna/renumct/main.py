"""

CT Renumbering Module
========================================================================


"""

from collections import defaultdict
from logging import getLogger
from pathlib import Path

from click import command

from ..core import path
from ..core.arg import (CMD_RENUMCT,
                        opt_ct_pos_5,
                        opt_inplace,
                        opt_out_dir,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.rna import renumber_ct as renumber_ct
from ..core.run import run_func
from ..core.task import dispatch

logger = getLogger(__name__)


@run_func(logger.critical)
def run(*,
        ct_pos_5: tuple[tuple[str, int], ...],
        inplace: bool,
        out_dir: str,
        force: bool,
        max_procs: int,
        parallel: bool):
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
                    max_procs,
                    parallel,
                    args=args,
                    kwargs=dict(force=force or inplace),
                    pass_n_procs=False)


params = [
    opt_ct_pos_5,
    opt_inplace,
    opt_out_dir,
    opt_force,
    opt_max_procs,
    opt_parallel
]


@command(CMD_RENUMCT, params=params)
def cli(*args, **kwargs):
    """ Renumber connectivity table (CT) files given a 5' position. """
    return run(*args, **kwargs)

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
