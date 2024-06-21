"""

FASTA Cleaning Module
========================================================================


"""

from logging import getLogger

from click import command

from .cleanfa import clean_fasta
from ..core import path
from ..core.arg import (CMD_CLEANFA,
                        arg_input_path,
                        opt_inplace,
                        opt_out_dir,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.run import run_func
from ..core.task import dispatch

logger = getLogger(__name__)


@run_func(logger.critical)
def run(input_path: tuple[str, ...], *,
        inplace: bool,
        out_dir: str,
        force: bool,
        max_procs: int,
        parallel: bool):
    """ Clean the names and sequences in FASTA files. """
    # List all the files to clean.
    input_files = list(path.find_files_chain(input_path, [path.FastaSeg]))
    # Determine the output path of each cleaned file.
    if inplace:
        # If modifying in-place, then the output and input paths match.
        output_files = input_files
    else:
        # Otherwise, determine the output paths using out_dir.
        out_dir = path.sanitize(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        output_files = path.transpaths(out_dir, *input_files, strict=True)
    # Generate the positional arguments for clean_fasta.
    args = list(zip(input_files, output_files, strict=True))
    # Clean the files; if modifying in-place, force must be True.
    return dispatch(clean_fasta,
                    max_procs,
                    parallel,
                    args=args,
                    kwargs=dict(force=force or inplace),
                    pass_n_procs=False)


params = [
    arg_input_path,
    opt_inplace,
    opt_out_dir,
    opt_force,
    opt_max_procs,
    opt_parallel
]


@command(CMD_CLEANFA, params=params)
def cli(*args, **kwargs):
    """ Clean the names and sequences in FASTA files. """
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
