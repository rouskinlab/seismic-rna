"""

CT Renumbering Module
========================================================================


"""

from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command

from ..core import path
from ..core.arg import (CMD_CTRENUM,
                        docdef,
                        opt_renumber_ct,
                        opt_out_dir,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.parallel import dispatch
from ..core.rna import renumber_ct as renumber_ct_func

logger = getLogger(__name__)

params = [
    opt_renumber_ct,
    opt_out_dir,
    opt_force,
    opt_max_procs,
    opt_parallel
]


@command(CMD_CTRENUM, params=params)
def cli(*args, **kwargs):
    """ Renumber connectivity table (CT) files. """
    return run(*args, **kwargs)


@docdef.auto()
def run(renumber_ct: tuple[tuple[str, int], ...],
        out_dir: str,
        force: bool,
        max_procs: int,
        parallel: bool):
    """ Renumber connectivity table (CT) files. """
    # For each start position, find all CT files to renumber.
    start_cts = {start: list(path.find_files(Path(cts), [path.ConnectTableSeg]))
                 for cts, start in renumber_ct}
    in_cts = list(chain(*start_cts.values()))
    # Make the output directory, if it does not exist.
    out_dir = path.sanitize(out_dir)
    if in_cts:
        out_dir.mkdir(parents=True, exist_ok=True)
    # Determine the output path of each renumbered CT file.
    ct_map = dict()
    for in_ct, out_ct in zip(in_cts,
                             path.transpaths(out_dir, *in_cts, strict=True),
                             strict=True):
        if in_ct in ct_map:
            logger.warning(f"Duplicate CT file to renumber: {in_ct}")
        else:
            ct_map[in_ct] = out_ct
    # Make arguments for renumber_ct.
    args = [(ct, ct_map[ct], start)
            for start, cts in start_cts.items()
            for ct in cts]
    # Renumber the CT files.
    return dispatch(renumber_ct_func,
                    max_procs,
                    parallel,
                    args=args,
                    kwargs=dict(force=force),
                    pass_n_procs=False)

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
