from collections import defaultdict
from logging import getLogger
from pathlib import Path

from click import command

from .meta import parse_refs_metadata, parse_samples_metadata
from .web import export_sample
from ..core.arg import (CMD_EXPORT,
                        arg_input_path,
                        opt_samples_meta,
                        opt_refs_meta,
                        opt_all_pos,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.run import run_func
from ..core.task import dispatch
from ..table.base import (MaskTable,
                          ClustTable,
                          ClustFreqTable)
from ..table.load import load_all_tables

logger = getLogger(__name__)


@run_func(logger.critical)
def run(input_path: tuple[str, ...], *,
        samples_meta: str,
        refs_meta: str,
        all_pos: bool,
        force: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Export a file of each sample for the seismic-graph web app. """
    tables = defaultdict(list)
    samples_metadata = (parse_samples_metadata(Path(samples_meta))
                        if samples_meta
                        else dict())
    refs_metadata = (parse_refs_metadata(Path(refs_meta))
                     if refs_meta
                     else dict())
    for table in load_all_tables(input_path):
        if isinstance(table, (MaskTable, ClustTable, ClustFreqTable)):
            tables[(table.top, table.sample)].append(table)
    return list(dispatch(export_sample,
                         max_procs,
                         parallel,
                         pass_n_procs=False,
                         args=list(tables.items()),
                         kwargs=dict(samples_metadata=samples_metadata,
                                     refs_metadata=refs_metadata,
                                     all_pos=all_pos,
                                     force=force)))


params = [
    arg_input_path,
    opt_samples_meta,
    opt_refs_meta,
    opt_all_pos,
    opt_force,
    opt_max_procs,
    opt_parallel,
]


@command(CMD_EXPORT, params=params)
def cli(*args, **kwargs):
    """ Export each sample to SEISMICgraph (https://seismicrna.org). """
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
