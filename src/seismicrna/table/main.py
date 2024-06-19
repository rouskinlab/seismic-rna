from itertools import chain
from logging import getLogger

from click import command

from .write import write
from ..cluster.data import load_cluster_dataset
from ..core import path
from ..core.arg import (CMD_TABLE,
                        arg_input_path,
                        opt_table_pos,
                        opt_table_read,
                        opt_table_clust,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..mask.data import load_mask_dataset
from ..pool.data import load_relate_dataset

logger = getLogger(__name__)


@run_func(logger.critical)
def run(input_path: tuple[str, ...], *,
        table_pos: bool,
        table_read: bool,
        table_clust: bool,
        force: bool,
        max_procs: int,
        parallel: bool):
    """ Count mutations for each read and position; output tables. """
    relate_reports = path.find_files_chain(
        input_path, load_relate_dataset.report_path_seg_types
    )
    if relate_reports:
        logger.debug(f"Found relate report files: {relate_reports}")
    mask_reports = path.find_files_chain(
        input_path, load_mask_dataset.report_path_seg_types
    )
    if mask_reports:
        logger.debug(f"Found mask report files: {mask_reports}")
    clust_reports = path.find_files_chain(
        input_path, load_cluster_dataset.report_path_seg_types
    )
    if clust_reports:
        logger.debug(f"Found cluster report files: {clust_reports}")
    tasks = as_list_of_tuples(chain(relate_reports,
                                    mask_reports,
                                    clust_reports))
    return list(chain(*dispatch(write,
                                max_procs,
                                parallel,
                                args=tasks,
                                kwargs=dict(table_pos=table_pos,
                                            table_read=table_read,
                                            table_clust=table_clust,
                                            force=force),
                                pass_n_procs=False)))


params = [arg_input_path,
          opt_table_pos,
          opt_table_read,
          opt_table_clust,
          opt_force,
          opt_max_procs,
          opt_parallel]


@command(CMD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Count mutations for each read and position; output tables. """
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
