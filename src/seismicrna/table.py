from itertools import chain
from pathlib import Path

from click import command

from .cluster.data import load_cluster_dataset
from .cluster.table import ClusterDatasetTabulator
from .core.arg import (CMD_TABLE,
                       arg_input_path,
                       opt_relate_pos_table,
                       opt_relate_read_table,
                       opt_mask_pos_table,
                       opt_mask_read_table,
                       opt_cluster_pos_table,
                       opt_cluster_abundance_table,
                       opt_verify_times,
                       opt_max_procs,
                       opt_force)
from .core.data import MutsDataset, load_datasets
from .core.run import run_func
from .core.table import DatasetTabulator
from .core.task import dispatch
from .mask.data import load_mask_dataset
from .mask.table import MaskDatasetTabulator
from .relate.data import load_relate_dataset
from .relate.table import RelateDatasetTabulator


def tabulate(dataset: MutsDataset,
             tabulator: type[DatasetTabulator],
             pos_table: bool,
             read_table: bool,
             clust_table: bool,
             force: bool,
             n_procs: int):
    files = tabulator(dataset=dataset,
                      count_pos=pos_table,
                      count_read=read_table,
                      max_procs=n_procs).write_tables(pos=pos_table,
                                                      read=read_table,
                                                      clust=clust_table,
                                                      force=force)
    return list({file.parent for file in files})


@run_func(CMD_TABLE)
def run(input_path: tuple[str, ...], *,
        relate_pos_table: bool,
        relate_read_table: bool,
        mask_pos_table: bool,
        mask_read_table: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        max_procs: int,
        force: bool) -> list[Path]:
    """ Tabulate counts of relationships per read and position. """
    # Load the datasets from the report files.
    args = list()
    for dataset in load_datasets(input_path,
                                 load_relate_dataset):
        args.append((dataset,
                     RelateDatasetTabulator,
                     relate_pos_table,
                     relate_read_table,
                     False))
    for dataset in load_datasets(input_path,
                                 load_mask_dataset,
                                 verify_times=verify_times):
        args.append((dataset,
                     MaskDatasetTabulator,
                     mask_pos_table,
                     mask_read_table,
                     False))
    for dataset in load_datasets(input_path,
                                 load_cluster_dataset,
                                 verify_times=verify_times):
        args.append((dataset,
                     ClusterDatasetTabulator,
                     cluster_pos_table,
                     False,
                     cluster_abundance_table))
    return list(chain(*dispatch(tabulate,
                                max_procs=max_procs,
                                pass_n_procs=True,
                                args=args,
                                kwargs=dict(force=force))))


params = [
    arg_input_path,
    # Table outputs
    opt_relate_pos_table,
    opt_relate_read_table,
    opt_mask_pos_table,
    opt_mask_read_table,
    opt_cluster_pos_table,
    opt_cluster_abundance_table,
    # Validation
    opt_verify_times,
    # Parallelization
    opt_max_procs,
    # Effort
    opt_force,
]


@command(CMD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Tabulate counts of relationships per read and position. """
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
