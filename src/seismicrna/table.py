from itertools import chain
from pathlib import Path
from typing import Iterable

from click import command

from .cluster.data import (ClusterCountTabulator,
                           ClusterDatasetTabulator,
                           load_cluster_dataset)
from .core.arg import (CMD_TABLE,
                       arg_input_path,
                       opt_relate_pos_table,
                       opt_relate_read_table,
                       opt_mask_pos_table,
                       opt_mask_read_table,
                       opt_cluster_pos_table,
                       opt_cluster_abundance_table,
                       opt_verify_times,
                       opt_num_cpus,
                       opt_force)
from .core.dataset import Dataset, MutsDataset
from .core.run import run_func
from .core.table import DatasetTabulator
from .core.task import dispatch
from .mask.dataset import load_mask_dataset
from .mask.table import MaskCountTabulator, MaskDatasetTabulator
from .relate.dataset import load_relate_dataset
from .relate.table import RelateCountTabulator, RelateDatasetTabulator


def tabulate(dataset: MutsDataset,
             tabulator_type: type[DatasetTabulator],
             pos_table: bool,
             read_table: bool,
             clust_table: bool,
             force: bool,
             num_cpus: int):
    files = tabulator_type(dataset=dataset,
                           count_pos=(pos_table or clust_table),
                           count_read=read_table,
                           num_cpus=num_cpus).write_tables(pos=pos_table,
                                                           read=read_table,
                                                           clust=clust_table,
                                                           force=force)
    return list({file.parent for file in files})


def load_all_datasets(input_path: Iterable[str | Path], **kwargs):
    """ Load datasets from all steps of the workflow. """
    if not isinstance(input_path, (tuple, list, set, dict)):
        # Make sure input_path is not an iterator that will be exhausted
        # on the first run-through.
        input_path = list(input_path)
    for load_func in [load_relate_dataset,
                      load_mask_dataset,
                      load_cluster_dataset]:
        yield from load_func.iterate(input_path, **kwargs)


def get_tabulator_type(dataset_type: type[Dataset], count: bool = False):
    """ Determine which class of Tabulator can process the dataset. """
    if issubclass(dataset_type, load_relate_dataset.dataset_types):
        return RelateCountTabulator if count else RelateDatasetTabulator
    if issubclass(dataset_type, load_mask_dataset.dataset_types):
        return MaskCountTabulator if count else MaskDatasetTabulator
    if issubclass(dataset_type, load_cluster_dataset.dataset_types):
        return ClusterCountTabulator if count else ClusterDatasetTabulator
    raise ValueError(f"No tabulator class exists for {dataset_type.__name__}")


def get_dataset_flags(dataset: MutsDataset,
                      relate_pos_table: bool,
                      relate_read_table: bool,
                      mask_pos_table: bool,
                      mask_read_table: bool,
                      cluster_pos_table: bool,
                      cluster_abundance_table: bool):
    """ Return the tabulator and table flags for a dataset. """
    if isinstance(dataset, load_relate_dataset.dataset_types):
        return relate_pos_table, relate_read_table, False
    if isinstance(dataset, load_mask_dataset.dataset_types):
        return mask_pos_table, mask_read_table, False
    if isinstance(dataset, load_cluster_dataset.dataset_types):
        return cluster_pos_table, False, cluster_abundance_table
    raise TypeError(dataset)


@run_func(CMD_TABLE)
def run(input_path: Iterable[str | Path], *,
        relate_pos_table: bool,
        relate_read_table: bool,
        mask_pos_table: bool,
        mask_read_table: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        num_cpus: int,
        force: bool) -> list[Path]:
    """ Tabulate counts of relationships per read and position. """
    # Load the datasets from the report files.
    args = list()
    for dataset in load_all_datasets(input_path, verify_times=verify_times):
        args.append((dataset,
                     get_tabulator_type(type(dataset)),
                     *get_dataset_flags(dataset,
                                        relate_pos_table,
                                        relate_read_table,
                                        mask_pos_table,
                                        mask_read_table,
                                        cluster_pos_table,
                                        cluster_abundance_table)))
    return list(chain(*dispatch(tabulate,
                                num_cpus=num_cpus,
                                pass_num_cpus=True,
                                as_list=False,
                                ordered=False,
                                raise_on_error=False,
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
    opt_num_cpus,
    # Effort
    opt_force,
]


@command(CMD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Tabulate counts of relationships per read and position. """
    return run(*args, **kwargs)
