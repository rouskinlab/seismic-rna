from itertools import chain
from pathlib import Path
from typing import Iterable, Literal, overload

from click import command

from .cluster.data import (
    ClusterCountTabulator,
    ClusterDatasetTabulator,
    load_cluster_dataset,
)
from .core.arg.cmd import CMD_TABLE
from .core.arg.cli import (
    arg_input_path,
    opt_idmut_pos_table,
    opt_idmut_read_table,
    opt_filter_pos_table,
    opt_filter_read_table,
    opt_cluster_pos_table,
    opt_cluster_abundance_table,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
)
from .core.dataset import Dataset, MutsDataset
from .core.run import run_func
from .core.table.write import CountTabulator, DatasetTabulator
from .core.task import dispatch
from .filter.dataset import load_filter_dataset
from .filter.table import FilterCountTabulator, FilterDatasetTabulator
from .idmut.dataset import load_idmut_dataset
from .idmut.table import IDmutCountTabulator, IDmutDatasetTabulator


def tabulate(
    dataset: MutsDataset,
    tabulator_type: type[DatasetTabulator],
    pos_table: bool,
    read_table: bool,
    clust_table: bool,
    force: bool,
    num_cpus: int,
):
    """Write tables for a dataset using the appropriate tabulator.

    Parameters
    ----------
    dataset: MutsDataset
        Dataset to tabulate (from the idmut, filter, or cluster step).
    tabulator_type: type[DatasetTabulator]
        Tabulator class that can process this dataset type.
    pos_table: bool
        If True, write a per-position table.
    read_table: bool
        If True, write a per-read table.
    clust_table: bool
        If True, write a cluster abundance table.
    force: bool
        If True, overwrite existing table files.
    num_cpus: int
        Number of CPUs to use for computation.
    """
    files = tabulator_type(
        dataset=dataset,
        count_pos=(pos_table or clust_table),
        count_read=read_table,
        num_cpus=num_cpus,
    ).write_tables(pos=pos_table, read=read_table, clust=clust_table, force=force)
    return list({file.parent for file in files})


def load_all_datasets(input_path: Iterable[str | Path], **kwargs):
    """Load datasets from all steps of the workflow."""
    if not isinstance(input_path, (tuple, list, set, dict)):
        # Make sure input_path is not an iterator that will be exhausted
        # on the first run-through.
        input_path = list(input_path)
    for load_func in [load_idmut_dataset, load_filter_dataset, load_cluster_dataset]:
        yield from load_func.iterate(input_path, **kwargs)


@overload
def get_tabulator_type(
    dataset_type: type[Dataset], count: Literal[False] = False
) -> type[DatasetTabulator]: ...
@overload
def get_tabulator_type(
    dataset_type: type[Dataset], count: Literal[True]
) -> type[CountTabulator]: ...
def get_tabulator_type(
    dataset_type: type[Dataset], count: bool = False
) -> type[DatasetTabulator] | type[CountTabulator]:
    """Determine which class of Tabulator can process the dataset."""
    if issubclass(dataset_type, load_idmut_dataset.dataset_types):
        return IDmutCountTabulator if count else IDmutDatasetTabulator
    if issubclass(dataset_type, load_filter_dataset.dataset_types):
        return FilterCountTabulator if count else FilterDatasetTabulator
    if issubclass(dataset_type, load_cluster_dataset.dataset_types):
        return ClusterCountTabulator if count else ClusterDatasetTabulator
    raise ValueError(f"No tabulator class exists for {dataset_type.__name__}")


def get_dataset_flags(
    dataset: MutsDataset,
    idmut_pos_table: bool,
    idmut_read_table: bool,
    filter_pos_table: bool,
    filter_read_table: bool,
    cluster_pos_table: bool,
    cluster_abundance_table: bool,
):
    """Return the tabulator and table flags for a dataset."""
    if isinstance(dataset, load_idmut_dataset.dataset_types):
        return idmut_pos_table, idmut_read_table, False
    if isinstance(dataset, load_filter_dataset.dataset_types):
        return filter_pos_table, filter_read_table, False
    if isinstance(dataset, load_cluster_dataset.dataset_types):
        return cluster_pos_table, False, cluster_abundance_table
    raise TypeError(dataset)


@run_func(CMD_TABLE)
def run(
    input_path: Iterable[str | Path],
    *,
    idmut_pos_table: bool,
    idmut_read_table: bool,
    filter_pos_table: bool,
    filter_read_table: bool,
    cluster_pos_table: bool,
    cluster_abundance_table: bool,
    verify_times: bool,
    num_cpus: int,
    force: bool,
) -> list[Path]:
    """Tabulate counts of relationships per read and position."""
    # Load the datasets from the report files.
    args = list()
    for dataset in load_all_datasets(input_path, verify_times=verify_times):
        args.append(
            (
                dataset,
                get_tabulator_type(type(dataset)),
                *get_dataset_flags(
                    dataset,
                    idmut_pos_table,
                    idmut_read_table,
                    filter_pos_table,
                    filter_read_table,
                    cluster_pos_table,
                    cluster_abundance_table,
                ),
            )
        )
    return list(
        chain(
            *dispatch(
                tabulate,
                num_cpus=num_cpus,
                pass_num_cpus=True,
                as_list=False,
                ordered=False,
                raise_on_error=False,
                args=args,
                kwargs=dict(force=force),
            )
        )
    )


params = [
    arg_input_path,
    # Table outputs
    opt_idmut_pos_table,
    opt_idmut_read_table,
    opt_filter_pos_table,
    opt_filter_read_table,
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
    """Tabulate counts of relationships per read and position."""
    return run(*args, **kwargs)
