from pathlib import Path

from .core.table import TableWriter
from .core.dataset import FailedToLoadDatasetError
from .relate.dataset import load_relate_dataset, MutsDataset
from .mask.dataset import load_mask_dataset
from .cluster.data import load_cluster_dataset
from .table import get_tabulator_type


def dataset_from_report(report_path: str | Path,
                        verify_times: bool = True) -> MutsDataset:
    """Load a dataset from a report file.

    Parameters
    ----------
    report_path: str | Path
        The path to a report json file from the relate, mask, or cluster steps.
    verify_times: bool = True
        Ensure that the report file does not have a timestamp
        that is earlier than that of one of its constituents.

    Returns
    -------
    RelateMutsDataset | MaskMutsDataset | ClusterMutsDataset
        The type of MutsDataset returned depends on the report file.
    """
    if isinstance(report_path, str):
        report = Path(report_path)

    dataset = None
    errors = dict()
    for load_func in [load_relate_dataset,
                      load_mask_dataset,
                      load_cluster_dataset]:
        try:
            dataset = load_func(report_path, verify_times=verify_times)
        except Exception as error:
            errors[load_func] = error
            pass
    if dataset is None:
        errmsg = "\n".join(f"{type_name}: {error}"
                           for type_name, error in errors.items())
        raise FailedToLoadDatasetError(
            f"Failed to load {report_path}:\n{errmsg}")
    return dataset


def table_from_dataset(dataset: MutsDataset,
                       table: str = "pos") -> TableWriter:
    """Tabulate a dataset to generate a TableWriter

    Parameters
    ----------
    dataset: RelateMutsDataset | MaskMutsDataset | ClusterMutsDataset
        A dataset from the Relate, Mask, or Cluster steps.
    table: str = "pos"
        The type of table to generate. Valid options include 
        "pos" for per-position table,
        "read" for per-read table, and
        "abundance" for a cluster abundance table.

    Returns
    -------
    TableWriter
        The type of TableWriter returned depends on the Dataset type.
    """
    tabulator_type = get_tabulator_type(type(dataset))
    pos_table = False
    read_table = False
    clust_table = False
    if table == "pos":
        pos_table = True
    elif table == "read":
        read_table = True
    elif table == "abundance":
        if not isinstance(dataset, load_cluster_dataset.dataset_types):
            raise Exception('Abundance tables can only be produced for clustered datasets')
        clust_table = True
    else:
        raise Exception('table kwarg must be "pos", "read" or "abundance"')
    table = next(tabulator_type(dataset=dataset,
                                count_pos=pos_table,
                                count_read=read_table).generate_tables(pos=pos_table,
                                                                       read=read_table,
                                                                       clust=clust_table))
    return table