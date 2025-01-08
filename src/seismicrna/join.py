from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

from click import command

from .cluster.dataset import load_cluster_dataset
from .cluster.table import ClusterDatasetTabulator
from .core.arg import (CMD_JOIN,
                       arg_input_path,
                       opt_joined,
                       opt_join_clusts,
                       opt_mask_pos_table,
                       opt_mask_read_table,
                       opt_cluster_pos_table,
                       opt_cluster_abundance_table,
                       opt_verify_times,
                       opt_max_procs,
                       opt_force)
from .core.dataset import load_datasets
from .core.join.cluster import parse_join_clusts_file
from .core.join.dataset import JoinMutsDataset
from .core.join.report import JoinMaskReport, JoinClusterReport
from .core.logs import logger
from .core.run import run_func
from .core.task import dispatch
from .core.write import need_write
from .mask.dataset import load_mask_dataset
from .mask.report import MaskReport
from .mask.table import MaskDatasetTabulator
from .table import tabulate

DEFAULT_JOIN = "joined"


def write_report(report_type: type[JoinMaskReport | JoinClusterReport],
                 out_dir: Path,
                 **kwargs):
    report = report_type(ended=datetime.now(), **kwargs)
    return report.save(out_dir, force=True)


def join_regions(out_dir: Path,
                 name: str,
                 sample: str,
                 ref: str,
                 regs: Iterable[str],
                 clustered: bool, *,
                 clusts: dict[str, dict[int, dict[int, int]]],
                 mask_pos_table: bool,
                 mask_read_table: bool,
                 cluster_pos_table: bool,
                 cluster_abundance_table: bool,
                 verify_times: bool,
                 n_procs: int,
                 force: bool):
    """ Join one or more regions (horizontally).

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    name: str
        Name of the joined region.
    sample: str
        Name of the sample.
    ref: str
        Name of the reference.
    regs: Iterable[str]
        Names of the regions being joined.
    clustered: bool
        Whether the dataset is clustered.
    clusts: dict[str, dict[int, dict[int, int]]]
        For each region, for each number of clusters, the cluster from
        the original region to use as the cluster in the joined region
        (ignored if `clustered` is False).
    mask_pos_table: bool
        Tabulate relationships per position for mask data.
    mask_read_table: bool
        Tabulate relationships per read for mask data
    cluster_pos_table: bool
        Tabulate relationships per position for cluster data.
    cluster_abundance_table: bool
        Tabulate number of reads per cluster for cluster data.
    verify_times: bool
        Verify that report files from later steps have later timestamps.
    n_procs: bool
        Number of processors to use.
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Join report file.
    """
    began = datetime.now()
    # Deduplicate and sort the regions.
    reg_counts = Counter(regs)
    if max(reg_counts.values()) > 1:
        logger.warning(f"Joined region {repr(name)} of sample {repr(sample)}, "
                       f"reference {repr(ref)} in {out_dir} got duplicate "
                       f"regions: {reg_counts}")
    regs = sorted(reg_counts)
    report_kwargs = dict(sample=sample,
                         ref=ref,
                         reg=name,
                         joined_regions=regs,
                         began=began)
    # Determine whether the dataset is clustered.
    if clustered:
        report_kwargs |= dict(joined_clusters=clusts)
        tabulator_type = ClusterDatasetTabulator
        pos_table = cluster_pos_table
        read_table = False
        clust_table = cluster_abundance_table
    else:
        tabulator_type = MaskDatasetTabulator
        pos_table = mask_pos_table
        read_table = mask_read_table
        clust_table = False
    load_function = tabulator_type.load_function()
    _, dataset_type = load_function.dataset_types
    report_type = dataset_type.get_report_type()
    # Determine the output report file.
    report_file = report_type.build_path(top=out_dir,
                                         sample=sample,
                                         ref=ref,
                                         reg=name)
    if need_write(report_file, force):
        # Because a Join report file has the same name as a Mask/Cluster
        # report, it would be possible to overwrite the latter with a
        # Join report, rendering its datasets unusable; prevent this.
        if report_file.is_file():
            # Check whether the report file contains a Join report.
            try:
                report_type.load(report_file)
            except ValueError:
                # The report file does not contain a Join report.
                raise TypeError(
                    f"Overwriting {report_file} with {report_type.__name__} "
                    "would cause data loss"
                )
        # Write the report file with a placeholder time.
        write_report(report_type, out_dir, **report_kwargs)
        # Tabulate the joined dataset.
        dataset = load_function(report_file, verify_times=verify_times)
        tabulate(dataset,
                 tabulator_type,
                 pos_table=pos_table,
                 read_table=read_table,
                 clust_table=clust_table,
                 n_procs=n_procs,
                 force=True)
        # Rewrite the report file with the updated time.
        write_report(report_type, out_dir, **report_kwargs)
    return report_file


@run_func(CMD_JOIN)
def run(input_path: tuple[str, ...], *,
        joined: str,
        join_clusts: str | None,
        # Tabulation
        mask_pos_table: bool,
        mask_read_table: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        # Parallelization
        max_procs: int,
        # Effort
        force: bool) -> list[Path]:
    """ Merge regions (horizontally) from the Mask or Cluster step. """
    if not joined:
        # Exit immediately if no joined name was given.
        return list()
    if join_clusts is not None:
        clusts = parse_join_clusts_file(join_clusts)
    else:
        clusts = dict()
    # Group the datasets by output directory, sample, and reference.
    joins = defaultdict(list)
    load_funcs = {False: load_mask_dataset,
                  True: load_cluster_dataset}
    for clustered, load_func in load_funcs.items():
        for dataset in load_datasets(input_path, load_func):
            # Check whether the dataset was joined.
            if isinstance(dataset, JoinMutsDataset):
                # If so, then use all joined regions.
                regs = dataset.region_names
            else:
                # Otherwise, use just the region of the dataset.
                regs = [dataset.region.name]
            joins[(dataset.top,
                   dataset.sample,
                   dataset.ref,
                   clustered)].extend(regs)
            # If clustered, then check if the corresponding Mask dataset
            # has been joined.
            mask_report_file = MaskReport.build_path(top=dataset.top,
                                                     sample=dataset.sample,
                                                     ref=dataset.ref,
                                                     reg=joined)
            if not mask_report_file.is_file():
                # If not, then also join the Mask dataset.
                mask_joins = joins[(dataset.top,
                                    dataset.sample,
                                    dataset.ref,
                                    False)]
                mask_joins.extend([reg for reg in regs
                                   if reg not in mask_joins])
    # Create each joined region.
    return dispatch(join_regions,
                    max_procs=max_procs,
                    pass_n_procs=True,
                    args=[(out_dir, joined, sample, ref, regs, clustered)
                          for (out_dir, sample, ref, clustered), regs
                          in joins.items()],
                    kwargs=dict(clusts=clusts,
                                mask_pos_table=mask_pos_table,
                                mask_read_table=mask_read_table,
                                cluster_pos_table=cluster_pos_table,
                                cluster_abundance_table=cluster_abundance_table,
                                verify_times=verify_times,
                                force=force))


params = [
    arg_input_path,
    # Joining
    opt_joined,
    opt_join_clusts,
    # Tabulation
    opt_mask_pos_table,
    opt_mask_read_table,
    opt_cluster_pos_table,
    opt_cluster_abundance_table,
    opt_verify_times,
    # Parallelization
    opt_max_procs,
    # Effort
    opt_force,
]


@command(CMD_JOIN, params=params)
def cli(*args, joined: str, **kwargs):
    """ Merge regions (horizontally) from the Mask or Cluster step. """
    if not joined:
        logger.warning(f"{CMD_JOIN} expected a name via --joined, but got "
                       f"{repr(joined)}; defaulting to {repr(DEFAULT_JOIN)}")
        joined = DEFAULT_JOIN
    return run(*args, joined=joined, **kwargs)

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
