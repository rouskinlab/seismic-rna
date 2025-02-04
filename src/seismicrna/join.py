from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd
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
                       opt_tmp_pfx,
                       opt_keep_tmp,
                       opt_max_procs,
                       opt_force,
                       extra_defaults)
from .core.dataset import load_datasets
from .core.error import InconsistentValueError
from .core.header import ClustHeader, parse_header
from .core.join import JoinMutsDataset, JoinReport
from .core.logs import logger
from .core.report import JoinedRegionsF
from .core.run import run_func
from .core.task import dispatch
from .core.tmp import release_to_out, with_tmp_dir
from .core.write import need_write
from .mask.dataset import load_mask_dataset
from .mask.report import JoinMaskReport
from .mask.table import MaskDatasetTabulator
from .table import tabulate


def joined_mask_report_exists(top: Path,
                              sample: str,
                              ref: str,
                              joined: str,
                              regs: Iterable[str]):
    """ Return whether a mask report for the joined region exists. """
    mask_report_file = JoinMaskReport.build_path(top=top,
                                                 sample=sample,
                                                 ref=ref,
                                                 reg=joined)
    if not mask_report_file.is_file():
        logger.detail(
            f"Joined mask report {mask_report_file} does not already exist"
        )
        return False
    mask_report = JoinMaskReport.load(mask_report_file)
    report_regs = sorted(mask_report.get_field(JoinedRegionsF))
    cluster_regs = sorted(regs)
    if report_regs != cluster_regs:
        raise InconsistentValueError(
            f"The regions in the mask report file {mask_report_file} "
            f"are {report_regs}, which differ from the regions that "
            f"are being joined for the cluster step: {cluster_regs} "
            "(use --force to suppress this error and overwrite)"
        )
    logger.detail(f"Joined mask report {mask_report_file} already exists")
    return True


def parse_join_clusts_file(file: str | Path):
    """ Parse a file of joined clusters. """
    n_cols = len(ClustHeader.level_names())
    clusts_df = pd.read_csv(file, index_col=list(range(n_cols)))
    header = parse_header(clusts_df.index)
    # Verify the index: use type() rather than isinstance() so that
    # subclasses of ClustHeader will yield False, not True.
    if type(header) is not ClustHeader:
        raise TypeError(f"Expected first {n_cols} of {file} to be a valid "
                        f"cluster header, but got {header}")
    # Rearrange the DataFrame into a dict.
    clusts_dict = {reg: {k: dict() for k in header.ks}
                   for reg in clusts_df.columns}
    for reg, clusts in clusts_df.items():
        for (k, clust), reg_clust in clusts.items():
            if not 1 <= reg_clust <= k:
                raise ValueError(f"Region {repr(reg)} k {k} got a "
                                 f"cluster number out of range: {reg_clust}")
            if reg_clust in clusts_dict[reg][k].values():
                raise ValueError(f"Region {repr(reg)} k={k} got a "
                                 f"repeated cluster number: {reg_clust}")
            clusts_dict[reg][k][clust] = reg_clust
    return clusts_dict


def write_report(report_type: type[JoinReport],
                 out_dir: Path,
                 **kwargs):
    report = report_type(ended=datetime.now(), **kwargs)
    return report.save(out_dir, force=True)


@with_tmp_dir(pass_keep_tmp=False)
def join_regions(out_dir: Path,
                 name: str,
                 sample: str,
                 ref: str,
                 regs: Iterable[str],
                 clustered: bool, *,
                 tmp_dir: Path,
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
    tmp_dir: Path
        Temporary directory.
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
        # To be able to load, the joined dataset must have access to the
        # original mask/cluster dataset(s) in the temporary directory.
        for reg in regs:
            load_function(
                report_type.build_path(top=out_dir,
                                       sample=sample,
                                       ref=ref,
                                       reg=reg),
                verify_times=verify_times
            ).link_data_dirs_to_tmp(tmp_dir)
        # Tabulate the joined dataset.
        dataset = load_function(write_report(report_type,
                                             tmp_dir,
                                             **report_kwargs),
                                verify_times=verify_times)
        tabulate(dataset,
                 tabulator_type,
                 pos_table=pos_table,
                 read_table=read_table,
                 clust_table=clust_table,
                 n_procs=n_procs,
                 force=True)
        # Rewrite the report file with the updated time.
        if clustered:
            # Update joined_clusters in case clusts was initially empty,
            # in which case the dataset would have determined the best
            # way to join the clusters.
            report_kwargs |= dict(joined_clusters=dataset.joined_clusts)
        release_to_out(out_dir,
                       tmp_dir,
                       write_report(report_type,
                                    tmp_dir,
                                    **report_kwargs).parent)
    return report_file.parent


@run_func(CMD_JOIN, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        joined: str,
        join_clusts: str | None,
        mask_pos_table: bool,
        mask_read_table: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        max_procs: int,
        force: bool) -> list[Path]:
    """ Merge regions (horizontally) from the Mask or Cluster step. """
    if not joined:
        raise ValueError(
            "No name for the joined region was given via --joined"
        )
    if join_clusts is not None:
        clusts = parse_join_clusts_file(join_clusts)
    else:
        clusts = dict()
    # Group the datasets by output directory, sample, and reference.
    joins = defaultdict(list)
    load_funcs = {False: load_mask_dataset,
                  True: load_cluster_dataset}
    for clustered, load_func in load_funcs.items():
        for dataset in load_datasets(input_path,
                                     load_func,
                                     verify_times=verify_times):
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
            logger.detail(f"Added regions {regs} for {dataset}")
    # For every joined cluster dataset that does not have a joined mask
    # dataset, also create the joined mask dataset.
    for (top, sample, ref, clustered), regs in list(joins.items()):
        if clustered and (force or not joined_mask_report_exists(top,
                                                                 sample,
                                                                 ref,
                                                                 joined,
                                                                 regs)):
            mask_joins = joins[top, sample, ref, False]
            mask_joins.extend(reg for reg in regs if reg not in mask_joins)
    # Join the masked regions first, then the clustered regions, because
    # the clustered regions require the masked regions.
    results = list()
    for use_clustered in [False, True]:
        args = [(out_dir, joined, sample, ref, regs, clustered)
                for (out_dir, sample, ref, clustered), regs
                in joins.items()
                if clustered == use_clustered]
        kwargs = dict(clusts=clusts,
                      mask_pos_table=mask_pos_table,
                      mask_read_table=mask_read_table,
                      cluster_pos_table=cluster_pos_table,
                      cluster_abundance_table=cluster_abundance_table,
                      verify_times=verify_times,
                      tmp_pfx=tmp_pfx,
                      keep_tmp=keep_tmp,
                      force=force)
        results.extend(dispatch(join_regions,
                                max_procs=max_procs,
                                pass_n_procs=True,
                                args=args,
                                kwargs=kwargs))
    return results


params = [
    arg_input_path,
    opt_joined,
    opt_join_clusts,
    opt_mask_pos_table,
    opt_mask_read_table,
    opt_cluster_pos_table,
    opt_cluster_abundance_table,
    opt_verify_times,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_max_procs,
    opt_force,
]


@command(CMD_JOIN, params=params)
def cli(*args, **kwargs):
    """ Merge regions (horizontally) from the Mask or Cluster step. """
    return run(*args, **kwargs)
