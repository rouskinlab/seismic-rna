from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

from click import command

from .cluster.data import load_cluster_dataset
from .cluster.report import ClusterReport
from .core.arg import (CMD_JOIN,
                       arg_input_path,
                       opt_joined,
                       opt_join_clusts,
                       opt_max_procs,
                       opt_force)
from .core.data import load_datasets
from .core.join.cluster import parse_join_clusts_file
from .core.join.data import JoinMutsDataset
from .core.join.report import JoinMaskReport, JoinClusterReport
from .core.logs import logger
from .core.run import run_func
from .core.task import dispatch
from .core.write import need_write
from .mask.data import load_mask_dataset
from .mask.report import MaskReport

DEFAULT_JOIN = "joined"


def join_regions(out_dir: Path,
                 name: str,
                 sample: str,
                 ref: str,
                 regs: Iterable[str],
                 clustered: bool, *,
                 clusts: dict[str, dict[int, dict[int, int]]],
                 force: bool):
    """ Join one or more regions.

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
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Pool report file.
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
                         joined_regions=regs)
    # Determine whether the dataset is clustered.
    if clustered:
        report_kwargs |= dict(joined_clusters=clusts)
        join_type = JoinClusterReport
        part_type = ClusterReport
    else:
        join_type = JoinMaskReport
        part_type = MaskReport
    # Determine the output report file.
    report_file = join_type.build_path(top=out_dir,
                                       sample=sample,
                                       ref=ref,
                                       reg=name)
    if need_write(report_file, force):
        # Because Join report files have the same name as Mask/Cluster
        # reports, it would be possible to overwrite the latter with a
        # Join report, rendering their datasets unusable; prevent this.
        if report_file.is_file():
            # Check if the report file contains a Mask/Cluster report.
            try:
                part_type.load(report_file)
            except ValueError:
                # The file does not contain a Mask/Cluster report.
                pass
            else:
                # The file contains a Mask/Cluster report.
                raise TypeError(f"Overwriting {part_type.__name__} in "
                                f"{report_file} with {join_type.__name__} "
                                f"would cause data loss")
            # Check whether the report file contains a Join report.
            try:
                join_type.load(report_file)
            except ValueError:
                # The report file does not contain a Join report.
                raise TypeError(f"Overwriting {report_file} with "
                                f"{join_type.__name__} would cause data loss")
        ended = datetime.now()
        report = join_type(**report_kwargs, began=began, ended=ended)
        report.save(out_dir, force=True)
    return report_file


@run_func(CMD_JOIN)
def run(input_path: tuple[str, ...], *,
        joined: str,
        join_clusts: str | None,
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
                regs = dataset.regs
            else:
                # Otherwise, use just the region of the dataset.
                regs = [dataset.reg]
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
    # Make each joined region.
    return dispatch(join_regions,
                    max_procs=max_procs,
                    pass_n_procs=False,
                    args=[(out_dir, joined, sample, ref, regs, clustered)
                          for (out_dir, sample, ref, clustered), regs
                          in joins.items()],
                    kwargs=dict(clusts=clusts, force=force))


params = [
    arg_input_path,
    # Joining
    opt_joined,
    opt_join_clusts,
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
