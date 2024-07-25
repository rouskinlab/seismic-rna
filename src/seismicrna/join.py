from collections import Counter, defaultdict
from datetime import datetime
from logging import getLogger
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
                       opt_parallel,
                       opt_force)
from .core.data import load_datasets
from .core.run import run_func
from .core.task import dispatch
from .core.write import need_write
from .joinbase.cluster import parse_join_clusts_file
from .joinbase.data import JoinMutsDataset
from .joinbase.report import JoinMaskReport, JoinClusterReport
from .mask.data import load_mask_dataset
from .mask.report import MaskReport

logger = getLogger(__name__)

DEFAULT_JOIN = "joined"


def join_sections(out_dir: Path,
                  name: str,
                  sample: str,
                  ref: str,
                  sects: Iterable[str],
                  clustered: bool, *,
                  clusts: dict[str, dict[int, dict[int, int]]],
                  force: bool):
    """ Join one or more sections.

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    name: str
        Name of the joined section.
    sample: str
        Name of the sample.
    ref: str
        Name of the reference.
    sects: Iterable[str]
        Names of the sections being joined.
    clustered: bool
        Whether the dataset is clustered.
    clusts: dict[str, dict[int, dict[int, int]]]
        For each section, for each number of clusters, the cluster from
        the original section to use as the cluster in the joined section
        (ignored if `clustered` is False).
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Pool report file.
    """
    began = datetime.now()
    # Deduplicate and sort the sections.
    sect_counts = Counter(sects)
    if max(sect_counts.values()) > 1:
        logger.warning(f"Joined section {repr(name)} of sample {repr(sample)}, "
                       f"reference {repr(ref)} in {out_dir} got duplicate "
                       f"sections: {sect_counts}")
    sects = sorted(sect_counts)
    report_kwargs = dict(sample=sample,
                         ref=ref,
                         sect=name,
                         joined_sections=sects)
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
                                       sect=name)
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
        logger.info(f"Began joining sections {sects} into {repr(name)} with "
                    f"sample {repr(sample)}, reference {repr(ref)} in output "
                    f"directory {out_dir}")
        ended = datetime.now()
        report = join_type(**report_kwargs, began=began, ended=ended)
        report.save(out_dir, force=True)
        logger.info(f"Ended joining sections {sects} into {repr(name)} with "
                    f"sample {repr(sample)}, reference {repr(ref)} in output "
                    f"directory {out_dir}")
    return report_file


@run_func(logger.critical)
def run(input_path: tuple[str, ...], *,
        joined: str,
        join_clusts: str | None,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        force: bool) -> list[Path]:
    """ Merge sections (horizontally) from the Mask or Cluster step. """
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
                # If so, then use all joined sections.
                sects = dataset.sects
            else:
                # Otherwise, use just the section of the dataset.
                sects = [dataset.sect]
            joins[(dataset.top,
                   dataset.sample,
                   dataset.ref,
                   clustered)].extend(sects)
            # If clustered, then check if the corresponding Mask dataset
            # has been joined.
            mask_report_file = MaskReport.build_path(top=dataset.top,
                                                     sample=dataset.sample,
                                                     ref=dataset.ref,
                                                     sect=joined)
            if not mask_report_file.is_file():
                # If not, then also join the Mask dataset.
                mask_joins = joins[(dataset.top,
                                    dataset.sample,
                                    dataset.ref,
                                    False)]
                mask_joins.extend([sect for sect in sects
                                   if sect not in mask_joins])
    # Make each joined section.
    return dispatch(join_sections,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=[(out_dir, joined, sample, ref, sects, clustered)
                          for (out_dir, sample, ref, clustered), sects
                          in joins.items()],
                    kwargs=dict(clusts=clusts, force=force))


params = [
    arg_input_path,
    # Joining
    opt_joined,
    opt_join_clusts,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_force,
]


@command(CMD_JOIN, params=params)
def cli(*args, joined: str, **kwargs):
    """ Merge sections (horizontally) from the Mask or Cluster step. """
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
