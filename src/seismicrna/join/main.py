from collections import Counter, defaultdict
from datetime import datetime
from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .data import JoinMaskMutsDataset
from .report import JoinMaskReport
from ..core.arg import (CMD_JOIN,
                        docdef,
                        arg_input_path,
                        opt_joined,
                        opt_join_clusts,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.data import LoadFunction, load_datasets
from ..core.parallel import dispatch
from ..core.write import need_write
from ..mask.data import MaskMutsDataset
from ..mask.report import MaskReport

logger = getLogger(__name__)

DEFAULT_JOIN = "joined"


def join_sections(out_dir: Path,
                  name: str,
                  sample: str,
                  ref: str,
                  sects: Iterable[str], *,
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
    # Determine the output report file.
    report_file = JoinMaskReport.build_path(top=out_dir,
                                            sample=sample,
                                            ref=ref,
                                            sect=name)
    if need_write(report_file, force):
        # Because Mask and Join report files have the same name, it
        # would be possible to overwrite a Mask report with a Join
        # report, rendering the Mask dataset unusable; prevent this.
        if report_file.is_file():
            # Check whether the report file contains a Mask report.
            try:
                MaskReport.load(report_file)
            except ValueError:
                # The report file does not contain a Mask report.
                pass
            else:
                # The report file contains a Mask report.
                raise TypeError(f"Overwriting {MaskReport.__name__} in "
                                f"{report_file} with {JoinMaskReport.__name__} "
                                f"would cause data loss")
            # Check whether the report file contains a Join report.
            try:
                JoinMaskReport.load(report_file)
            except ValueError:
                # The report file does not contain a Pool report.
                raise TypeError(f"Overwriting {report_file} with "
                                f"{JoinMaskReport.__name__} would cause data loss")
        logger.info(f"Began joining sections {sects} into {repr(name)} with "
                    f"sample {repr(sample)}, reference {repr(ref)} in output "
                    f"directory {out_dir}")
        ended = datetime.now()
        report = JoinMaskReport(sample=sample,
                                ref=ref,
                                sect=name,
                                joined_sections=sects,
                                began=began,
                                ended=ended)
        report.save(out_dir, force=force)
        logger.info(f"Ended joining sections {sects} into {repr(name)} with "
                    f"sample {repr(sample)}, reference {repr(ref)} in output "
                    f"directory {out_dir}")
    return report_file


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        joined: str,
        join_clusts: dict,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        force: bool) -> list[Path]:
    """ Join sections (horizontally) from the Mask or Cluster step. """
    if not joined:
        # Exit immediately if no joined name was given.
        return list()
    # Group the datasets by output directory, sample, and reference.
    joins = defaultdict(list)
    for dataset in load_datasets(input_path,
                                 LoadFunction(MaskMutsDataset,
                                              JoinMaskMutsDataset)):
        # Check whether the dataset was joined.
        if isinstance(dataset, JoinMaskMutsDataset):
            # If so, then use all joined sections.
            sects = dataset.sects
        else:
            # Otherwise, use just the section of the dataset.
            sects = [dataset.sect]
        joins[dataset.top, dataset.sample, dataset.ref].extend(sects)
    # Make each joined section.
    return dispatch(join_sections,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=[(out_dir, joined, sample, ref, sects)
                          for (out_dir, sample, ref), sects
                          in joins.items()],
                    kwargs=dict(force=force))


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
    """ Join sections (horizontally) from the Mask or Cluster step. """
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
