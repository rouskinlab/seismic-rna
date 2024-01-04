from collections import Counter, defaultdict
from datetime import datetime
from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .data import PoolDataset
from .report import PoolReport
from ..core.arg import (CMD_POOL,
                        docdef,
                        arg_input_path,
                        opt_pool,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.data import load_datasets
from ..core.parallel import dispatch
from ..core.write import need_write

logger = getLogger(__name__)

DEFAULT_POOL = "pooled"

params = [
    arg_input_path,
    # Pooling
    opt_pool,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_force,
]


@command(CMD_POOL, params=params)
def cli(*args, pool: str, **kwargs):
    """ Combine samples from the Relate step. """
    if not pool:
        logger.warning(f"{repr(CMD_POOL)} expected a name via --pool, but got "
                       f"{repr(pool)}; defaulting to {repr(DEFAULT_POOL)}")
        pool = DEFAULT_POOL
    return run(*args, pool=pool, **kwargs)


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        pool: str,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        force: bool) -> list[Path]:
    """ Combine samples from the Relate step. """
    if not pool:
        # Exit immediately if no pool name was given.
        return list()
    # Group the datasets by output directory and reference name.
    pools = defaultdict(list)
    for dataset in load_datasets(input_path,
                                 PoolDataset.get_dataset_load_func()):
        pools[dataset.top, dataset.ref].append(dataset.sample)
    # Make each pool of samples.
    return dispatch(make_pool,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=[(out_dir, pool, ref, samples)
                          for (out_dir, ref), samples in pools.items()],
                    kwargs=dict(force=force))


def make_pool(out_dir: Path,
              name: str,
              ref: str,
              samples: Iterable[str], *,
              force: bool):
    """ Combine one or more samples into a pooled sample.

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    name: str
        Name of the pool.
    ref: str
        Name of the reference
    samples: Iterable[str]
        Names of the samples in the pool.
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Pool report file.
    """
    began = datetime.now()
    # Deduplicate and sort the samples.
    sample_counts = Counter(samples)
    if dups := sorted(sample
                      for sample, count in sample_counts.items()
                      if count > 1):
        logger.warning(f"Pool {repr(name)} for output directory {out_dir}, "
                       f"reference {repr(ref)} got duplicate samples: {dups}")
    samples = sorted(sample_counts)
    # Determine the output report file.
    report_file = PoolReport.build_path(top=out_dir, sample=name, ref=ref)
    if need_write(report_file, force):
        logger.info(f"Began pooling samples {samples} into {repr(name)} with "
                    f"reference {repr(ref)} in output directory {out_dir}")
        ended = datetime.now()
        report = PoolReport(sample=name,
                            ref=ref,
                            pooled_samples=samples,
                            began=began,
                            ended=ended)
        report.save(out_dir, force=force)
        logger.info(f"Ended pooling samples {samples} into {repr(name)} with "
                    f"reference {repr(ref)} in output directory {out_dir}")
    return report_file

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
