from logging import getLogger
from pathlib import Path

from click import command

from .write import cluster
from ..core import path
from ..core.arg import (CMD_CLUSTER,
                        arg_input_path,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_max_clusters,
                        opt_em_runs,
                        opt_em_thresh,
                        opt_min_em_iter,
                        opt_max_em_iter,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs,
                        opt_force)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..mask.data import load_mask_dataset

logger = getLogger(__name__)

DEFAULT_ORDER = 2


@run_func(logger.critical, with_tmp=True)
def run(input_path: tuple[str, ...], *,
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        brotli_level: int,
        max_procs: int,
        parallel: bool,
        force: bool,
        tmp_dir: Path) -> list[Path]:
    """ Infer alternative structures by clustering reads' mutations. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    # Find the mask report files.
    report_files = path.find_files_chain(
        input_path, load_mask_dataset.report_path_seg_types
    )
    # Cluster each mask dataset.
    return dispatch(cluster,
                    max_procs,
                    parallel,
                    pass_n_procs=True,
                    args=as_list_of_tuples(report_files),
                    kwargs=dict(max_order=max_clusters,
                                n_runs=em_runs,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                em_thresh=em_thresh,
                                brotli_level=brotli_level,
                                force=force,
                                tmp_dir=tmp_dir))


params = [
    # Input files
    arg_input_path,
    # Clustering options
    opt_max_clusters,
    opt_em_runs,
    opt_em_thresh,
    opt_min_em_iter,
    opt_max_em_iter,
    # Compression
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_force,
    opt_tmp_pfx,
    opt_keep_tmp,
]


@command(CMD_CLUSTER, params=params)
def cli(*args, max_clusters: int, **kwargs):
    """ Infer alternative structures by clustering reads' mutations. """
    # When cluster is called via the command "cluster" (instead of via
    # the run() function), assume that clustering is intentional. Thus,
    # override the default max_clusters == 0 (which disables clustering)
    # by setting it to 2 (the minimum non-trivial order of clustering).
    if max_clusters <= 0:
        logger.warning(f"{repr(CMD_CLUSTER)} expected --max-clusters to be ≥ 1, "
                       f"but got {max_clusters}; defaulting to {DEFAULT_ORDER}")
        max_clusters = DEFAULT_ORDER
    return run(*args, max_clusters=max_clusters, **kwargs)

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
