from logging import getLogger
from pathlib import Path

from click import command

from .write import cluster
from ..core import path
from ..core.arg import (CMD_CLUSTER,
                        arg_input_path,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_min_clusters,
                        opt_max_clusters,
                        opt_em_runs,
                        opt_em_thresh,
                        opt_min_em_iter,
                        opt_max_em_iter,
                        opt_max_pearson_run,
                        opt_min_nrmsd_run,
                        opt_max_loglike_vs_best,
                        opt_min_pearson_vs_best,
                        opt_max_nrmsd_vs_best,
                        opt_try_all_ks,
                        opt_write_all_ks,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs,
                        opt_force)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..mask.data import load_mask_dataset

logger = getLogger(__name__)

DEFAULT_MIN_CLUSTERS = 1


@run_func(logger.critical, with_tmp=True)
def run(input_path: tuple[str, ...], *,
        min_clusters: int,
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        min_nrmsd_run: float,
        max_pearson_run: float,
        max_loglike_vs_best: float,
        min_pearson_vs_best: float,
        max_nrmsd_vs_best: float,
        try_all_ks: bool,
        write_all_ks: bool,
        brotli_level: int,
        max_procs: int,
        parallel: bool,
        force: bool,
        tmp_dir: Path) -> list[Path]:
    """ Infer alternative structures by clustering reads' mutations. """
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
                    kwargs=dict(min_clusters=min_clusters,
                                max_clusters=max_clusters,
                                em_runs=em_runs,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                em_thresh=em_thresh,
                                min_nrmsd_run=min_nrmsd_run,
                                max_pearson_run=max_pearson_run,
                                max_loglike_vs_best=max_loglike_vs_best,
                                min_pearson_vs_best=min_pearson_vs_best,
                                max_nrmsd_vs_best=max_nrmsd_vs_best,
                                try_all_ks=try_all_ks,
                                write_all_ks=write_all_ks,
                                brotli_level=brotli_level,
                                force=force,
                                tmp_dir=tmp_dir))


params = [
    # Input files
    arg_input_path,
    # Clustering options
    opt_min_clusters,
    opt_max_clusters,
    opt_em_runs,
    opt_em_thresh,
    opt_min_em_iter,
    opt_max_em_iter,
    opt_max_pearson_run,
    opt_min_nrmsd_run,
    opt_max_loglike_vs_best,
    opt_min_pearson_vs_best,
    opt_max_nrmsd_vs_best,
    opt_try_all_ks,
    opt_write_all_ks,
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
def cli(*args, **kwargs):
    """ Infer alternative structures by clustering reads' mutations. """
    return run(*args, **kwargs)

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
