from logging import getLogger
from pathlib import Path

from click import command

from .write import cluster
from ..core import path
from ..core.arg import (CMD_CLUST,
                        docdef,
                        arg_input_path,
                        opt_max_clusters,
                        opt_min_nmut_read,
                        opt_em_runs,
                        opt_em_thresh,
                        opt_min_em_iter,
                        opt_max_em_iter,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs,
                        opt_force)
from ..core.parallel import as_list_of_tuples, dispatch

logger = getLogger(__name__)

DEFAULT_ORDER = 2

params = [
    # Input files
    arg_input_path,
    # Clustering options
    opt_max_clusters,
    opt_min_nmut_read,
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
]


@command(CMD_CLUST, params=params)
def cli(*args, max_clusters: int, **kwargs):
    """ Cluster reads from 'mask' using Expectation-Maximization to find
    alternative structures in the RNA ensemble. """
    # When cluster is called via the command "cluster" (instead of via
    # the run() function), assume that clustering is intentional. Thus,
    # override the default max_clusters == 0 (which disables clustering)
    # by setting it to 2 (the minimum non-trivial order of clustering).
    if max_clusters <= 0:
        logger.warning(f"Command '{CMD_CLUST}' got a maximum clustering "
                       f"order of {max_clusters}: setting to {DEFAULT_ORDER}")
        max_clusters = DEFAULT_ORDER
    return run(*args, max_clusters=max_clusters, **kwargs)


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        max_clusters: int,
        em_runs: int,
        min_nmut_read: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        brotli_level: int,
        max_procs: int,
        parallel: bool,
        force: bool) -> list[Path]:
    """ Run the clustering module. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    # Run clustering on each set of called mutations.
    files = path.find_files_chain(map(Path, input_path), [path.MaskRepSeg])
    return dispatch(cluster,
                    max_procs,
                    parallel,
                    args=as_list_of_tuples(files),
                    kwargs=dict(max_order=max_clusters,
                                n_runs=em_runs,
                                min_muts=min_nmut_read,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                conv_thresh=em_thresh,
                                brotli_level=brotli_level,
                                force=force))

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
