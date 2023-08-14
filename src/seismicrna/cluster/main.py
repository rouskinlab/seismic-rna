from logging import getLogger
from pathlib import Path

from click import command

from .krun import cluster
from ..core import docdef, path
from ..core.cli import (arg_input_file, opt_max_clusters,
                        opt_min_nmut_read, opt_em_runs, opt_em_thresh,
                        opt_min_em_iter, opt_max_em_iter,
                        opt_parallel, opt_max_procs, opt_rerun)
from ..core.parallel import as_list_of_tuples, dispatch

logger = getLogger(__name__)

DEFAULT_ORDER = 2

params = [
    # Input files
    arg_input_file,
    # Clustering options
    opt_max_clusters,
    opt_min_nmut_read,
    opt_em_runs,
    opt_em_thresh,
    opt_min_em_iter,
    opt_max_em_iter,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_rerun,
]


@command(path.MOD_CLUST, params=params)
def cli(*args, max_clusters: int, **kwargs):
    """ Cluster reads from 'mask' using Expectation-Maximization to find
    alternative structures in the RNA ensemble. """
    # When cluster is called via the command "cluster" (instead of via
    # the run() function), assume that clustering is intentional. Thus,
    # override the default max_clusters == 0 (which disables clustering)
    # by setting it to 2 (the minimum non-trivial order of clustering).
    if max_clusters <= 0:
        logger.warning(f"Command '{path.MOD_CLUST}' got a maximum clustering "
                       f"order of {max_clusters}: setting to {DEFAULT_ORDER}")
        max_clusters = DEFAULT_ORDER
    return run(*args, max_clusters=max_clusters, **kwargs)


@docdef.auto()
def run(input_file: tuple[str, ...], *,
        max_clusters: int,
        min_nmut_read: int,
        em_runs: int,
        em_thresh: float,
        min_em_iter: int,
        max_em_iter: int,
        max_procs: int,
        parallel: bool,
        rerun: bool) -> list[Path]:
    """ Run the clustering module. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    # Run clustering on each set of called mutations.
    files = path.find_files_chain(map(Path, input_file), [path.MaskRepSeg])
    return dispatch(cluster, max_procs, parallel,
                    args=as_list_of_tuples(files),
                    kwargs=dict(max_order=max_clusters,
                                min_muts=min_nmut_read,
                                n_runs=em_runs,
                                conv_thresh=em_thresh,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                rerun=rerun))
