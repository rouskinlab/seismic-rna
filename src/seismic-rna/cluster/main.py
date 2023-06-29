from logging import getLogger
from pathlib import Path

from click import command

from .krun import cluster
from ..core import docdef, path
from ..core.cli import (opt_report, opt_max_clusters, opt_em_runs,
                        opt_min_em_iter, opt_max_em_iter, opt_em_thresh,
                        opt_parallel, opt_max_procs, opt_rerun)
from ..core.parallel import as_list_of_tuples, dispatch

logger = getLogger(__name__)

params = [
    # Input files
    opt_report,
    # Clustering options
    opt_max_clusters,
    opt_em_runs,
    opt_min_em_iter,
    opt_max_em_iter,
    opt_em_thresh,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_rerun,
]


@command(path.MOD_CLUST, params=params)
def cli(*args, **kwargs):
    """ Cluster reads from 'mask' using Expectation-Maximization to find
    alternative structures in the RNA ensemble. """
    return run(*args, **kwargs)


@docdef.auto()
def run(report: tuple[str, ...], *,
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        max_procs: int,
        parallel: bool,
        rerun: bool) -> list[Path]:
    """ Run the clustering module. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    # Run clustering on each set of called mutations.
    files = path.find_files_chain(map(Path, report), [path.MaskRepSeg])
    return dispatch(cluster, max_procs, parallel,
                    args=as_list_of_tuples(files),
                    kwargs=dict(max_order=max_clusters,
                                n_runs=em_runs,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                conv_thresh=em_thresh,
                                rerun=rerun))
