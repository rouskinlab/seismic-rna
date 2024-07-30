import os
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (opt_ct_file,
                        opt_clust_conc,
                        opt_force,
                        opt_parallel,
                        opt_max_procs)
from ..core.header import ClustHeader
from ..core.rna import from_ct
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

rng = np.random.default_rng()

PROPORTION = "Proportion"


def sim_pclust(num_clusters: int,
               concentration: float | None = None,
               sort: bool = True):
    """ Simulate proportions of clusters using a Dirichlet distribution.

    Parameters
    ----------
    num_clusters: int
        Number of clusters to simulate; must be ≥ 1.
    concentration: float | None = None
        Concentration parameter for Dirichlet distribution; defaults to
        1 / (`num_clusters` - 1); must be > 0.
    sort: bool = False
        Sort the cluster proportions from greatest to least.

    Returns
    -------
    pd.Series
        Simulated proportion of each cluster.
    """
    if num_clusters < 1:
        raise ValueError(f"num_clusters must be ≥ 1, but got {num_clusters}")
    if num_clusters == 1:
        props = np.ones(num_clusters)
    else:
        if concentration is None:
            concentration = 1. / (num_clusters - 1.)
        props = rng.dirichlet(np.full(num_clusters, concentration))
        if sort:
            props = np.sort(props)[::-1]
    return pd.Series(props,
                     index=ClustHeader(ks=[num_clusters]).index,
                     name=PROPORTION)


def sim_pclust_ct(ct_file: Path, *,
                  concentration: float,
                  force: bool):
    pclust_file = ct_file.with_suffix(path.PARAM_CLUSTS_EXT)
    if need_write(pclust_file, force):
        num_structures = sum(1 for _ in from_ct(ct_file))
        pclust = sim_pclust(num_structures, concentration)
        pclust.to_csv(pclust_file)
    return pclust_file


def load_pclust(pclust_file: Path):
    """ Load cluster proportions from a file. """
    return pd.read_csv(
        pclust_file,
        index_col=list(range(ClustHeader.num_levels()))
    )[PROPORTION]


@run_func(logger.critical)
def run(*,
        ct_file: tuple[str, ...],
        clust_conc: float,
        force: bool,
        parallel: bool,
        max_procs: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(sim_pclust_ct,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(concentration=(clust_conc if clust_conc
                                               else None),
                                force=force))


params = [
    opt_ct_file,
    opt_clust_conc,
    opt_force,
    opt_parallel,
    opt_max_procs
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the proportions of 5' and 3' end coordinates. """
    run(*args, **kwargs)
