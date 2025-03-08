import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (opt_ct_file,
                        opt_clust_conc,
                        opt_force,
                        opt_num_cpus)
from ..core.header import ClustHeader
from ..core.rna import from_ct
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..core.validate import require_atleast
from ..core.write import need_write

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
    concentration: float | None
        Concentration parameter for Dirichlet distribution; defaults to
        1 / (`num_clusters` - 1); must be > 0.
    sort: bool
        Sort the cluster proportions from greatest to least.

    Returns
    -------
    pd.Series
        Simulated proportion of each cluster.
    """
    require_atleast("num_clusters", num_clusters, 1)
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
        index_col=list(range(ClustHeader.get_num_levels()))
    )[PROPORTION]


@run_func(COMMAND)
def run(*,
        ct_file: Iterable[str | Path],
        clust_conc: float,
        force: bool,
        num_cpus: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(sim_pclust_ct,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(concentration=(clust_conc if clust_conc
                                               else None),
                                force=force))


params = [
    opt_ct_file,
    opt_clust_conc,
    opt_force,
    opt_num_cpus
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the proportions of 5' and 3' end coordinates. """
    run(*args, **kwargs)
