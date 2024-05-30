import numpy as np
import pandas as pd

from ..core.header import index_order_clusts

rng = np.random.default_rng()


def sim_pclust(order: int, sort: bool = True):
    """ Simulate the proportions of clusters.

    Parameters
    ----------
    order: int
        Number of clusters to simulate; must be ≥ 1.
    sort: bool = False
        Sort the cluster proportions from greatest to least.

    Returns
    -------
    pd.Series
        Simulated proportion of each cluster.
    """
    if order < 1:
        raise ValueError(f"order must be ≥ 1, but got {order}")
    # Simulate cluster proportions with a Dirichlet distribution.
    props = rng.dirichlet(np.ones(order))
    if sort:
        props = np.sort(props)[::-1]
    return pd.Series(props, index=index_order_clusts(order))
