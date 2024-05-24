from logging import getLogger
from shutil import rmtree

import numpy as np
import pandas as pd

from ..core import path
from ..core.header import format_clust_name, index_order_clusts
from ..core.rel import MATCH, NOCOV, REL_TYPE
from ..core.rna import RNAProfile, from_ct
from ..core.seq import DNA, Section
from ..core.stats import calc_beta_params, calc_dirichlet_params
from ..fold.rnastructure import fold

logger = getLogger(__name__)

rng = np.random.default_rng()


def sim_paired(refseq: DNA,
               structures: int,
               use_fold: bool = True,
               f_paired: float = 0.5):
    """ Simulate whether each base in paired in one or more structures.

    Parameters
    ----------
    refseq: DNA
        Reference sequence.
    structures: int
        Number of structures to simulate; must be ≥ 1.
    use_fold: bool = True
        Use RNAstructure Fold to predict the structure(s); on failure,
        default to False.
    f_paired: float = 0.5
        If `fold` is False, the fraction of bases to make paired.

    Returns
    -------
    pd.DataFrame
        Whether each base at each position is paired in each structure.
    """
    ref = path.randname(8)
    section = Section(ref, refseq)
    is_paired = pd.DataFrame(index=section.range, dtype=bool)
    if use_fold:
        temp_dir = None
        try:
            sample = path.randname(8)
            data_name = path.randname(8)
            section = Section(ref, refseq)
            data = pd.Series(np.nan, index=is_paired.index)
            profile = RNAProfile(sample=sample,
                                 section=section,
                                 data_sect=section.name,
                                 data_name=data_name,
                                 data=data)
            temp_dir = path.randdir()
            ct_file = fold(profile, out_dir=temp_dir, temp_dir=temp_dir)
            for number, structure in zip(range(structures),
                                         from_ct(ct_file),
                                         strict=False):
                name = format_clust_name(structures, number + 1)
                is_paired[name] = structure.is_paired
        except Exception as error:
            logger.warning(f"Failed to simulate {refseq} with use_fold=True; "
                           f"defaulting to use_fold=False:\n{error}")
            return sim_paired(refseq, structures, False, f_paired)
        finally:
            # Always delete the temporary directory, if it exists.
            if temp_dir is not None:
                rmtree(temp_dir, ignore_errors=True)
    # If use_fold is False or an insufficient number of structures were
    # modeled, then add random structures.
    while (number := is_paired.columns.size) < structures:
        name = format_clust_name(structures, number + 1)
        is_paired[name] = pd.Series(rng.random(is_paired.index.size) < f_paired,
                                    index=is_paired.index)
    return is_paired


def sim_pmut(paired: pd.Series,
             pm: pd.Series,
             um: pd.Series,
             pv: float,
             uv: float):
    """ Simulate mutation rates using two beta distributions for the
    paired and unpaired bases.

    Parameters
    ----------
    paired: pd.Series
        Whether each base is paired.
    pm: pd.Series
        Mean of the mutation rates for paired bases.
    um: pd.Series
        Mean of the mutation rates for unpaired bases.
    pv: float
        Variance of the mutation rates for paired bases, as a fraction
        of its supremum.
    uv: float
        Variance of the mutation rates for unpaired bases, as a fraction
        of its supremum.

    Returns
    -------
    pd.DataFrame
        Mutation rates, with the same index as
    """
    if not isinstance(paired, pd.Series):
        raise TypeError(f"paired must be a Series, "
                        f"but got {type(paired).__name__}")
    if not isinstance(pm, pd.Series):
        raise TypeError(f"pm must be a Series, but got {type(pm).__name__}")
    if not isinstance(um, pd.Series):
        raise TypeError(f"um must be a Series, but got {type(um).__name__}")
    if pm.min() < 0.:
        raise ValueError(f"All pm must be ≥ 0, but got {pm[pm < 0.]}")
    if um.min() < 0.:
        raise ValueError(f"All um must be ≥ 0, but got {um[um < 0.]}")
    if (pm_sum := pm.sum()) > 1.:
        raise ValueError(f"The sum of pm must be ≤ 1, but got {pm_sum}")
    if (um_sum := um.sum()) > 1.:
        raise ValueError(f"The sum of um must be ≤ 1, but got {um_sum}")
    # Count paired/unpaired bases.
    paired = paired.astype(bool, copy=False)
    positions = paired.index
    num_paired = np.count_nonzero(paired)
    num_unpaired = paired.size - num_paired
    # Determine the types of relationships.
    rels = pd.Index.union(pm.index, um.index).astype(REL_TYPE, copy=False)
    if MATCH in rels:
        raise ValueError(f"Matches cannot be a relationship to simulate")
    if NOCOV in rels:
        raise ValueError(f"No coverage cannot be a relationship to simulate")
    rels = rels.append(pd.Index([MATCH]))
    # Copy the mean mutation rates to prevent the originals from being
    # modified, and set their indexes to that of all relationships.
    pm = pm.reindex(rels, fill_value=0.)
    um = um.reindex(rels, fill_value=0.)
    # Simulate matches as the default category.
    pm[MATCH] = 1. - pm_sum
    um[MATCH] = 1. - um_sum
    # Determine which mean mutation rates are not zero.
    pm_nonzero = pm[pm != 0.]
    um_nonzero = um[um != 0.]
    # Simulate the mutation rates for the paired/unpaired bases.
    pmut = pd.DataFrame(0., index=paired.index, columns=rels)
    if pv > 0.:
        pv_nonzero = pv * (pm_nonzero * (1. - pm_nonzero))
        if pm_nonzero.size > 1:
            ppmut = rng.dirichlet(calc_dirichlet_params(pm_nonzero.values,
                                                        pv_nonzero.values),
                                  size=num_paired)
        else:
            ppmut = rng.beta(*calc_beta_params(pm_nonzero.values[0],
                                               pv_nonzero.values[0]),
                             size=num_paired)
    else:
        ppmut = np.broadcast_to(pm_nonzero.values[np.newaxis, :],
                                (num_paired, pm_nonzero.size))
    pmut.loc[positions[paired], pm_nonzero.index] = ppmut
    if uv > 0.:
        uv_nonzero = uv * (um_nonzero * (1. - um_nonzero))
        if um_nonzero.size > 1:
            upmut = rng.dirichlet(calc_dirichlet_params(um_nonzero.values,
                                                        uv_nonzero.values),
                                  size=num_unpaired)
        else:
            upmut = rng.beta(*calc_beta_params(um_nonzero.values[0],
                                               uv_nonzero.values[0]),
                             size=num_unpaired)
    else:
        upmut = np.broadcast_to(um_nonzero.values[np.newaxis, :],
                                (num_unpaired, um_nonzero.size))
    pmut.loc[positions[~paired], um_nonzero.index] = upmut
    # Delete the column for matches.
    return pmut.drop(columns=MATCH)


def sim_props(order: int, sort: bool = True):
    """ Simulate the proportions of `n` clusters.

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
    props = rng.dirichlet(1. - rng.random(order))
    if sort:
        props = np.sort(props)[::-1]
    return pd.Series(props, index=index_order_clusts(order))
