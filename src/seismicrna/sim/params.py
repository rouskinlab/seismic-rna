from logging import getLogger
from shutil import rmtree
from typing import Any

import numpy as np
import pandas as pd

from ..core import path
from ..core.header import format_clust_name, index_order_clusts
from ..core.rel import (MATCH,
                        NOCOV,
                        DELET,
                        INDEL,
                        SUB_A,
                        SUB_C,
                        SUB_G,
                        SUB_T,
                        ANY_B,
                        ANY_D,
                        ANY_H,
                        ANY_V,
                        ANY_N,
                        REL_TYPE)
from ..core.rna import RNAProfile, from_ct
from ..core.seq import (BASE_NAME,
                        BASEA,
                        BASEC,
                        BASEG,
                        BASET,
                        BASEN,
                        DNA,
                        Section)
from ..core.stats import calc_beta_params, calc_dirichlet_params
from ..fold.rnastructure import fold

logger = getLogger(__name__)

rng = np.random.default_rng()


def verify_proportions(p: Any):
    """ Verify that `p` is a valid set of proportions:

    - Every element of `p` must be ≥ 0 and ≤ 1.
    - The sum of `p` must equal 1.

    Parameters
    ----------
    p: Any
        Proportions to verify; must be a NumPy array or convertable into
        a NumPy array.
    """
    arr = np.asarray_chkfinite(p)
    if not np.isclose(sum_p := arr.sum(), 1.):
        raise ValueError(f"Proportions must sum to 1, but got {sum_p}")
    if (min_p := arr.min()) < 0.:
        raise ValueError(f"Every proportion must be ≥ 0, "
                         f"but minimum is {min_p}")
    if (max_p := arr.max()) > 1.:
        raise ValueError(f"Every proportion must be ≤ 1, "
                         f"but maximum is {max_p}")


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


def make_pmut_means(*,
                    ploq: float = 0.05,
                    pam: float,
                    pac: float = 0.30,
                    pag: float = 0.16,
                    pat: float = 0.50,
                    pcm: float,
                    pca: float = 0.32,
                    pcg: float = 0.32,
                    pct: float = 0.32,
                    pgm: float,
                    pga: float = 0.32,
                    pgc: float = 0.32,
                    pgt: float = 0.32,
                    ptm: float,
                    pta: float = 0.32,
                    ptc: float = 0.32,
                    ptg: float = 0.32,
                    pnm: float = 0.00,
                    pnd: float = 0.04):
    """ Generate mean mutation rates.

    Mutations are assumed to behave as follows:

    -   A base `n` mutates with probability `pnm`.
        -   If it mutates, then it is a substitution with probability
            (`pna` + `pnc` + `png` + `pnt`).
            -   If it is a substitution, then it is high-quality with
                probability (1 - `ploq`).
            -   Otherwise, it is low-quality.
        -   Otherwise, it is a deletion.
    -   Otherwise, it is low-quality with probability `ploq`.

    So the overall probability of being low-quailty is the probability
    given a mutation, `pnm` * (`pna` + `pnc` + `png` + `pnt`) * `ploq`,
    plus the probability given no mutation, (1 - `pnm`) * `ploq`, which
    equals `ploq` * (1 - `pam` * (1 - (`pna` + `pnc` + `png` + `pnt`))).

    Parameters
    ----------
    ploq: float
        Probability that a base is low-quality.
    pam: float
        Probability that an A is mutated.
    pac: float
        Probability that a mutated A is a substitution to C.
    pag: float
        Probability that a mutated A is a substitution to G.
    pat: float
        Probability that a mutated A is a substitution to T.
    pcm: float
        Probability that a C is mutated.
    pca: float
        Probability that a mutated C is a substitution to A.
    pcg: float
        Probability that a mutated C is a substitution to G.
    pct: float
        Probability that a mutated C is a substitution to T.
    pgm: float
        Probability that a G is mutated.
    pga: float
        Probability that a mutated G is a substitution to A.
    pgc: float
        Probability that a mutated G is a substitution to C.
    pgt: float
        Probability that a mutated G is a substitution to T.
    ptm: float
        Probability that a T is mutated.
    pta: float
        Probability that a mutated T is a substitution to A.
    ptc: float
        Probability that a mutated T is a substitution to C.
    ptg: float
        Probability that a mutated T is a substitution to G.
    pnm: float
        Probability that an N is mutated.
    pnd: float
        Probability that a mutated N is a deletion.
    """
    if not 0. <= ploq <= 1.:
        raise ValueError(f"ploq must be ≥ 0 and ≤ 1, but got {ploq}")
    phiq = 1. - ploq
    # Mutations at A bases.
    pas = [pac, pag, pat]
    verify_proportions(pas + [(pad := 1. - sum(pas))])
    a = pd.Series({SUB_C: pam * phiq * pac,
                   SUB_G: pam * phiq * pag,
                   SUB_T: pam * phiq * pat,
                   DELET: pam * pad,
                   ANY_B: (1. - pam * pad) * ploq})
    # Mutations at C bases.
    pcs = [pca, pcg, pct]
    verify_proportions(pas + [(pcd := 1. - sum(pcs))])
    c = pd.Series({SUB_A: pcm * phiq * pca,
                   SUB_G: pcm * phiq * pcg,
                   SUB_T: pcm * phiq * pct,
                   DELET: pcm * pcd,
                   ANY_D: (1. - pcm * pcd) * ploq})
    # Mutations at G bases.
    pgs = [pga, pgc, pgt]
    verify_proportions(pgs + [(pgd := 1. - sum(pgs))])
    g = pd.Series({SUB_A: pgm * phiq * pga,
                   SUB_C: pgm * phiq * pgc,
                   SUB_T: pgm * phiq * pgt,
                   DELET: pgm * pgd,
                   ANY_H: (1. - pgm * pgd) * ploq})
    # Mutations at T bases.
    pts = [pta, ptc, ptg]
    verify_proportions(pgs + [(ptd := 1. - sum(pts))])
    t = pd.Series({SUB_A: ptm * phiq * pta,
                   SUB_C: ptm * phiq * ptc,
                   SUB_G: ptm * phiq * ptg,
                   DELET: ptm * ptd,
                   ANY_V: (1. - ptm * ptd) * ploq})
    # Mutations at N bases.
    n = pd.Series({DELET: pnm * ptd,
                   ANY_N: 1. - pnm * pnd})
    return pd.DataFrame.from_dict({BASEA: a,
                                   BASEC: c,
                                   BASEG: g,
                                   BASET: t,
                                   BASEN: n}).fillna(0.)


def make_pmut_means_paired(pam: float = 0.004,
                           pcm: float = 0.003,
                           pgm: float = 0.005,
                           ptm: float = 0.002,
                           pnm: float = 0.005,
                           **kwargs):
    """ Generate mean mutation rates for paired bases. """
    return make_pmut_means(pam=pam,
                           pcm=pcm,
                           pgm=pgm,
                           ptm=ptm,
                           pnm=pnm,
                           **kwargs)


def make_pmut_means_unpaired(pam: float = 0.040,
                             pcm: float = 0.030,
                             pgm: float = 0.005,
                             ptm: float = 0.002,
                             pnm: float = 0.005,
                             **kwargs):
    """ Generate mean mutation rates for unpaired bases. """
    return make_pmut_means(pam=pam,
                           pcm=pcm,
                           pgm=pgm,
                           ptm=ptm,
                           pnm=pnm,
                           **kwargs)


def sim_pmut(is_paired: pd.Series,
             pm: pd.DataFrame,
             um: pd.DataFrame,
             pv: float,
             uv: float):
    """ Simulate mutation rates using two Dirichlet distributions for
    the paired and unpaired bases.

    Parameters
    ----------
    is_paired: pd.Series
        Whether each base is paired.
    pm: pd.DataFrame
        Mean of the mutation rates for each type of paired base.
    um: pd.DataFrame
        Mean of the mutation rates for each type of unpaired base.
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
    if not isinstance(is_paired, pd.Series):
        raise TypeError(f"is_paired must be a Series, "
                        f"but got {type(is_paired).__name__}")
    if not isinstance(pm, pd.DataFrame):
        raise TypeError(f"pm must be a DataFrame, but got {type(pm).__name__}")
    if not isinstance(um, pd.DataFrame):
        raise TypeError(f"um must be a DataFrame, but got {type(um).__name__}")
    if pm.values.min() < 0.:
        raise ValueError(f"All pm must be ≥ 0, but got {pm}")
    if um.values.min() < 0.:
        raise ValueError(f"All um must be ≥ 0, but got {um}")
    is_paired = is_paired.astype(bool, copy=False)
    bases = is_paired.index.get_level_values(BASE_NAME)
    # Determine the types of relationships.
    rels = pd.Index.union(pm.index, um.index).astype(REL_TYPE, copy=False)
    if MATCH in rels:
        raise ValueError(f"Matches cannot be a relationship to simulate")
    if NOCOV in rels:
        raise ValueError(f"No coverage cannot be a relationship to simulate")
    rels = rels.append(pd.Index([MATCH]))
    # Copy the mean mutation rates to prevent the originals from being
    # modified, and set their indexes to that of all relationships.
    pm = pm.reindex(index=rels, columns=DNA.alph(), fill_value=0.)
    um = um.reindex(index=rels, columns=DNA.alph(), fill_value=0.)
    # Simulate matches as the default category.
    if ((pm_sum := pm.sum(axis=0)) > 1.).any():
        raise ValueError(f"The sum of pm must be ≤ 1, but got {pm_sum}")
    if ((um_sum := um.sum(axis=0)) > 1.).any():
        raise ValueError(f"The sum of um must be ≤ 1, but got {um_sum}")
    pm.loc[MATCH] = 1. - pm_sum
    um.loc[MATCH] = 1. - um_sum
    # Simulate mutation rates for paired/unpaired bases of each kind.
    pmut = pd.DataFrame(0., index=is_paired.index, columns=rels)
    for base in DNA.alph():
        # Determine which bases of this kind are paired, and count them.
        base_is_paired = is_paired.loc[bases == base]
        if base_is_paired.size == 0:
            continue
        base_index = base_is_paired.index
        base_num_paired = np.count_nonzero(base_is_paired)
        base_num_unpaired = base_is_paired.size - base_num_paired
        # Determine which mean mutation rates are not zero.
        pm_nonzero = pm.loc[pm[base] != 0., base]
        um_nonzero = um.loc[um[base] != 0., base]
        pv_nonzero = pv * (pm_nonzero * (1. - pm_nonzero))
        uv_nonzero = uv * (um_nonzero * (1. - um_nonzero))
        # Simulate the mutation rates for the paired/unpaired bases.
        if (num_pv_nz := np.count_nonzero(pv_nonzero)) > 1:

            ppmut = rng.dirichlet(calc_dirichlet_params(pm_nonzero.values,
                                                        pv_nonzero.values),
                                  size=base_num_paired)
        elif num_pv_nz == 1:
            ppmut = rng.beta(*calc_beta_params(pm_nonzero.values[0],
                                               pv_nonzero.values[0]),
                             size=base_num_paired)
        else:
            ppmut = np.broadcast_to(pm_nonzero.values[np.newaxis, :],
                                    (base_num_paired, pm_nonzero.size))
        pmut.loc[base_index[base_is_paired], pm_nonzero.index] = ppmut
        if (num_uv_nz := np.count_nonzero(uv_nonzero)) > 1:

            upmut = rng.dirichlet(calc_dirichlet_params(um_nonzero.values,
                                                        uv_nonzero.values),
                                  size=base_num_unpaired)
        elif num_uv_nz == 1:
            upmut = rng.beta(*calc_beta_params(um_nonzero.values[0],
                                               uv_nonzero.values[0]),
                             size=base_num_unpaired)
        else:
            upmut = np.broadcast_to(um_nonzero.values[np.newaxis, :],
                                    (base_num_unpaired, um_nonzero.size))
        pmut.loc[base_index[~base_is_paired], um_nonzero.index] = upmut
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
