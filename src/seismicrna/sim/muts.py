import os
from logging import getLogger
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (docdef,
                        opt_ct_file,
                        opt_pmut_paired,
                        opt_pmut_unpaired,
                        opt_vmut_paired,
                        opt_vmut_unpaired,
                        opt_force,
                        opt_parallel,
                        opt_max_procs)
from ..core.header import RelClustHeader, make_header
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rel import (MATCH,
                        NOCOV,
                        DELET,
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
from ..core.rna import from_ct
from ..core.seq import (BASE_NAME,
                        BASEA,
                        BASEC,
                        BASEG,
                        BASET,
                        BASEN,
                        DNA)
from ..core.stats import calc_beta_params, calc_dirichlet_params
from ..core.write import need_write

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

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


def get_paired(ct_file: Path):
    """ Determine whether each base in paired in one or more structures.

    Parameters
    ----------
    ct_file: DNA
        File of RNA structures in connectivity table (CT) format.

    Returns
    -------
    pd.DataFrame
        Whether each base at each position is paired in each structure.
    """
    return pd.DataFrame.from_dict({
        number: structure.is_paired
        for number, structure in enumerate(from_ct(ct_file), start=1)
    })


def make_pmut_means(*,
                    ploq: float = 0.04,
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
    pmut_means = pd.DataFrame.from_dict({BASEA: a,
                                         BASEC: c,
                                         BASEG: g,
                                         BASET: t,
                                         BASEN: n}).fillna(0.)
    pmut_means.loc[MATCH] = 1. - pmut_means.sum(axis=0)
    return pmut_means


def make_pmut_means_paired(pam: float = 0.004,
                           pcm: float = 0.003,
                           pgm: float = 0.003,
                           ptm: float = 0.001,
                           pnm: float = 0.002,
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
                             pgm: float = 0.003,
                             ptm: float = 0.001,
                             pnm: float = 0.002,
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
    if pm.values.min(initial=1.) < 0.:
        raise ValueError(f"All pm must be ≥ 0, but got {pm}")
    if um.values.min(initial=1.) < 0.:
        raise ValueError(f"All um must be ≥ 0, but got {um}")
    is_paired = is_paired.astype(bool, copy=False)
    bases = is_paired.index.get_level_values(BASE_NAME)
    # Determine the types of relationships.
    rels = pd.Index.union(pm.index, um.index).astype(REL_TYPE, copy=False)
    if MATCH not in rels:
        raise ValueError(f"Relationships omit matches ({MATCH}): {rels}")
    if NOCOV in rels:
        raise ValueError(f"Relationships include no coverage ({NOCOV}): {rels}")
    # Copy the mean mutation rates to prevent the originals from being
    # modified, and set their indexes to that of all relationships.
    pm = pm.reindex(index=rels, columns=DNA.alph(), fill_value=0.)
    um = um.reindex(index=rels, columns=DNA.alph(), fill_value=0.)
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
    return pmut


def _make_pmut_means_kwargs(pmut: tuple[tuple[str, float], ...]):
    """ Make keyword arguments for `make_pmut_means`. """
    return {f"p{mut}": p for mut, p in pmut}


def run_struct(ct_file: Path,
               pmut_paired: tuple[tuple[str, float], ...],
               pmut_unpaired: tuple[tuple[str, float], ...],
               vmut_paired: float,
               vmut_unpaired: float,
               force: bool):
    pmut_file = ct_file.with_suffix(path.PARAM_MUTS_EXT)
    if need_write(pmut_file, force):
        is_paired = get_paired(ct_file)
        num_structs = is_paired.columns.size
        if num_structs == 0:
            raise ValueError(f"No structures in {ct_file}")
        pm = make_pmut_means_paired(**_make_pmut_means_kwargs(pmut_paired))
        um = make_pmut_means_unpaired(**_make_pmut_means_kwargs(pmut_unpaired))
        # Simulate mutation rates for each structure.
        pmut = {i: sim_pmut(is_paired_i, pm, um, vmut_paired, vmut_unpaired)
                for i, is_paired_i in is_paired.items()}
        # Assemble the mutation rates into one DataFrame.
        rels = list(pmut.values())[0].columns
        header = make_header(rels=map(str, rels),
                             max_order=num_structs,
                             min_order=num_structs)
        pmut_whole = pd.DataFrame(np.nan, is_paired.index, header.index)
        for rel in rels:
            for i, pmut_number in pmut.items():
                pmut_whole[str(rel), num_structs, i] = pmut_number[rel]
        pmut_whole.to_csv(pmut_file)
    return pmut_file


def load_pmut(pmut_file: Path):
    """ Load mutation rates from a file. """
    pmut = pd.read_csv(pmut_file,
                       index_col=list(range(2)),
                       header=list(range(RelClustHeader.num_levels())))
    # Convert the columns from strings to integers.
    pmut.columns = pd.MultiIndex.from_arrays(
        [pmut.columns.get_level_values(level).astype(int)
         for level in pmut.columns.names],
        names=pmut.columns.names
    )
    return pmut


@docdef.auto()
def run(ct_file: tuple[str, ...],
        pmut_paired: tuple[tuple[str, float], ...],
        pmut_unpaired: tuple[tuple[str, float], ...],
        vmut_paired: float,
        vmut_unpaired: float,
        force: bool,
        parallel: bool,
        max_procs: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(run_struct,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(pmut_paired=pmut_paired,
                                pmut_unpaired=pmut_unpaired,
                                vmut_paired=vmut_paired,
                                vmut_unpaired=vmut_unpaired,
                                force=force))


params = [
    opt_ct_file,
    opt_pmut_paired,
    opt_pmut_unpaired,
    opt_vmut_paired,
    opt_vmut_unpaired,
    opt_force,
    opt_parallel,
    opt_max_procs
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the rate of each kind of mutation at each position. """
    run(*args, **kwargs)
