import os
from itertools import chain
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (opt_ct_file,
                        opt_region_coords,
                        opt_region_primers,
                        opt_pmut_paired,
                        opt_pmut_unpaired,
                        opt_probe,
                        opt_vmut_paired,
                        opt_vmut_unpaired,
                        opt_force,
                        opt_num_cpus,
                        opt_seed,
                        PROBE_DMS,
                        PROBE_ETC,
                        PROBE_NONE,
                        PROBE_SHAPE,
                        PROBES)
from ..core.error import NoDataError
from ..core.header import RelClustHeader, list_clusts, make_header
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
from ..core.rna import UNPAIRED, find_enclosing_pairs, from_ct
from ..core.run import run_func
from ..core.seq import (BASE_NAME,
                        BASEA,
                        BASEC,
                        BASEG,
                        BASET,
                        BASEN,
                        DNA,
                        POS_NAME,
                        RegionFinder,
                        get_shared_index)
from ..core.random import get_random_integer_generator
from ..core.stats import calc_beta_params, calc_dirichlet_params
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]


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


def make_pmut_means(*,
                    ploq: float,
                    pam: float,
                    pac: float,
                    pag: float,
                    pat: float,
                    pcm: float,
                    pca: float,
                    pcg: float,
                    pct: float,
                    pgm: float,
                    pga: float,
                    pgc: float,
                    pgt: float,
                    ptm: float,
                    pta: float,
                    ptc: float,
                    ptg: float,
                    pnm: float):
    """ Generate mean mutation rates.

    Mutations are assumed to behave as follows:

    - A base ``n`` mutates with probability ``pnm``.

      - If it mutates, then it is a substitution with probability
        (``pna`` + ``pnc`` + ``png`` + ``pnt``).

        - If it is a substitution, then it is high-quality with
          probability (1 - ``ploq``).
        - Otherwise, it is low-quality.

      - Otherwise, it is a deletion.

    - Otherwise, it is low-quality with probability ``ploq``.

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

    Returns
    -------
    pd.DataFrame
        Mean rate of each type of mutation (column) and each base (row).
    """
    if not 0. <= ploq <= 1.:
        raise ValueError(f"ploq must be ≥ 0 and ≤ 1, but got {ploq}")
    # Probability that a base is high-quality.
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
    verify_proportions(pcs + [(pcd := 1. - sum(pcs))])
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
    verify_proportions(pts + [(ptd := 1. - sum(pts))])
    t = pd.Series({SUB_A: ptm * phiq * pta,
                   SUB_C: ptm * phiq * ptc,
                   SUB_G: ptm * phiq * ptg,
                   DELET: ptm * ptd,
                   ANY_V: (1. - ptm * ptd) * ploq})
    # Mutations at N bases.
    n = pd.Series({DELET: pnm,
                   ANY_N: (1. - pnm) * ploq})
    pmut_means = pd.DataFrame.from_dict({BASEA: a,
                                         BASEC: c,
                                         BASEG: g,
                                         BASET: t,
                                         BASEN: n}).fillna(0.)
    pmut_means.loc[MATCH] = 1. - pmut_means.sum(axis=0)
    return pmut_means


_PMUT_MEANS_PAIRED_DEFAULTS = {
    PROBE_DMS: dict(ploq=0.002,
                    pam=0.010, pac=0.32, pag=0.32, pat=0.32,
                    pcm=0.010, pca=0.32, pcg=0.32, pct=0.32,
                    pgm=0.001, pga=0.32, pgc=0.32, pgt=0.32,
                    ptm=0.001, pta=0.32, ptc=0.32, ptg=0.32,
                    pnm=0.001),
    PROBE_SHAPE: dict(ploq=0.002,
                    pam=0.010, pac=0.30, pag=0.30, pat=0.30,
                    pcm=0.010, pca=0.30, pcg=0.30, pct=0.30,
                    pgm=0.010, pga=0.30, pgc=0.30, pgt=0.30,
                    ptm=0.010, pta=0.30, ptc=0.30, ptg=0.30,
                    pnm=0.001),
    PROBE_ETC: dict(ploq=0.002,
                    pam=0.001, pac=0.30, pag=0.30, pat=0.30,
                    pcm=0.001, pca=0.30, pcg=0.30, pct=0.30,
                    pgm=0.010, pga=0.30, pgc=0.30, pgt=0.30,
                    ptm=0.010, pta=0.30, ptc=0.30, ptg=0.30,
                    pnm=0.001),
    PROBE_NONE: dict(ploq=0.002,
                    pam=0.001, pac=0.32, pag=0.32, pat=0.32,
                    pcm=0.001, pca=0.32, pcg=0.32, pct=0.32,
                    pgm=0.001, pga=0.32, pgc=0.32, pgt=0.32,
                    ptm=0.001, pta=0.32, ptc=0.32, ptg=0.32,
                    pnm=0.001),
}


def make_pmut_means_paired(probe: str, **kwargs: float):
    """ Generate mean mutation rates for paired bases. """
    if probe not in PROBES:
        raise ValueError(f"probe must be one of {PROBES}, but got {probe!r}")
    return make_pmut_means(**(_PMUT_MEANS_PAIRED_DEFAULTS[probe] | kwargs))


_PMUT_MEANS_UNPAIRED_DEFAULTS = {
    PROBE_DMS: dict(ploq=0.002,
                    pam=0.045, pac=0.32, pag=0.32, pat=0.32,
                    pcm=0.045, pca=0.32, pcg=0.32, pct=0.32,
                    pgm=0.001, pga=0.32, pgc=0.32, pgt=0.32,
                    ptm=0.001, pta=0.32, ptc=0.32, ptg=0.32,
                    pnm=0.001),
    PROBE_SHAPE: dict(ploq=0.002,
                    pam=0.020, pac=0.30, pag=0.30, pat=0.30,
                    pcm=0.020, pca=0.30, pcg=0.30, pct=0.30,
                    pgm=0.020, pga=0.30, pgc=0.30, pgt=0.30,
                    ptm=0.020, pta=0.30, ptc=0.30, ptg=0.30,
                    pnm=0.001),
    PROBE_ETC: dict(ploq=0.002,
                    pam=0.001, pac=0.30, pag=0.30, pat=0.30,
                    pcm=0.001, pca=0.30, pcg=0.30, pct=0.30,
                    pgm=0.045, pga=0.30, pgc=0.30, pgt=0.30,
                    ptm=0.045, pta=0.30, ptc=0.30, ptg=0.30,
                    pnm=0.001),
    PROBE_NONE: dict(ploq=0.002,
                    pam=0.001, pac=0.32, pag=0.32, pat=0.32,
                    pcm=0.001, pca=0.32, pcg=0.32, pct=0.32,
                    pgm=0.001, pga=0.32, pgc=0.32, pgt=0.32,
                    ptm=0.001, pta=0.32, ptc=0.32, ptg=0.32,
                    pnm=0.001),
}


def make_pmut_means_unpaired(probe: str, **kwargs: float):
    """ Generate mean mutation rates for unpaired bases. """
    if probe not in PROBES:
        raise ValueError(f"probe must be one of {PROBES}, but got {probe!r}")
    return make_pmut_means(**(_PMUT_MEANS_UNPAIRED_DEFAULTS[probe] | kwargs))


_VMUT_PAIRED_DEFAULTS = {
    PROBE_DMS:   dict(a=0.004, c=0.004, g=0.001, t=0.001, n=0.001),
    PROBE_SHAPE: dict(a=0.004, c=0.004, g=0.004, t=0.004, n=0.001),
    PROBE_ETC:   dict(a=0.001, c=0.001, g=0.004, t=0.004, n=0.001),
    PROBE_NONE:  dict(a=0.001, c=0.001, g=0.001, t=0.001, n=0.001),
}

_VMUT_UNPAIRED_DEFAULTS = {
    PROBE_DMS:   dict(a=0.025, c=0.025, g=0.001, t=0.001, n=0.001),
    PROBE_SHAPE: dict(a=0.025, c=0.025, g=0.025, t=0.025, n=0.001),
    PROBE_ETC:   dict(a=0.001, c=0.001, g=0.025, t=0.025, n=0.001),
    PROBE_NONE:  dict(a=0.001, c=0.001, g=0.001, t=0.001, n=0.001),
}


def make_vmut_paired(probe: str, **kwargs: float):
    """ Generate per-base relative variances for paired bases. """
    if probe not in PROBES:
        raise ValueError(f"probe must be one of {PROBES}, but got {probe!r}")
    return _VMUT_PAIRED_DEFAULTS[probe] | kwargs


def make_vmut_unpaired(probe: str, **kwargs: float):
    """ Generate per-base relative variances for unpaired bases. """
    if probe not in PROBES:
        raise ValueError(f"probe must be one of {PROBES}, but got {probe!r}")
    return _VMUT_UNPAIRED_DEFAULTS[probe] | kwargs


def sim_pmut(positions: pd.Index,
             mean: pd.DataFrame,
             relative_variance: dict[str, float],
             end5: int | None,
             end3: int | None,
             seed: int | None):
    """ Simulate mutation rates using a Dirichlet distribution.

    Parameters
    ----------
    positions: pd.Index
        Index of positions and bases.
    mean: pd.DataFrame
        Mean of the mutation rates for each type of base.
    relative_variance: dict[str, float]
        Variance of the mutation rates for each base (keyed by lower-case
        base letter), as a fraction of its supremum.

    Returns
    -------
    pd.DataFrame
        Mutation rates, with the same index as
    """
    rng = np.random.default_rng(seed)
    if not isinstance(positions, pd.MultiIndex):
        raise TypeError(f"positions must be a MultiIndex, "
                        f"but got {type(positions).__name__}")
    if not isinstance(mean, pd.DataFrame):
        raise TypeError(f"mean must be a DataFrame, "
                        f"but got {type(mean).__name__}")
    if mean.values.min(initial=1.) < 0.:
        raise ValueError(f"All mean mutation rates must be ≥ 0, but got {mean}")
    # Determine the types of relationships.
    rels = mean.index.astype(REL_TYPE, copy=False)
    if MATCH not in rels:
        raise ValueError(f"Relationships omit matches ({MATCH}): {rels}")
    if NOCOV in rels:
        raise ValueError(f"Relationships include no coverage ({NOCOV}): {rels}")    
    # Simulate mutation rates for each kind of base.
    pmut = pd.DataFrame(0., index=positions, columns=rels)
    for base in DNA.alph():
        base_pos = positions[positions.get_level_values(BASE_NAME) == base]
        # Determine which mean mutation rates are not zero.
        mean_nonzero = mean.loc[mean[base] != 0., base]
        rel_var = relative_variance[base.lower()]
        var_nonzero = rel_var * (mean_nonzero * (1. - mean_nonzero))
        # Simulate the mutation rates.
        num_nonzero = np.count_nonzero(mean_nonzero)
        if num_nonzero > 1:
            pmut_base = rng.dirichlet(calc_dirichlet_params(mean_nonzero.values,
                                                            var_nonzero.values),
                                      size=base_pos.size)
        elif num_nonzero == 1:
            pmut_base = rng.beta(*calc_beta_params(mean_nonzero.values[0],
                                                   var_nonzero.values[0]),
                                 size=base_pos.size)
        else:
            pmut_base = np.broadcast_to(mean_nonzero.values[np.newaxis, :],
                                        (base_pos.size, mean_nonzero.size))
        pmut.loc[base_pos, mean_nonzero.index] = pmut_base
    if end5 is not None and end3 is not None:
        outside = positions[
            (positions.get_level_values(POS_NAME) < end5) |
            (positions.get_level_values(POS_NAME) > end3)
        ]
        pmut.loc[outside] = 0.
        pmut.loc[outside, MATCH] = 1.
    return pmut


def _make_pmut_means_kwargs(pmut: Iterable[tuple[str, float]]):
    """ Make keyword arguments for `make_pmut_means`. """
    return {f"p{mut}": p for mut, p in pmut}


def _make_vmut_kwargs(vmut: Iterable[tuple[str, float]]):
    """ Make keyword arguments for `make_vmut_paired` / `make_vmut_unpaired`. """
    return {base.lower(): v for base, v in vmut}


def run_struct(ct_file: Path,
               pmut_paired: Iterable[tuple[str, float]],
               pmut_unpaired: Iterable[tuple[str, float]],
               vmut_paired: Iterable[tuple[str, float]],
               vmut_unpaired: Iterable[tuple[str, float]],
               probe: str,
               force: bool,
               seed: int | None,
               region_coords: Iterable[tuple[str, int, int]] = (),
               region_primers: Iterable[tuple[str, DNA, DNA]] = ()):
    """
    Simulate per-position mutation rates for a CT file and write them.

    For each structure in the CT file, mutation rates are simulated
    using a Dirichlet distribution, with separate mean rates for paired
    and unpaired bases.

    Parameters
    ----------
    ct_file: Path
        Path to the connectivity table (CT) file defining structures
        and base-pairing.
    pmut_paired: Iterable[tuple[str, float]]
        Mutation-type/probability pairs for paired bases, passed to
        `make_pmut_means_paired`.
    pmut_unpaired: Iterable[tuple[str, float]]
        Mutation-type/probability pairs for unpaired bases, passed to
        `make_pmut_means_unpaired`.
    vmut_paired: Iterable[tuple[str, float]]
        Per-base relative variance of mutation rates for paired bases
        (each value between 0 and 1), passed to `make_vmut_paired`.
    vmut_unpaired: Iterable[tuple[str, float]]
        Per-base relative variance of mutation rates for unpaired bases
        (each value between 0 and 1), passed to `make_vmut_unpaired`.
    force: bool
        Whether to overwrite an existing output file.
    seed: int | None
        Random seed for reproducibility; None for no fixed seed.

    Returns
    -------
    Path
        Path of the written mutation-rate CSV file.
    """
    pmut_file = ct_file.with_suffix(path.PARAM_MUTS_EXT)
    if need_write(pmut_file, force):
        # Calculate mean mutation rates.
        pm = make_pmut_means_paired(probe, **_make_pmut_means_kwargs(pmut_paired))
        um = make_pmut_means_unpaired(probe, **_make_pmut_means_kwargs(pmut_unpaired))
        # Build per-base variance dicts.
        vp = make_vmut_paired(probe, **_make_vmut_kwargs(vmut_paired))
        vu = make_vmut_unpaired(probe, **_make_vmut_kwargs(vmut_unpaired))
        # Load the structures.
        structures = list(from_ct(ct_file))
        if not structures:
            raise NoDataError(f"{ct_file} contains 0 structures")
        num_structures = len(structures)
        index = get_shared_index(structure.table.index
                                 for structure in structures)
        # Compute the coordinate limits from region_coords/region_primers if given.
        region = structures[0].region
        ref = region.ref
        seq = region.seq
        seq5 = region.end5
        reg_end5 = None
        reg_end3 = None
        region_coords = list(region_coords)
        region_primers = list(region_primers)
        ref_coords = [(e5, e3) for r, e5, e3 in region_coords if r == ref]
        ref_primers = [(fwd, rev) for r, fwd, rev in region_primers if r == ref]
        if ref_coords or ref_primers:
            regs = ([RegionFinder(ref, seq, seq5=seq5, end5=e5, end3=e3)
                     for e5, e3 in ref_coords]
                    + [RegionFinder(ref, seq, seq5=seq5, fwd=fwd, rev=rev,
                                    exclude_primers=True)
                       for fwd, rev in ref_primers])
            reg_end5 = max(r.end5 for r in regs)
            reg_end3 = min(r.end3 for r in regs)
        # For every unique base pair, simulate paired/unpaied mutation
        # rates for the bases enclosed by the pair.
        mu_paired = dict()
        mu_unpaired = dict()
        seeds = get_random_integer_generator(seed)

        def update_mus(pair_: tuple[int, int]):
            """ Simulate mutation rates for paired/unpaired bases that
            are enclosed by a given base pair. """
            end5, end3 = pair_
            if end5 == UNPAIRED == end3:
                use_index = index
            else:
                use_index = index[np.logical_and(
                    index.get_level_values(POS_NAME) >= end5,
                    index.get_level_values(POS_NAME) <= end3
                )]
            mu_paired[pair_] = sim_pmut(use_index, pm, vp,
                                        reg_end5, reg_end3, seed=next(seeds))
            mu_unpaired[pair_] = sim_pmut(use_index, um, vu,
                                          reg_end5, reg_end3, seed=next(seeds))

        unpair = UNPAIRED, UNPAIRED
        update_mus(unpair)
        for pair in set(chain(*[structure.pairs for structure in structures])):
            update_mus(pair)
        # Assemble mutation rates for each structure.
        rels = mu_paired[unpair].columns
        header = make_header(rels=map(str, rels),
                             ks=[num_structures])
        pmut = pd.DataFrame(np.nan, index, header.index)
        for structure, cluster in zip(structures,
                                      list_clusts(num_structures),
                                      strict=True):
            # Find the base pair that encloses each position.
            enclosing = find_enclosing_pairs(structure.table)
            for rel in rels:
                # For each position, select the simulated mutation rates
                # for its enclosing base pair.
                for position, paired in structure.is_paired.items():
                    mu = mu_paired if paired else mu_unpaired
                    pmut.at[position, (str(rel), num_structures, cluster)] = (
                        mu[tuple(enclosing.loc[position])].at[position, rel]
                    )
        pmut.to_csv(pmut_file)
    return pmut_file


def load_pmut(pmut_file: Path):
    """ Load mutation rates from a file. """
    pmut = pd.read_csv(pmut_file,
                       index_col=list(range(2)),
                       header=list(range(RelClustHeader.get_num_levels())))
    # Convert the columns from strings to integers.
    pmut.columns = pd.MultiIndex.from_arrays(
        [pmut.columns.get_level_values(level).astype(int)
         for level in pmut.columns.names],
        names=pmut.columns.names
    )
    return pmut


@run_func(COMMAND)
def run(*,
        ct_file: Iterable[str | Path],
        pmut_paired: Iterable[tuple[str, float]],
        pmut_unpaired: Iterable[tuple[str, float]],
        vmut_paired: Iterable[tuple[str, float]],
        vmut_unpaired: Iterable[tuple[str, float]],
        probe: str,
        region_coords: Iterable[tuple[str, int, int]],
        region_primers: Iterable[tuple[str, DNA, DNA]],
        force: bool,
        num_cpus: int,
        seed: int | None):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(run_struct,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(pmut_paired=pmut_paired,
                                pmut_unpaired=pmut_unpaired,
                                vmut_paired=vmut_paired,
                                vmut_unpaired=vmut_unpaired,
                                probe=probe,
                                region_coords=region_coords,
                                region_primers=region_primers,
                                force=force,
                                seed=seed))


params = [
    opt_ct_file,
    opt_pmut_paired,
    opt_pmut_unpaired,
    opt_vmut_paired,
    opt_vmut_unpaired,
    opt_probe,
    opt_region_coords,
    opt_region_primers,
    opt_force,
    opt_num_cpus,
    opt_seed,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the rate of each kind of mutation at each position. """
    run(*args, **kwargs)
