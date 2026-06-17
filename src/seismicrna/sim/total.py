import os
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from . import (
    fastq as fastq_mod,
    fold as fold_mod,
    params as params_mod,
    ref as ref_mod,
    idmut as idmut_mod,
)
from .muts import load_pmut
from ..core import path
from ..core.arg import (
    PROBE_NONE,
    arg_fasta,
    arg_input_path,
    opt_branch,
    opt_batch_size,
    opt_brotli_level,
    opt_ct_file,
    opt_mask_a,
    opt_mask_c,
    opt_mask_g,
    opt_mask_polya,
    opt_mask_u,
    opt_max_fraction_ident,
    opt_max_pearson_sim,
    opt_max_tries,
    opt_min_marcd_sim,
    opt_param_dir,
    opt_write_read_names,
    merge_params,
)
from ..core.logs import logger
from ..core.random import get_random_integer_generator
from ..core.mu.compare import calc_mean_arcsine_distance, calc_pearson
from ..core.rel import RelPattern
from ..core.rna import from_ct
from ..core.run import run_func
from ..core.seq import DNA, parse_fasta, RefRegions
from .idmut import set_sim_mut_params
from ..filter.main import set_mask_acgu, set_mask_polya
from ..idmut.sim import calc_pmut_pattern

COMMAND = __name__.split(os.path.extsep)[-1]


def _clusters_distinct(
    ct_file: Path,
    fasta: str | Path,
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    mask_a: bool,
    mask_c: bool,
    mask_g: bool,
    mask_u: bool,
    mask_polya: int,
    max_fraction_ident: float,
    max_pearson: float,
    min_marcd: float,
) -> bool:
    """Return True if all cluster pairs satisfy the similarity thresholds."""
    try:
        region = next(from_ct(ct_file)).region
    except StopIteration:
        raise ValueError(f"CT file {ct_file} contains 0 structures") from None
    # Mask positions from the coordinates and primers.
    region_coords = list(region_coords)
    region_primers = list(region_primers)
    select_regs = (
        list(
            RefRegions(
                parse_fasta(Path(fasta), DNA),
                ends=region_coords,
                primers=region_primers,
                exclude_primers=True,
            ).regions
        )
        if (region_coords or region_primers)
        else []
    )
    for reg in select_regs:
        if reg.ref == region.ref:
            region.mask_list(range(region.end5, reg.end5 + 1))
            region.mask_list(range(reg.end3, region.end3 + 1))
    # Apply per-nucleotide masking.
    if mask_a:
        region.mask_a()
    if mask_c:
        region.mask_c()
    if mask_g:
        region.mask_g()
    if mask_u:
        region.mask_t()
    region.mask_polya(mask_polya)
    # Compute per-cluster mutation rates over unmasked positions.
    pmut_path = ct_file.with_suffix(path.PARAM_MUTS_EXT)
    mus = calc_pmut_pattern(
        load_pmut(pmut_path), RelPattern.from_counts(count_del=True, count_ins=True)
    )
    mus = mus.loc[region.unmasked_int]
    if not mus.empty:
        # Pairwise checks.
        for col1, col2 in combinations(mus.columns, 2):
            mu1 = mus[col1]
            mu2 = mus[col2]
            if np.mean(np.isclose(mu1, mu2)) >= max_fraction_ident:
                return False
            if float(calc_pearson(mu1, mu2)) >= max_pearson:
                return False
            if (
                min_marcd > 0.0
                and float(calc_mean_arcsine_distance(mu1, mu2)) < min_marcd
            ):
                return False
    return True


@run_func(COMMAND, default=None)
def run(
    *,
    sim_dir: str | Path,
    tmp_pfx: str | Path,
    sample: str,
    refs: str,
    ref: str,
    reflen: int,
    profile_name: str,
    fold_backend: str,
    pseudoknots: bool,
    fold_coords: Iterable[tuple[str, int, int]],
    fold_primers: Iterable[tuple[str, DNA, DNA]],
    fold_regions_file: str | None,
    fold_constraint: str | None,
    fold_temp: float,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_min: int,
    fold_percent: float,
    fold_edelta: float,
    pmut_paired: Iterable[tuple[str, float]],
    pmut_unpaired: Iterable[tuple[str, float]],
    vmut_paired: Iterable[tuple[str, float]],
    vmut_unpaired: Iterable[tuple[str, float]],
    center_fmean: float,
    center_fvar: float,
    length_fmean: float,
    length_fvar: float,
    clust_conc: float,
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    mask_a: bool | None,
    mask_c: bool | None,
    mask_g: bool | None,
    mask_u: bool | None,
    mask_polya: int | None,
    max_fraction_ident: float,
    max_pearson_sim: float,
    min_marcd_sim: float,
    max_tries: int,
    paired_end: bool,
    read_length: int,
    reverse_fraction: float,
    probe: str,
    min_mut_gap_weights: str | None,
    injected_mut_probs: str | None,
    fq_gzip: bool,
    num_reads: int,
    keep_tmp: bool,
    force: bool,
    num_cpus: int,
    seed: int | None,
):
    """Simulate FASTQ files from scratch."""
    seeds = get_random_integer_generator(seed)
    mask_a, mask_c, mask_g, mask_u = set_mask_acgu(
        probe, mask_a, mask_c, mask_g, mask_u
    )
    mask_polya = set_mask_polya(probe, mask_polya)
    min_mut_gap_weights, injected_mut_probs = set_sim_mut_params(
        probe, min_mut_gap_weights, injected_mut_probs
    )
    for attempt in range(1, max(max_tries, 1) + 1):
        logger.debug("Began simulation attempt {} of up to {}", attempt, max_tries)
        # Simulate the reference sequence.
        fasta = str(
            ref_mod.run(
                sim_dir=sim_dir,
                refs=refs,
                ref=ref,
                reflen=reflen,
                force=force,
                seed=next(seeds),
            )
        )
        # Simulate the structures.
        ct_files = fold_mod.run(
            fasta=fasta,
            sim_dir=sim_dir,
            tmp_pfx=tmp_pfx,
            profile_name=profile_name,
            fold_backend=fold_backend,
            pseudoknots=pseudoknots,
            probe=probe,
            fold_coords=fold_coords,
            fold_primers=fold_primers,
            fold_regions_file=fold_regions_file,
            fold_constraint=fold_constraint,
            fold_temp=fold_temp,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_min=fold_min,
            fold_percent=fold_percent,
            fold_edelta=fold_edelta,
            keep_tmp=keep_tmp,
            force=force,
            num_cpus=num_cpus,
        )
        # Check whether folding succeeded.
        if len(ct_files) == 0:
            # If not, delete the simulated FASTA and start over.
            Path(fasta).unlink(missing_ok=True)
            logger.warning(
                "Got fewer than fold_min={} clusters on attempt {} of up to {}",
                fold_min,
                attempt,
                max_tries,
            )
            # Increase fold_edelta in case that parameter is limiting
            # the number of structures.
            fold_edelta += 1.0  # kcal/mol
            continue
        elif len(ct_files) > 1:
            raise ValueError(f"Folding produced {ct_files} CT files")
        params_mod.run(
            ct_file=ct_files,
            pmut_paired=pmut_paired,
            pmut_unpaired=pmut_unpaired,
            vmut_paired=vmut_paired,
            vmut_unpaired=vmut_unpaired,
            probe=probe,
            region_coords=region_coords,
            region_primers=region_primers,
            center_fmean=center_fmean,
            center_fvar=center_fvar,
            length_fmean=length_fmean,
            length_fvar=length_fvar,
            clust_conc=clust_conc,
            force=force,
            num_cpus=num_cpus,
            seed=next(seeds),
        )
        # Check if the clusters are sufficiently distinct, unless probe
        # is none, in which case there should be no distinct clusters.
        ct_file = ct_files[0]
        if probe != PROBE_NONE and not _clusters_distinct(
            ct_file,
            fasta,
            region_coords,
            region_primers,
            mask_a,
            mask_c,
            mask_g,
            mask_u,
            mask_polya,
            max_fraction_ident,
            max_pearson_sim,
            min_marcd_sim,
        ):
            # If so, delete the simulated files and start over.
            Path(fasta).unlink(missing_ok=True)
            ct_file.unlink(missing_ok=True)
            for suffix in [
                path.PARAM_MUTS_EXT,
                path.PARAM_ENDS_EXT,
                path.PARAM_CLUSTS_EXT,
            ]:
                ct_file.with_suffix(suffix).unlink(missing_ok=True)
            logger.warning(
                "Clusters were too similar on attempt {} of up to {}",
                attempt,
                max_tries,
            )
            continue
        # Successfully generated structures and parameters.
        break
    else:
        # Failed to generate structures and parameters in max_attempts.
        raise RuntimeError(
            "Failed to simulate references, structures, and parameters within "
            f"{max_tries} attempt(s)"
        )
    return fastq_mod.run(
        input_path=(),
        param_dir=[ct_file.parent],
        profile_name=profile_name,
        sample=sample,
        paired_end=paired_end,
        read_length=read_length,
        reverse_fraction=reverse_fraction,
        probe=probe,
        min_mut_gap_weights=min_mut_gap_weights,
        injected_mut_probs=injected_mut_probs,
        fq_gzip=fq_gzip,
        num_reads=num_reads,
        force=force,
        num_cpus=num_cpus,
        seed=next(seeds),
    )


params = merge_params(
    [
        opt_max_tries,
        opt_mask_a,
        opt_mask_c,
        opt_mask_g,
        opt_mask_u,
        opt_mask_polya,
        opt_max_fraction_ident,
        opt_max_pearson_sim,
        opt_min_marcd_sim,
    ],
    fastq_mod.params,
    fold_mod.params,
    params_mod.params,
    ref_mod.params,
    idmut_mod.params,
    exclude=[
        arg_fasta,
        arg_input_path,
        opt_branch,
        opt_batch_size,
        opt_brotli_level,
        opt_ct_file,
        opt_param_dir,
        opt_write_read_names,
    ],
)


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Simulate FASTQ files from scratch."""
    return run(*args, **kwargs)
