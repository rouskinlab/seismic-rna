from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Iterable

from click import command

from .write import filter_region
from ..core.arg import (
    CMD_FILTER,
    PROBES,
    PROBE_DMS,
    PROBE_ETC,
    MUT_COLLISIONS_AUTO,
    DEFAULT_MIN_MUT_GAPS,
    DEFAULT_MUT_COLLISIONS,
    arg_input_path,
    opt_tmp_pfx,
    opt_branch,
    opt_keep_tmp,
    opt_probe,
    opt_region_coords,
    opt_region_primers,
    opt_primer_gap,
    opt_regions_file,
    opt_max_filter_iter,
    opt_count_del,
    opt_count_ins,
    opt_no_mut,
    opt_only_mut,
    opt_mask_polya,
    opt_mask_a,
    opt_mask_c,
    opt_mask_g,
    opt_mask_u,
    opt_mask_pos,
    opt_mask_pos_file,
    opt_drop_read,
    opt_drop_read_file,
    opt_drop_discontig,
    opt_min_ncov_read,
    opt_min_fcov_read,
    opt_min_finfo_read,
    opt_max_fmut_read,
    opt_min_mut_gap,
    opt_mut_collisions,
    opt_min_ninfo_pos,
    opt_max_fmut_pos,
    opt_quick_unbias,
    opt_quick_unbias_thresh,
    opt_filter_pos_table,
    opt_filter_read_table,
    opt_brotli_level,
    opt_self_contained,
    opt_num_cpus,
    opt_force,
    optional_path,
)
from ..core.logs import logger
from ..core.run import run_func
from ..core.seq import DNA, RefRegions
from ..core.task import dispatch
from ..idmut.dataset import load_idmut_dataset


def set_mask_acgu(
    probe: str,
    mask_a: bool | None = None,
    mask_c: bool | None = None,
    mask_g: bool | None = None,
    mask_u: bool | None = None,
):
    """Resolve per-base masking flags based on the probe type.

    Parameters
    ----------
    probe: str
        Probe type (one of the values in ``PROBES``), used to set
        defaults when a flag is None.
    mask_a: bool or None, optional
        Whether to mask adenine positions; if None, inferred from
        ``probe``.
    mask_c: bool or None, optional
        Whether to mask cytosine positions; if None, inferred from
        ``probe``.
    mask_g: bool or None, optional
        Whether to mask guanine positions; if None, inferred from
        ``probe``.
    mask_u: bool or None, optional
        Whether to mask uracil/thymine positions; if None, inferred
        from ``probe``.

    Returns
    -------
    tuple[bool, bool, bool, bool]
        Resolved ``(mask_a, mask_c, mask_g, mask_u)`` flags.
    """
    if probe not in PROBES:
        raise ValueError(f"Invalid probe type: {repr(probe)}")
    if mask_a is None:
        mask_a = probe in [PROBE_ETC]
    if mask_c is None:
        mask_c = probe in [PROBE_ETC]
    if mask_g is None:
        mask_g = probe in [PROBE_DMS]
    if mask_u is None:
        mask_u = probe in [PROBE_DMS]
    return mask_a, mask_c, mask_g, mask_u


def set_mut_gap_params(
    probe: str,
    min_mut_gap: int | None = None,
    mut_collisions: str = MUT_COLLISIONS_AUTO,
):
    """Resolve mutation-gap and collision parameters based on the probe type.

    Parameters
    ----------
    probe: str
        Probe type (one of the values in ``PROBES``), used to set
        defaults when a parameter is ``None`` / ``MUT_COLLISIONS_AUTO``.
    min_mut_gap: int or None, optional
        Minimum gap (in nucleotides) between two mutations in the same
        read; if None, a probe-specific default is used.
    mut_collisions: str, optional
        How to handle reads with mutations closer than ``min_mut_gap``;
        if ``MUT_COLLISIONS_AUTO``, a probe-specific default is used.

    Returns
    -------
    tuple[int, str]
        Resolved ``(min_mut_gap, mut_collisions)`` values.
    """
    if min_mut_gap is None:
        min_mut_gap = DEFAULT_MIN_MUT_GAPS[probe]
        logger.trace(f"Auto-selected min_mut_gap={min_mut_gap} for probe {repr(probe)}")
    if mut_collisions == MUT_COLLISIONS_AUTO:
        mut_collisions = DEFAULT_MUT_COLLISIONS[probe]
        logger.trace(
            f"Auto-selected mut_collisions={repr(mut_collisions)} for probe "
            f"{repr(probe)}"
        )
    return min_mut_gap, mut_collisions


def load_regions(
    input_path: Iterable[str | Path],
    coords: Iterable[tuple[str, int, int]],
    primers: Iterable[tuple[str, DNA, DNA]],
    primer_gap: int,
    regions_file: Path | None = None,
):
    """Load regions of idmut reports."""
    # Load all datasets, grouped by their reference names.
    datasets = defaultdict(list)
    for dataset in load_idmut_dataset.iterate(input_path):
        try:
            datasets[dataset.ref].append(dataset)
        except Exception as error:
            logger.error(error)
    # Determine the regions for each reference in the datasets.
    regions = RefRegions(
        {
            (dataset.ref, dataset.refseq)
            for ref_datasets in datasets.values()
            for dataset in ref_datasets
        },
        regs_file=regions_file,
        ends=coords,
        primers=primers,
        primer_gap=primer_gap,
        exclude_primers=True,
    )
    return datasets, regions


@run_func(CMD_FILTER)
def run(
    input_path: Iterable[str | Path],
    *,
    branch: str,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    # Regions
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    primer_gap: int,
    regions_file: str | None,
    # Mutation counting
    count_del: bool,
    count_ins: bool,
    no_mut: Iterable[str],
    only_mut: Iterable[str],
    # Filtering
    probe: str,
    mask_a: bool | None,
    mask_c: bool | None,
    mask_g: bool | None,
    mask_u: bool | None,
    mask_polya: int,
    mask_pos: Iterable[tuple[str, int]],
    mask_pos_file: Iterable[str | Path],
    drop_read: Iterable[str],
    drop_read_file: Iterable[str | Path],
    drop_discontig: bool,
    min_ninfo_pos: int,
    max_fmut_pos: float,
    min_ncov_read: int,
    min_fcov_read: float,
    min_finfo_read: float,
    max_fmut_read: float,
    min_mut_gap: int | None,
    mut_collisions: str,
    # Observer bias correction
    quick_unbias: bool,
    quick_unbias_thresh: float,
    # Iteration
    max_filter_iter: int,
    # Table options
    filter_pos_table: bool,
    filter_read_table: bool,
    # Compression
    brotli_level: int,
    # Self-contained batches
    self_contained: bool,
    # Parallelization
    num_cpus: int,
    # Effort
    force: bool,
) -> list[Path]:
    """Define mutations and regions to filter reads and positions."""
    mask_a, mask_c, mask_g, mask_u = set_mask_acgu(
        probe, mask_a, mask_c, mask_g, mask_u
    )
    min_mut_gap, mut_collisions = set_mut_gap_params(probe, min_mut_gap, mut_collisions)
    datasets, regions = load_regions(
        input_path,
        coords=region_coords,
        primers=region_primers,
        primer_gap=primer_gap,
        regions_file=optional_path(regions_file),
    )
    return dispatch(
        filter_region,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=[
            (dataset, region)
            for ref, ref_datasets in datasets.items()
            for dataset, region in product(ref_datasets, regions.list(ref))
        ],
        kwargs=dict(
            branch=branch,
            tmp_pfx=tmp_pfx,
            keep_tmp=keep_tmp,
            count_del=count_del,
            count_ins=count_ins,
            no_mut=no_mut,
            only_mut=only_mut,
            mask_polya=mask_polya,
            mask_a=mask_a,
            mask_c=mask_c,
            mask_g=mask_g,
            mask_u=mask_u,
            mask_pos=list(mask_pos),
            mask_pos_file=list(mask_pos_file),
            drop_read=list(drop_read),
            drop_read_file=list(drop_read_file),
            drop_discontig=drop_discontig,
            min_ncov_read=min_ncov_read,
            min_fcov_read=min_fcov_read,
            min_finfo_read=min_finfo_read,
            max_fmut_read=max_fmut_read,
            min_mut_gap=min_mut_gap,
            mut_collisions=mut_collisions,
            probe=probe,
            min_ninfo_pos=min_ninfo_pos,
            max_fmut_pos=max_fmut_pos,
            quick_unbias=quick_unbias,
            quick_unbias_thresh=quick_unbias_thresh,
            max_filter_iter=max_filter_iter,
            filter_pos_table=filter_pos_table,
            filter_read_table=filter_read_table,
            brotli_level=brotli_level,
            self_contained=self_contained,
            force=force,
        ),
    )


params = [
    # Input/output paths
    arg_input_path,
    opt_branch,
    opt_tmp_pfx,
    opt_keep_tmp,
    # Regions
    opt_region_coords,
    opt_region_primers,
    opt_primer_gap,
    opt_regions_file,
    # Mutation counting
    opt_count_del,
    opt_count_ins,
    opt_no_mut,
    opt_only_mut,
    # Filtering
    opt_probe,
    opt_mask_a,
    opt_mask_c,
    opt_mask_g,
    opt_mask_u,
    opt_mask_polya,
    opt_mask_pos,
    opt_mask_pos_file,
    opt_min_ninfo_pos,
    opt_max_fmut_pos,
    opt_drop_read,
    opt_drop_read_file,
    opt_drop_discontig,
    opt_min_ncov_read,
    opt_min_fcov_read,
    opt_min_finfo_read,
    opt_max_fmut_read,
    opt_min_mut_gap,
    opt_mut_collisions,
    # Observer bias correction
    opt_quick_unbias,
    opt_quick_unbias_thresh,
    # Iteration
    opt_max_filter_iter,
    # Table options
    opt_filter_pos_table,
    opt_filter_read_table,
    # Compression
    opt_brotli_level,
    # Self-contained batches
    opt_self_contained,
    # Parallelization
    opt_num_cpus,
    # Effort
    opt_force,
]


@command(CMD_FILTER, params=params)
def cli(*args, **kwargs):
    """Define mutations and regions to filter reads and positions."""
    return run(*args, **kwargs)
