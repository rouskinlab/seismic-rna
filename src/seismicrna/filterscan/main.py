from pathlib import Path
from typing import Iterable

from click import command

from .write import filterscan
from .. import filter as filter_mod
from ..core import path
from ..core.arg.cmd import CMD_FILTERSCAN
from ..core.arg.cli import (
    merge_params,
    opt_tile_length,
    opt_tile_min_overlap,
    opt_erase_tiles,
    opt_band_width,
    opt_detect_fdr,
    opt_merge_fdr,
    opt_min_pair_coverage,
    opt_min_expect_both,
    opt_anticorr_only,
    opt_min_domain_length,
    opt_gap_mode,
)
from ..core.run import run_func
from ..core.seq.xna import DNA
from ..core.task import as_list_of_tuples, dispatch
from ..idmut.dataset import load_idmut_dataset


@run_func(CMD_FILTERSCAN)
def run(
    input_path: Iterable[str | Path],
    *,
    # General options
    branch: str,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    brotli_level: int,
    force: bool,
    num_cpus: int,
    # Domain-detection options
    tile_length: int,
    tile_min_overlap: float,
    erase_tiles: bool,
    band_width: int,
    detect_fdr: float,
    merge_fdr: float,
    min_pair_coverage: int,
    min_expect_both: float,
    anticorr_only: bool,
    min_domain_length: int,
    gap_mode: str,
    # Filter options
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    primer_gap: int,
    regions_file: str | None,
    count_del: bool,
    count_ins: bool,
    no_mut: Iterable[str],
    only_mut: Iterable[str],
    probe: str,
    mask_a: bool | None,
    mask_c: bool | None,
    mask_g: bool | None,
    mask_u: bool | None,
    mask_polya: int | None,
    mask_pos: Iterable[tuple[str, int]],
    mask_pos_file: Iterable[str | Path],
    drop_read: Iterable[str],
    drop_read_file: Iterable[str | Path],
    drop_discontig: bool,
    min_ncov_read: int,
    min_fcov_read: float,
    min_finfo_read: float,
    max_fmut_read: float,
    min_mut_gap: int | None,
    mut_collisions: str,
    min_ninfo_pos: int,
    max_fmut_pos: float,
    quick_unbias: bool,
    quick_unbias_thresh: float,
    max_filter_iter: int,
    filter_pos_table: bool,
    filter_read_table: bool,
    self_contained: bool,
):
    """Scan an RNA for domains of correlated base pairs."""
    seg_types = load_idmut_dataset.report_path_seg_types
    idmut_report_files = list(path.find_files_chain(input_path, seg_types))
    kwargs = dict(
        branch=branch,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        brotli_level=brotli_level,
        force=force,
        # Domain-detection options
        tile_length=tile_length,
        tile_min_overlap=tile_min_overlap,
        erase_tiles=erase_tiles,
        band_width=band_width,
        detect_fdr=detect_fdr,
        merge_fdr=merge_fdr,
        min_pair_coverage=min_pair_coverage,
        min_expect_both=min_expect_both,
        anticorr_only=anticorr_only,
        min_domain_length=min_domain_length,
        gap_mode=gap_mode,
        # Filter options
        region_coords=region_coords,
        region_primers=region_primers,
        primer_gap=primer_gap,
        regions_file=regions_file,
        count_del=count_del,
        count_ins=count_ins,
        no_mut=no_mut,
        only_mut=only_mut,
        probe=probe,
        mask_a=mask_a,
        mask_c=mask_c,
        mask_g=mask_g,
        mask_u=mask_u,
        mask_polya=mask_polya,
        mask_pos=mask_pos,
        mask_pos_file=mask_pos_file,
        drop_read=drop_read,
        drop_read_file=drop_read_file,
        drop_discontig=drop_discontig,
        min_ncov_read=min_ncov_read,
        min_fcov_read=min_fcov_read,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        mut_collisions=mut_collisions,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        quick_unbias=quick_unbias,
        quick_unbias_thresh=quick_unbias_thresh,
        max_filter_iter=max_filter_iter,
        filter_pos_table=filter_pos_table,
        filter_read_table=filter_read_table,
        self_contained=self_contained,
    )
    return dispatch(
        filterscan,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=as_list_of_tuples(idmut_report_files),
        kwargs=kwargs,
    )


params = merge_params(
    filter_mod.params,
    [
        opt_tile_length,
        opt_tile_min_overlap,
        opt_erase_tiles,
        opt_band_width,
        opt_detect_fdr,
        opt_merge_fdr,
        opt_min_pair_coverage,
        opt_min_expect_both,
        opt_anticorr_only,
        opt_min_domain_length,
        opt_gap_mode,
    ],
)


@command(CMD_FILTERSCAN, params=params)
def cli(*args, **kwargs):
    """Scan an RNA for domains of correlated base pairs."""
    return run(*args, **kwargs)
