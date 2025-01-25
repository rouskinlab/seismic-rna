from collections import defaultdict, namedtuple
from itertools import groupby
from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from . import (mask as mask_mod,
               cluster as cluster_mod,
               join as join_mod)
from .cluster.report import ClusterReport
from .core import path
from .core.arg import (CMD_ENSEMBLES,
                       CMD_JOIN,
                       merge_params,
                       extra_defaults,
                       opt_join_clusts,
                       opt_region_length,
                       opt_region_min_overlap)
from .core.logs import logger
from .core.report import KsWrittenF, End5F, End3F
from .core.run import run_func
from .core.seq import DNA, DuplicateReferenceNameError
from .core.task import dispatch
from .join import join_regions
from .mask.main import load_regions
from .mask.report import MaskReport


def calc_regions(total_end5: int,
                 total_end3: int,
                 region_length: int,
                 region_min_overlap: float):
    if not isinstance(total_end5, int):
        raise TypeError(total_end5)
    if not isinstance(total_end3, int):
        raise TypeError(total_end3)
    if not 1 <= total_end5 <= total_end3:
        raise ValueError("Must have 1 ≤ total_end5 ≤ total_end3, "
                         f"but got total_end5={total_end5} "
                         f"and total_end3={total_end3}")
    total_length = total_end3 - total_end5 + 1
    assert total_length >= 1
    if not isinstance(region_length, int):
        raise TypeError(region_length)
    if region_length < 1:
        raise ValueError(f"region_length must be ≥ 1, but got {region_length}")
    if region_length > total_length:
        logger.warning(f"region_length ({region_length}) is greater than "
                       f"total length of region ({total_length}): "
                       f"using region_length of {total_length}")
        return [(total_end5, total_end3)]
    assert 1 <= region_length <= total_length <= total_end3
    if not isinstance(region_min_overlap, float):
        raise TypeError(region_min_overlap)
    if not 0. < region_min_overlap < 1.:
        raise ValueError("region_min_overlap must be > 0 and < 1, "
                         f"but got {region_min_overlap}")
    max_step_size = int(region_length * (1. - region_min_overlap))
    assert 0 <= max_step_size < region_length
    if max_step_size == 0:
        raise ValueError(f"Cannot have region_length={region_length} "
                         f"with region_min_overlap={region_min_overlap}")
    num_regions = 1 + ceil((total_length - region_length) / max_step_size)
    region_end5s = np.asarray(
        np.round(np.linspace(total_end5,
                             total_end3 - region_length + 1,
                             num_regions)),
        dtype=int
    )
    region_end3s = region_end5s + (region_length - 1)
    return [(int(end5), int(end3))
            for end5, end3 in zip(region_end5s, region_end3s, strict=True)]


def generate_regions(input_path: Iterable[str | Path],
                     coords: Iterable[tuple[str, int, int]],
                     primers: Iterable[tuple[str, DNA, DNA]],
                     primer_gap: int,
                     regions_file: str | None,
                     region_length: int,
                     region_min_overlap: float):
    """ For each reference, list the regions over which to mask. """
    _, total_regions = load_regions(input_path,
                                    coords,
                                    primers,
                                    primer_gap,
                                    regions_file)
    mask_regions = list()
    for ref, ref_total_regions in total_regions.dict.items():
        assert len(ref_total_regions) > 0
        if len(ref_total_regions) > 1:
            raise DuplicateReferenceNameError(ref)
        ref_total_region = ref_total_regions[0]
        ref_regions = calc_regions(ref_total_region.end5,
                                   ref_total_region.end3,
                                   region_length,
                                   region_min_overlap)
        mask_regions.extend((ref, end5, end3) for end5, end3 in ref_regions)
    return mask_regions


RegionInfo = namedtuple("RegionInfo", ["reg", "end5", "end3", "ks"])


def group_clusters(cluster_dirs: Iterable[Path]):
    # List the fields and 5'/3' ends of each clustered dataset.
    clusters = defaultdict(list)
    for cluster_dir in cluster_dirs:
        path_fields = path.parse(cluster_dir, *path.REG_DIR_SEGS)
        report_file = ClusterReport.build_path(**path_fields)
        report = ClusterReport.load(report_file)
        ks = tuple(sorted(report.get_field(KsWrittenF)))
        if ks:
            mask_report = MaskReport.load(MaskReport.build_path(
                **(path_fields | {path.CMD: path.MASK_STEP})
            ))
            top = path_fields[path.TOP]
            sample = path_fields[path.SAMP]
            ref = path_fields[path.REF]
            key = top, sample, ref
            reg = path_fields[path.REG]
            end5 = mask_report.get_field(End5F)
            end3 = mask_report.get_field(End3F)
            clusters[key].append(RegionInfo(reg, end5, end3, ks))
        else:
            logger.warning(f"Skipped {report_file} with no clusters written")
    # Group the consecutive regions with the same number of clusters.
    for key in clusters:
        clusters[key].sort(key=lambda info: (info.end5, info.end3))
    return {key: [[info.reg for info in group]
                  for ks, group in groupby(clusters[key],
                                           key=lambda info: info.ks)]
            for key in clusters}


@run_func(CMD_ENSEMBLES, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        # General options
        tmp_pfx: str | Path,
        keep_tmp: bool,
        brotli_level: int,
        force: bool,
        max_procs: int,
        # Mask options
        mask_coords: Iterable[tuple[str, int, int]],
        mask_primers: Iterable[tuple[str, DNA, DNA]],
        primer_gap: int,
        mask_regions_file: str | None,
        mask_del: bool,
        mask_ins: bool,
        mask_mut: Iterable[str],
        mask_polya: int,
        mask_gu: bool,
        mask_pos: Iterable[tuple[str, int]],
        mask_pos_file: str | None,
        mask_read: Iterable[str],
        mask_read_file: str | None,
        mask_discontig: bool,
        min_ncov_read: int,
        min_finfo_read: float,
        max_fmut_read: float,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        quick_unbias: bool,
        quick_unbias_thresh: float,
        max_mask_iter: int,
        mask_pos_table: bool,
        mask_read_table: bool,
        # Cluster options
        min_clusters: int,
        max_clusters: int,
        em_runs: int,
        jackpot: bool,
        jackpot_conf_level: float,
        max_jackpot_quotient: float,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        min_nrmsd_run: float,
        max_pearson_run: float,
        max_loglike_vs_best: float,
        min_pearson_vs_best: float,
        max_nrmsd_vs_best: float,
        try_all_ks: bool,
        write_all_ks: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        # Join options
        joined: str,
        region_length: int,
        region_min_overlap: float):
    """ Infer independent structure ensembles along an entire RNA. """
    if not joined:
        raise ValueError(
            "No prefix for joined regions was given via --joined"
        )
    # Since input_path is used twice, ensure it is not an exhaustible
    # generator.
    input_path = list(input_path)
    mask_regions = generate_regions(input_path,
                                    coords=mask_coords,
                                    primers=mask_primers,
                                    primer_gap=primer_gap,
                                    regions_file=mask_regions_file,
                                    region_length=region_length,
                                    region_min_overlap=region_min_overlap)
    mask_dirs = mask_mod.run(
        input_path=input_path,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        mask_coords=tuple(mask_regions),
        mask_primers=(),
        primer_gap=0,
        mask_regions_file=None,
        mask_del=mask_del,
        mask_ins=mask_ins,
        mask_mut=mask_mut,
        mask_polya=mask_polya,
        mask_gu=mask_gu,
        mask_pos=mask_pos,
        mask_pos_file=mask_pos_file,
        mask_read=mask_read,
        mask_read_file=mask_read_file,
        mask_discontig=mask_discontig,
        min_ncov_read=min_ncov_read,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        quick_unbias=quick_unbias,
        quick_unbias_thresh=quick_unbias_thresh,
        max_mask_iter=max_mask_iter,
        mask_pos_table=mask_pos_table,
        mask_read_table=mask_read_table,
        brotli_level=brotli_level,
        max_procs=max_procs,
        force=force,
    )
    cluster_dirs = cluster_mod.run(
        input_path=mask_dirs,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        em_runs=em_runs,
        jackpot=jackpot,
        jackpot_conf_level=jackpot_conf_level,
        max_jackpot_quotient=max_jackpot_quotient,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
        min_nrmsd_run=min_nrmsd_run,
        max_pearson_run=max_pearson_run,
        max_loglike_vs_best=max_loglike_vs_best,
        min_pearson_vs_best=min_pearson_vs_best,
        max_nrmsd_vs_best=max_nrmsd_vs_best,
        try_all_ks=try_all_ks,
        write_all_ks=write_all_ks,
        cluster_pos_table=cluster_pos_table,
        cluster_abundance_table=cluster_abundance_table,
        verify_times=verify_times,
        brotli_level=brotli_level,
        max_procs=max_procs,
        force=force,
    )
    logger.status(f"Began {CMD_JOIN}")
    cluster_groups = group_clusters(cluster_dirs)
    join_dirs = list()
    for clustered in [False, True]:
        args = list()
        for key, groups in cluster_groups.items():
            top, sample, ref = key
            for group_num, regs in enumerate(groups, start=1):
                joined_region = f"{joined}{group_num}"
                args.append((top, joined_region, sample, ref, regs, clustered))
        kwargs = dict(clusts=dict(),
                      mask_pos_table=mask_pos_table,
                      mask_read_table=mask_read_table,
                      cluster_pos_table=cluster_pos_table,
                      cluster_abundance_table=cluster_abundance_table,
                      verify_times=verify_times,
                      tmp_pfx=tmp_pfx,
                      keep_tmp=keep_tmp,
                      force=force)
        join_dirs.extend(dispatch(join_regions,
                                  max_procs=max_procs,
                                  pass_n_procs=True,
                                  args=args,
                                  kwargs=kwargs))
    logger.status(f"Ended {CMD_JOIN}")
    return join_dirs


params = merge_params(mask_mod.params,
                      cluster_mod.params,
                      join_mod.params,
                      [opt_region_length,
                       opt_region_min_overlap],
                      exclude=[opt_join_clusts])


@command(CMD_ENSEMBLES, params=params)
def cli(*args, **kwargs):
    """ Infer independent structure ensembles along an entire RNA. """
    return run(*args, **kwargs)
