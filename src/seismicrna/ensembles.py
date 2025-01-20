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
                       opt_join_clusts)
from .core.logs import logger
from .core.report import KsWrittenF, End5F, End3F
from .core.run import run_func
from .core.seq import DNA, DuplicateReferenceNameError
from .core.task import dispatch
from .join import join_regions
from .mask.main import load_regions
from .mask.report import MaskReport


def as_tuple_str(items: Iterable):
    return tuple(map(str, items))


def calc_windows(total_end5: int,
                 total_end3: int,
                 window_size: int,
                 min_overlap: float):
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
    if not isinstance(window_size, int):
        raise TypeError(window_size)
    if window_size < 1:
        raise ValueError(f"window_size must be ≥ 1, but got {window_size}")
    if window_size > total_length:
        logger.warning(f"Window size ({window_size}) is greater than "
                       f"total length of region ({total_length}): "
                       f"using window size of {total_length}")
        return [(total_end5, total_end3)]
    assert 1 <= window_size <= total_length <= total_end3
    if not isinstance(min_overlap, float):
        raise TypeError(min_overlap)
    if not 0. < min_overlap < 1.:
        raise ValueError(
            f"min_overlap must be > 0 and < 1, but got {min_overlap}"
        )
    max_step_size = int(window_size * (1. - min_overlap))
    assert 0 <= max_step_size < window_size
    if max_step_size == 0:
        raise ValueError(f"Cannot have window_size={window_size} "
                         f"with min_overlap={min_overlap}")
    num_windows = 1 + ceil((total_length - window_size) / max_step_size)
    window_end5s = np.asarray(np.round(np.linspace(total_end5,
                                                   total_end3 - window_size + 1,
                                                   num_windows)),
                              dtype=int)
    window_end3s = window_end5s + (window_size - 1)
    return [(int(end5), int(end3))
            for end5, end3 in zip(window_end5s, window_end3s, strict=True)]


def generate_windows(input_path: Iterable[str | Path],
                     coords: tuple[tuple[str, int, int], ...],
                     primers: tuple[tuple[str, DNA, DNA], ...],
                     primer_gap: int,
                     regions_file: str | None,
                     window_size: int,
                     min_overlap: float):
    """ For each reference, list the windows over which to mask. """
    # Load all datasets, grouped by their reference names.
    _, regions = load_regions(input_path,
                              coords,
                              primers,
                              primer_gap,
                              regions_file)
    windows = list()
    for ref, ref_regions in regions.dict.items():
        assert len(ref_regions) > 0
        if len(ref_regions) > 1:
            raise DuplicateReferenceNameError(ref)
        ref_region = ref_regions[0]
        ref_windows = calc_windows(ref_region.end5,
                                   ref_region.end3,
                                   window_size,
                                   min_overlap)
        windows.extend((ref, end5, end3) for end5, end3 in ref_windows)
    return windows


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
def run(input_path: tuple[str, ...], *,
        # General options
        tmp_pfx: str,
        keep_tmp: bool,
        brotli_level: int,
        force: bool,
        max_procs: int,
        # Mask options
        mask_coords: tuple[tuple[str, int, int], ...],
        mask_primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        mask_regions_file: str | None,
        mask_del: bool,
        mask_ins: bool,
        mask_mut: tuple[str, ...],
        mask_polya: int,
        mask_gu: bool,
        mask_pos: tuple[tuple[str, int], ...],
        mask_pos_file: str | None,
        mask_read: tuple[str, ...],
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
        joined: str):
    """ Infer independent structure ensembles along an entire RNA. """
    if not joined:
        raise ValueError(
            "No prefix for joined regions was given via --joined"
        )
    mask_windows = generate_windows(input_path,
                                    coords=mask_coords,
                                    primers=mask_primers,
                                    primer_gap=primer_gap,
                                    regions_file=mask_regions_file,
                                    window_size=180,
                                    min_overlap=2 / 3)
    mask_dirs = mask_mod.run(
        input_path=input_path,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        mask_coords=tuple(mask_windows),
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
        input_path=as_tuple_str(mask_dirs),
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
    results = list()
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
        results.extend(dispatch(join_regions,
                                max_procs=max_procs,
                                pass_n_procs=True,
                                args=args,
                                kwargs=kwargs))
    logger.status(f"Ended {CMD_JOIN}")
    return results


params = merge_params(mask_mod.params,
                      cluster_mod.params,
                      join_mod.params,
                      exclude=[opt_join_clusts])


@command(CMD_ENSEMBLES, params=params)
def cli(*args, **kwargs):
    """ Infer independent structure ensembles along an entire RNA. """
    return run(*args, **kwargs)

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
