from collections import defaultdict
from functools import cached_property
from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from . import (mask as mask_mod,
               cluster as cluster_mod,
               join as join_mod)
from .cluster.data import ClusterMutsDataset, get_clust_params
from .cluster.emk import calc_mean_arcsine_distance_clusters
from .cluster.report import ClusterReport
from .core import path
from .core.arg import (CMD_ENSEMBLES,
                       CMD_JOIN,
                       merge_params,
                       extra_defaults,
                       opt_join_clusts,
                       opt_region_length,
                       opt_region_min_overlap,
                       opt_max_marcd_join)
from .core.dataset import MutsDataset
from .core.error import (IncompatibleValuesError,
                         InconsistentValueError,
                         OutOfBoundsError)
from .core.logs import logger
from .core.rel import RelPattern
from .core.report import KsWrittenF, End5F, End3F
from .core.run import run_func
from .core.seq import DNA, DuplicateReferenceNameError
from .core.task import dispatch
from .join import joined_mask_report_exists, join_regions
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
        raise IncompatibleValuesError("Must have 1 ≤ total_end5 ≤ total_end3, "
                                      f"but got total_end5={total_end5} "
                                      f"and total_end3={total_end3}")
    total_length = total_end3 - total_end5 + 1
    assert total_length >= 1
    if not isinstance(region_length, int):
        raise TypeError(region_length)
    if region_length < 1:
        raise OutOfBoundsError(
            f"region_length must be ≥ 1, but got {region_length}"
        )
    if region_length > total_length:
        logger.warning(f"region_length ({region_length}) is greater than "
                       f"total length of region ({total_length}): "
                       f"using region_length of {total_length}")
        return [(total_end5, total_end3)]
    assert 1 <= region_length <= total_length <= total_end3
    if not isinstance(region_min_overlap, float):
        raise TypeError(region_min_overlap)
    if not 0. < region_min_overlap < 1.:
        raise OutOfBoundsError("region_min_overlap must be > 0 and < 1, "
                               f"but got {region_min_overlap}")
    max_step_size = int(region_length * (1. - region_min_overlap))
    assert 0 <= max_step_size < region_length
    if max_step_size == 0:
        raise IncompatibleValuesError(
            f"Cannot have region_length={region_length} "
            f"with region_min_overlap={region_min_overlap}"
        )
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


class CalcRefRegionLengthError(ValueError):
    """ Error when calculating mutation densities. """


def calc_ref_region_length(datasets: Iterable[MutsDataset],
                           pattern: RelPattern,
                           mask_discontig: bool,
                           min_mut_gap: int):
    logger.routine("Began calculating optimal region length")
    ref = None
    total_muts = 0
    total_cover = 0
    for dataset in datasets:
        if ref is None:
            ref = dataset.ref
        elif dataset.ref != ref:
            raise InconsistentValueError(
                f"Got multiple references: {repr(ref)} ≠ {repr(dataset.ref)}"
            )
        for batch in dataset.iter_batches():
            description = f"Reference {repr(ref)} batch {batch.batch}"
            # Ignore reads that are discontiguous or have two mutations
            # too close.
            logger.detail(f"{description} has {batch.num_reads} reads")
            use_reads = np.ones(batch.num_reads, dtype=bool)
            if mask_discontig:
                logger.detail(f"Dropped {batch.num_discontiguous} "
                              f"discontiguous reads from {description}")
                use_reads = np.logical_and(use_reads, batch.contiguous)
            if min_mut_gap > 0:
                min_mut_dist = batch.calc_min_mut_dist(pattern)
                valid_mut_dist = np.logical_or(min_mut_dist == 0,
                                               min_mut_dist > min_mut_gap)
                num_invalid_mut_dist = np.count_nonzero(
                    np.logical_and(use_reads, ~valid_mut_dist)
                )
                logger.detail(f"Dropped {num_invalid_mut_dist} reads with two "
                              f"mutations too close from {description}")
                use_reads = np.logical_and(use_reads, valid_mut_dist)
            logger.detail(f"Using {np.count_nonzero(use_reads)} reads "
                          f"in {description}")
            # Count the covered and mutated positions among all reads.
            _, muts = batch.count_per_read(pattern)
            batch_muts = int(muts.values[use_reads].sum())
            batch_cover = int(batch.cover_per_read.values[use_reads].sum())
            total_muts += batch_muts
            total_cover += batch_cover
            logger.detail(f"{description} has {batch_muts} mutations "
                          f"among {batch_cover} total bases")
    logger.detail(f"Reference {repr(ref)} has {total_muts} mutations "
                  f"among {total_cover} total bases")
    if total_cover == 0:
        raise CalcRefRegionLengthError(
            f"Got 0 base calls for reference {repr(ref)}"
        )
    if total_muts == 0.:
        raise CalcRefRegionLengthError(
            f"Got 0 mutations for reference {repr(ref)}"
        )
    # The length of a read expected to contain 2 mutations equals 2
    # divided by the density of mutations (the number of mutations
    # divided by the number of base calls). Make sure it is ≥ 1.
    region_length = max(int(np.ceil(2. * total_cover / total_muts)), 1)
    logger.routine(f"Ended calculating optimal region length: {region_length}")
    return region_length


def generate_regions(input_path: Iterable[str | Path],
                     coords: Iterable[tuple[str, int, int]],
                     primers: Iterable[tuple[str, DNA, DNA]],
                     primer_gap: int,
                     regions_file: str | None,
                     region_length: int,
                     region_min_overlap: float,
                     mask_del: bool,
                     mask_ins: bool,
                     mask_mut: list[str],
                     mask_discontig: bool,
                     min_mut_gap: int):
    """ For each reference, list the regions over which to mask. """
    pattern = RelPattern.from_counts(not mask_del, not mask_ins, mask_mut)
    datasets, total_regions = load_regions(input_path,
                                           coords,
                                           primers,
                                           primer_gap,
                                           regions_file)
    mask_regions = list()
    for ref, ref_total_regions in total_regions.dict.items():
        # Get the region for this reference.
        assert len(ref_total_regions) > 0
        if len(ref_total_regions) > 1:
            raise DuplicateReferenceNameError(ref)
        ref_total_region = ref_total_regions[0]
        if region_length > 0:
            # Use a prespecified region length.
            ref_region_length = region_length
        else:
            # Set the region length to the smallest length such that
            # half of reads are expected to have at least 2 mutations.
            try:
                assert datasets.get(ref)
                ref_region_length = calc_ref_region_length(datasets[ref],
                                                           pattern,
                                                           mask_discontig,
                                                           min_mut_gap)
            except CalcRefRegionLengthError as error:
                logger.warning(error)
                ref_region_length = region_length
        # Calculate and add the regions for this reference.
        ref_regions = calc_regions(ref_total_region.end5,
                                   ref_total_region.end3,
                                   ref_region_length,
                                   region_min_overlap)
        mask_regions.extend((ref, end5, end3) for end5, end3 in ref_regions)
    return mask_regions


class RegionInfo(object):

    def __init__(self,
                 reg: str,
                 end5: int,
                 end3: int,
                 ks: Iterable[int],
                 report_file: Path,
                 verify_times: bool,
                 max_procs: int):
        self.reg = reg
        self.end5 = end5
        self.end3 = end3
        self.ks = sorted(ks)
        self.report_file = report_file
        self.verify_times = verify_times
        self.max_procs = max_procs

    @property
    def ends(self):
        return self.end5, self.end3

    @cached_property
    def clust_params(self):
        return get_clust_params(
            ClusterMutsDataset(self.report_file,
                               verify_times=self.verify_times),
            max_procs=self.max_procs
        )

    def __str__(self):
        return f"Region {repr(self.reg)} with Ks {self.ks}"

    def __repr__(self):
        return str(self)


def group_clusters(cluster_dirs: Iterable[Path],
                   max_marcd_join,
                   verify_times: bool,
                   max_procs: int):
    logger.routine("Began grouping regions")
    # List the fields and 5'/3' ends of each clustered dataset.
    regs_info = defaultdict(list)
    for cluster_dir in cluster_dirs:
        path_fields = path.parse(cluster_dir, path.REG_DIR_SEGS)
        report_file = ClusterReport.build_path(path_fields)
        report = ClusterReport.load(report_file)
        ks = report.get_field(KsWrittenF)
        if ks:
            mask_report = MaskReport.load(MaskReport.build_path(
                path_fields | {path.STEP: path.MASK_STEP}
            ))
            top = path_fields[path.TOP]
            sample = path_fields[path.SAMPLE]
            branches_flat = tuple(path_fields[path.BRANCHES])
            ref = path_fields[path.REF]
            key = top, sample, branches_flat, ref
            reg = path_fields[path.REG]
            end5 = mask_report.get_field(End5F)
            end3 = mask_report.get_field(End3F)
            reg_info = RegionInfo(reg,
                                  end5,
                                  end3,
                                  ks,
                                  report_file,
                                  verify_times,
                                  max_procs)
            regs_info[key].append(reg_info)
            logger.detail(f"Found reference {repr(ref)} {reg_info} "
                          f"for sample {repr(sample)} in {top}")
        else:
            logger.warning(f"Skipped {report_file} with no clusters written")
    # Group the consecutive regions with the same number of clusters
    # if the clusters are similar enough.
    groups = defaultdict(list)
    # Use list(regs_info) because the values of regs_info will change;
    # list(regs_info) does not mutate an object being iterated over.
    for key in list(regs_info):
        # Sort by 5' and 3' ends to ensure regions are joined in order.
        regs_info[key].sort(key=lambda reg_info_: reg_info_.ends)
        top, sample, branches_flat, ref = key
        logger.detail(f"Grouping regions of reference {repr(ref)} "
                      f"sample {repr(sample)} in {top}: {regs_info[key]}")
        # Add the first region to the first group.
        assert regs_info[key]
        reg_info = regs_info[key][0]
        groups[key].append([reg_info])
        logger.detail(f"Added {reg_info} to group {len(groups[key])} of {key}")
        for i in range(1, len(regs_info[key])):
            prev_reg_info = reg_info
            reg_info = regs_info[key][i]
            # Compare the numbers of clusters in the current and the
            # previous regions.
            if reg_info.ks == prev_reg_info.ks:
                # The regions have the same numbers of clusters, so
                # check if the cluster parameters are similar enough.
                # Even if the regions share no positions, their
                # indexes will both include the proportion (0, "p").
                overlap = reg_info.clust_params.index.intersection(
                    prev_reg_info.clust_params.index
                )
                assert overlap.size > 0
                logger.detail(f"{prev_reg_info} and {reg_info} share "
                              f"{overlap.size} parameter(s)")
                # Calculate the difference between the clusters.
                marcd = calc_mean_arcsine_distance_clusters(
                    reg_info.clust_params.loc[overlap].values,
                    prev_reg_info.clust_params.loc[overlap].values
                )
                logger.detail(f"{prev_reg_info} and {reg_info} have a mean "
                              f"arcsine distance of {marcd}")
                if marcd <= max_marcd_join:
                    # The clusters are similar enough, so put them in
                    # the same group.
                    groups[key][-1].append(reg_info)
                else:
                    # The clusters are not similar enough, so put the
                    # regions in different groups.
                    groups[key].append([reg_info])
            else:
                # The regions have different numbers of clusters, so
                # put them in different groups.
                groups[key].append([reg_info])
            logger.detail(
                f"Added {reg_info} to group {len(groups[key])} of {key}"
            )
    logger.routine("Ended grouping regions")
    return {key: [[reg_info.reg for reg_info in key_group]
                  for key_group in key_groups]
            for key, key_groups in groups.items()}


@run_func(CMD_ENSEMBLES, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        # General options
        branch: str,
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
        mask_pos_file: Iterable[str | Path],
        mask_read: Iterable[str],
        mask_read_file: Iterable[str | Path],
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
        min_marcd_run: float,
        max_pearson_run: float,
        max_loglike_vs_best: float,
        min_pearson_vs_best: float,
        max_marcd_vs_best: float,
        try_all_ks: bool,
        write_all_ks: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        # Join options
        joined: str,
        region_length: int,
        region_min_overlap: float,
        max_marcd_join: float):
    """ Infer independent structure ensembles along an entire RNA. """
    if not joined:
        raise ValueError(
            "No prefix for joined regions was given via --joined"
        )
    # Ensure iterable parameters are not exhaustible generators.
    input_path = list(input_path)
    mask_mut = list(mask_mut)
    mask_regions = generate_regions(input_path,
                                    coords=mask_coords,
                                    primers=mask_primers,
                                    primer_gap=primer_gap,
                                    regions_file=mask_regions_file,
                                    region_length=region_length,
                                    region_min_overlap=region_min_overlap,
                                    mask_del=mask_del,
                                    mask_ins=mask_ins,
                                    mask_mut=mask_mut,
                                    mask_discontig=mask_discontig,
                                    min_mut_gap=min_mut_gap)
    mask_dirs = mask_mod.run(
        input_path=input_path,
        branch=branch,
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
        min_marcd_run=min_marcd_run,
        max_pearson_run=max_pearson_run,
        max_loglike_vs_best=max_loglike_vs_best,
        min_pearson_vs_best=min_pearson_vs_best,
        max_marcd_vs_best=max_marcd_vs_best,
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
    cluster_groups = group_clusters(cluster_dirs,
                                    max_marcd_join=max_marcd_join,
                                    verify_times=verify_times,
                                    max_procs=max_procs)
    join_dirs = list()
    # Join the masked regions first, then the clustered regions, because
    # the clustered regions require the masked regions.
    for clustered in [False, True]:
        args = list()
        for key, groups in cluster_groups.items():
            top, sample, branches_flat, ref = key
            for group_num, regs in enumerate(groups, start=1):
                joined_region = f"{joined}{group_num}"
                if clustered or force or not joined_mask_report_exists(
                        top,
                        sample,
                        branches_flat,
                        ref,
                        joined_region,
                        regs
                ):
                    # Every cluster dataset (joined or not) needs a mask
                    # dataset of the same region.
                    args.append((top,
                                 joined_region,
                                 sample,
                                 branches_flat,
                                 ref,
                                 regs,
                                 clustered))
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
                       opt_region_min_overlap,
                       opt_max_marcd_join],
                      exclude=[opt_join_clusts])


@command(CMD_ENSEMBLES, params=params)
def cli(*args, **kwargs):
    """ Infer independent structure ensembles along an entire RNA. """
    return run(*args, **kwargs)
