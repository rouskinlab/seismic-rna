from collections import defaultdict
from itertools import chain
from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from click import command

from . import (mask as mask_mod,
               cluster as cluster_mod,
               join as join_mod)
from .core import path
from .core.arg import (CMD_ENSEMBLES,
                       merge_params,
                       extra_defaults,
                       opt_join_clusts,
                       opt_region_length,
                       opt_region_min_overlap,
                       opt_max_marcd_join)
from .core.batch import (accumulate_confusion_matrices,
                         calc_confusion_pvals,
                         calc_confusion_phi,
                         label_significant_pvals)
from .core.dataset import MutsDataset
from .core.error import (IncompatibleValuesError,
                         OutOfBoundsError)
from .core.logs import logger
from .core.run import run_func
from .core.seq import (POS_NAME,
                       DNA,
                       DuplicateReferenceNameError,
                       RefRegions,
                       unite)
from .core.task import as_list_of_tuples, dispatch
from .mask.dataset import MaskMutsDataset, load_mask_dataset
from .mask.main import load_regions
from .relate.dataset import (load_relate_dataset)


def _group_reports(input_path: Iterable[str | Path]):
    """ Group relate reports by sample and branches. """
    seg_types = load_relate_dataset.report_path_seg_types
    groups = defaultdict(list)
    for relate_report_file in path.find_files_chain(input_path, seg_types):
        fields = path.parse(relate_report_file, seg_types)
        sample = fields[path.SAMPLE]
        branches = tuple(fields[path.BRANCHES])
        groups[sample, branches].append(relate_report_file)
    return groups


def _calc_mask_regions(total_end5: int,
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


def _get_batch_read_lengths(batch_num: int,
                            dataset: MutsDataset):
    batch = dataset.get_batch(batch_num)
    return batch.read_lengths


def _get_dataset_read_lengths(dataset: MutsDataset, *,
                              num_cpus: int):
    read_lengths = dispatch(_get_batch_read_lengths,
                            num_cpus=num_cpus,
                            pass_num_cpus=False,
                            as_list=True,
                            ordered=False,
                            raise_on_error=True,
                            args=as_list_of_tuples(dataset.batch_nums),
                            kwargs=dict(dataset=dataset))
    if sum(a.size for a in read_lengths) == 0:
        raise ValueError(f"{dataset} has 0 reads")
    return np.concatenate(read_lengths, axis=0)


def _calc_ref_region_length(datasets: Iterable[MutsDataset],
                            num_cpus: int):
    """ Calculate twice the median read length. """
    logger.routine("Began calculating optimal region length")
    read_lengths = dispatch(_get_dataset_read_lengths,
                            num_cpus=num_cpus,
                            pass_num_cpus=True,
                            as_list=True,
                            ordered=False,
                            raise_on_error=True,
                            args=as_list_of_tuples(datasets))
    if sum(a.size for a in read_lengths) == 0:
        raise ValueError("Datasets have 0 reads")
    read_lengths = np.concatenate(read_lengths, axis=0)
    median_length = np.median(read_lengths)
    logger.detail(f"The median read length is {median_length}")
    region_length = round(2 * median_length)
    logger.routine(f"Ended calculating optimal region length: {region_length}")
    return region_length


def _calc_ref_mask_regions(ref: str, *,
                           datasets: dict[str, list[MutsDataset]],
                           total_regions: RefRegions,
                           region_length: int,
                           region_min_overlap: float,
                           num_cpus: int):
    logger.routine(f"Began calculating regions for reference {repr(ref)}")
    ref_total_regions = total_regions.dict[ref]
    # Get the region for this reference.
    assert len(ref_total_regions) > 0
    if len(ref_total_regions) > 1:
        raise DuplicateReferenceNameError(ref)
    ref_total_region = ref_total_regions[0]
    if region_length > 0:
        # Use a prespecified region length.
        ref_region_length = region_length
    else:
        # Set the region length to twice the median read length.
        ref_region_length = _calc_ref_region_length(datasets[ref],
                                                    num_cpus=num_cpus)
    # Calculate and add the regions for this reference.
    ref_regions = _calc_mask_regions(ref_total_region.end5,
                                     ref_total_region.end3,
                                     ref_region_length,
                                     region_min_overlap)
    logger.routine(f"Ended calculating regions for reference {repr(ref)}")
    return [(ref, end5, end3) for end5, end3 in ref_regions]


def _generate_tiled_regions(input_path: Iterable[str | Path],
                            coords: Iterable[tuple[str, int, int]],
                            primers: Iterable[tuple[str, DNA, DNA]],
                            primer_gap: int,
                            regions_file: str | None,
                            region_length: int,
                            region_min_overlap: float,
                            num_cpus: int):
    """ For each reference, list the regions over which to mask. """
    datasets, total_regions = load_regions(input_path,
                                           coords,
                                           primers,
                                           primer_gap,
                                           regions_file)
    mask_regions = dispatch(
        _calc_ref_mask_regions,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=False,
        ordered=False,
        raise_on_error=True,
        args=as_list_of_tuples(total_regions.refs),
        kwargs=dict(datasets=datasets,
                    total_regions=total_regions,
                    region_length=region_length,
                    region_min_overlap=region_min_overlap)
    )
    return [region
            for ref_mask_regions in mask_regions
            for region in ref_mask_regions]


def _find_correlated_pairs(dataset: MaskMutsDataset,
                           alpha: float,
                           num_cpus: int):
    # Calculate the confusion matrix for the region.
    n, a, b, ab = accumulate_confusion_matrices(
        dataset.get_batch,
        dataset.num_batches,
        dataset.pattern,
        dataset.region.unmasked.get_level_values(POS_NAME),
        None,
        num_cpus=num_cpus
    )
    # Determine which pairs of positions correlate significantly at a
    # false discovery rate of alpha.
    pvals = calc_confusion_pvals(n, a, b, ab)
    is_significant = label_significant_pvals(pvals, alpha)
    pairs = n.index[is_significant].to_list()
    # Save the confusion matrix as a table and graph.
    table_file = dataset.report_file.with_name("correlation-matrix.csv")
    phi = calc_confusion_phi(n, a, b, ab)
    confusion_matrix = pd.DataFrame.from_dict(
        {"Neither Mutated": n - (a + b - ab),
         "Only A Mutated": a - ab,
         "Only B Mutated": b - ab,
         "Both Mutated": ab,
         "Phi Correlation": phi,
         "Raw P-Value": pvals,
         f"Significant at alpha={alpha}": is_significant},
    )
    confusion_matrix.to_csv(table_file)
    return pairs


def _generate_cluster_intervals(datasets: list[MaskMutsDataset],
                                min_pair_fraction: float,
                                min_cluster_length: int,
                                num_cpus: int,
                                **kwargs):
    # Find pairs of bases that correlate significantly in any region.
    pairs = chain(*dispatch(_find_correlated_pairs,
                            num_cpus=num_cpus,
                            pass_num_cpus=True,
                            as_list=False,
                            ordered=False,
                            raise_on_error=True,
                            args=as_list_of_tuples(datasets),
                            kwargs=kwargs))
    # Count correlated pairs that overlap each position in the region.
    region = unite([dataset.region for dataset in datasets])
    combined_positions = region.unmasked_int
    num_pairs = pd.Series(0, index=combined_positions)
    for end5, end3 in pairs:
        num_pairs.loc[end5: end3] += 1
    # Calculate the maximum number of pairs that could overlap each
    # position (i.e. if all pairs were correlated).
    max_pairs = pd.Series(0, index=combined_positions)
    for dataset in datasets:
        dataset_positions = dataset.region.unmasked_int
        size = dataset_positions.size
        counts = np.arange(1, size + 1) * np.arange(size, 0, -1)
        max_pairs.loc[dataset_positions] += counts
    # Determine for which positions the ratio of correlated pairs to
    # maximum possible pairs is at least min_pair_fraction.
    pair_fraction = num_pairs / max_pairs
    adequate_pair_fraction = pair_fraction >= min_pair_fraction
    # Find the intervals with adequate fractions of correlated pairs.
    boundaries = np.flatnonzero(np.diff(np.concatenate([[False],
                                                        adequate_pair_fraction,
                                                        [False]])))
    assert boundaries.size % 2 == 0
    end5s = combined_positions[boundaries[::2]]
    end3s = combined_positions[boundaries[1::2] - 1]
    # Remove intervals that are too short.
    sufficient_length = end3s - end5s >= min_cluster_length - 1
    end5s = end5s[sufficient_length]
    end3s = end3s[sufficient_length]
    # The end5s and end3s arrays must be converted to Python integers,
    # otherwise they will cause an error during JSON serialization.
    return list(zip(map(int, end5s), map(int, end3s), strict=True))


def _generate_cluster_regions(mask_dirs: list[Path],
                              num_cpus: int,
                              **kwargs):
    # Group the datasets by reference.
    groups = defaultdict(list)
    for dataset in load_mask_dataset.iterate(mask_dirs):
        ref = dataset.ref
        groups[ref].append(dataset)
    # Generate intervals for each group.
    intervals = dispatch(_generate_cluster_intervals,
                         num_cpus=num_cpus,
                         pass_num_cpus=True,
                         as_list=False,
                         ordered=True,
                         raise_on_error=True,
                         args=as_list_of_tuples(groups.values()),
                         kwargs=kwargs)
    # Generate regions for each group.
    return [(ref, end5, end3)
            for ref, ref_intervals in zip(groups.keys(), intervals)
            for end5, end3 in ref_intervals]


def _run_group(relate_report_files: list[Path], *,
               # General options
               branch: str,
               tmp_pfx: str | Path,
               keep_tmp: bool,
               brotli_level: int,
               force: bool,
               num_cpus: int,
               # Ensembles options
               region_length: int,
               region_min_overlap: float,
               alpha: float,
               min_pair_fraction: float,
               min_cluster_length: int,
               # Mask options
               mask_coords: Iterable[tuple[str, int, int]],
               mask_primers: Iterable[tuple[str, DNA, DNA]],
               primer_gap: int,
               mask_regions_file: str | None,
               mask_del: bool,
               mask_ins: bool,
               mask_mut: Iterable[str],
               count_mut: Iterable[str],
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
               min_em_runs: int,
               max_em_runs: int,
               jackpot: bool,
               jackpot_conf_level: float,
               max_jackpot_quotient: float,
               min_em_iter: int,
               max_em_iter: int,
               em_thresh: float,
               min_marcd_run: float,
               max_pearson_run: float,
               max_arcd_vs_ens_avg: float,
               max_gini_run: float,
               max_loglike_vs_best: float,
               min_pearson_vs_best: float,
               max_marcd_vs_best: float,
               try_all_ks: bool,
               write_all_ks: bool,
               cluster_pos_table: bool,
               cluster_abundance_table: bool,
               verify_times: bool):
    """ Run a group of datasets. """
    # Divide the reference into overlapping regions.
    tiled_regions = _generate_tiled_regions(
        relate_report_files,
        coords=mask_coords,
        primers=mask_primers,
        primer_gap=primer_gap,
        regions_file=mask_regions_file,
        region_length=region_length,
        region_min_overlap=region_min_overlap,
        num_cpus=num_cpus
    )
    tiled_dirs = mask_mod.run(
        input_path=relate_report_files,
        branch=branch,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        mask_coords=tuple(tiled_regions),
        mask_primers=(),
        primer_gap=0,
        mask_regions_file=None,
        mask_del=mask_del,
        mask_ins=mask_ins,
        mask_mut=mask_mut,
        count_mut=count_mut,
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
        num_cpus=num_cpus,
        force=force
    )
    # Find regions spanned by correlated base pairs.
    correl_regions = _generate_cluster_regions(
        tiled_dirs,
        alpha=alpha,
        min_pair_fraction=min_pair_fraction,
        min_cluster_length=min_cluster_length,
        num_cpus=num_cpus
    )
    correl_dirs = mask_mod.run(
        input_path=relate_report_files,
        branch=branch,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        mask_coords=tuple(correl_regions),
        mask_primers=(),
        primer_gap=0,
        mask_regions_file=None,
        mask_del=mask_del,
        mask_ins=mask_ins,
        mask_mut=mask_mut,
        count_mut=count_mut,
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
        num_cpus=num_cpus,
        force=force
    )
    # Cluster the regions spanned by correlated base pairs.
    cluster_dirs = cluster_mod.run(
        input_path=correl_dirs,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        min_em_runs=min_em_runs,
        max_em_runs=max_em_runs,
        jackpot=jackpot,
        jackpot_conf_level=jackpot_conf_level,
        max_jackpot_quotient=max_jackpot_quotient,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
        min_marcd_run=min_marcd_run,
        max_pearson_run=max_pearson_run,
        max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
        max_gini_run=max_gini_run,
        max_loglike_vs_best=max_loglike_vs_best,
        min_pearson_vs_best=min_pearson_vs_best,
        max_marcd_vs_best=max_marcd_vs_best,
        try_all_ks=try_all_ks,
        write_all_ks=write_all_ks,
        cluster_pos_table=cluster_pos_table,
        cluster_abundance_table=cluster_abundance_table,
        verify_times=verify_times,
        brotli_level=brotli_level,
        num_cpus=num_cpus,
        force=force,
    )
    return cluster_dirs


@run_func(CMD_ENSEMBLES, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        # General options
        branch: str,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        brotli_level: int,
        force: bool,
        num_cpus: int,
        # Mask options
        mask_coords: Iterable[tuple[str, int, int]],
        mask_primers: Iterable[tuple[str, DNA, DNA]],
        primer_gap: int,
        mask_regions_file: str | None,
        mask_del: bool,
        mask_ins: bool,
        mask_mut: Iterable[str],
        count_mut: Iterable[str],
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
        min_em_runs: int,
        max_em_runs: int,
        jackpot: bool,
        jackpot_conf_level: float,
        max_jackpot_quotient: float,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        min_marcd_run: float,
        max_pearson_run: float,
        max_arcd_vs_ens_avg: float,
        max_gini_run: float,
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
    groups = _group_reports(input_path)
    kwargs = dict(branch=branch,
                  tmp_pfx=tmp_pfx,
                  keep_tmp=keep_tmp,
                  brotli_level=brotli_level,
                  force=force,
                  num_cpus=num_cpus,
                  # Ensembles options
                  region_length=region_length,
                  region_min_overlap=region_min_overlap,
                  alpha=0.05,  # TODO: make this an option
                  min_pair_fraction=0.001,  # TODO: make this an option
                  min_cluster_length=30,  # TODO: make this an option
                  # Mask options
                  mask_coords=mask_coords,
                  mask_primers=mask_primers,
                  primer_gap=primer_gap,
                  mask_regions_file=mask_regions_file,
                  mask_del=mask_del,
                  mask_ins=mask_ins,
                  mask_mut=mask_mut,
                  count_mut=count_mut,
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
                  # Cluster options
                  min_clusters=min_clusters,
                  max_clusters=max_clusters,
                  min_em_runs=min_em_runs,
                  max_em_runs=max_em_runs,
                  jackpot=jackpot,
                  jackpot_conf_level=jackpot_conf_level,
                  max_jackpot_quotient=max_jackpot_quotient,
                  min_em_iter=min_em_iter,
                  max_em_iter=max_em_iter,
                  em_thresh=em_thresh,
                  min_marcd_run=min_marcd_run,
                  max_pearson_run=max_pearson_run,
                  max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
                  max_gini_run=max_gini_run,
                  max_loglike_vs_best=max_loglike_vs_best,
                  min_pearson_vs_best=min_pearson_vs_best,
                  max_marcd_vs_best=max_marcd_vs_best,
                  try_all_ks=try_all_ks,
                  write_all_ks=write_all_ks,
                  cluster_pos_table=cluster_pos_table,
                  cluster_abundance_table=cluster_abundance_table,
                  verify_times=verify_times)
    return list(chain(*dispatch(_run_group,
                                num_cpus=num_cpus,
                                pass_num_cpus=True,
                                as_list=False,
                                ordered=False,
                                raise_on_error=False,
                                args=as_list_of_tuples(groups.values()),
                                kwargs=kwargs)))


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
