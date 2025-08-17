from collections import defaultdict
from itertools import chain
from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from . import (mask as mask_mod,
               cluster as cluster_mod,
               join as join_mod)
from .core import path
from .core.arg import (CMD_ENSEMBLES,
                       ENSEMBLES_GAP_FILL_NONE,
                       ENSEMBLES_GAP_FILL_INSERT,
                       ENSEMBLES_GAP_FILL_EXPAND,
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
from .mask.io import MaskBatchIO
from .mask.main import load_regions
from .mask.report import MaskReport
from .relate.dataset import load_relate_dataset
from .relate.report import RelateReport


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


def _ends_arrays_to_tuples(end5s: Iterable, end3s: Iterable):
    # The end5s and end3s arrays must be mapped to Python integers,
    # otherwise they will cause an error during JSON serialization.
    return list(zip(map(int, end5s), map(int, end3s), strict=True))


def _calc_tiles(total_end5: int,
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
    return _ends_arrays_to_tuples(region_end5s, region_end3s)


def _calc_ref_tiled_regions(ref: str, *,
                            datasets: dict[str, list[MutsDataset]],
                            total_regions: RefRegions,
                            region_length: int,
                            region_min_overlap: float,
                            num_cpus: int):
    """ Calculate the tiles for one reference. """
    logger.routine(f"Began calculating tiles for reference {repr(ref)}")
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
        logger.routine("Began calculating optimal tile length")
        read_lengths = dispatch(_get_dataset_read_lengths,
                                num_cpus=num_cpus,
                                pass_num_cpus=True,
                                as_list=True,
                                ordered=False,
                                raise_on_error=True,
                                args=as_list_of_tuples(datasets[ref]))
        if sum(a.size for a in read_lengths) == 0:
            raise ValueError("Datasets have 0 reads")
        read_lengths = np.concatenate(read_lengths, axis=0)
        median_read_length = np.median(read_lengths)
        if median_read_length < 1:
            raise ValueError("The median read length must be ≥ 1, "
                             f"but got {median_read_length}")
        logger.detail(f"The median read length is {median_read_length}")
        ref_region_length = round(2 * median_read_length)
        logger.routine(
            f"Ended calculating optimal tile length: {ref_region_length}"
        )
    # Calculate the tiles for this reference.
    tiles = _calc_tiles(ref_total_region.end5,
                        ref_total_region.end3,
                        ref_region_length,
                        region_min_overlap)
    logger.routine(f"Ended calculating tiles for reference {repr(ref)}")
    return [(ref, end5, end3) for end5, end3 in tiles]


def _calc_tiled_coords(relate_report_files: Iterable[str | Path],
                       coords: Iterable[tuple[str, int, int]],
                       primers: Iterable[tuple[str, DNA, DNA]],
                       primer_gap: int,
                       regions_file: str | None,
                       region_length: int,
                       region_min_overlap: float,
                       num_cpus: int):
    """ Calculate the tiled regions for all references. """
    datasets, total_regions = load_regions(relate_report_files,
                                           coords,
                                           primers,
                                           primer_gap,
                                           regions_file)
    tiled_regions = dispatch(
        _calc_ref_tiled_regions,
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
            for ref_tiled_regions in tiled_regions
            for region in ref_tiled_regions]


def _find_correlated_pairs(dataset: MaskMutsDataset,
                           pair_fdr: float,
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
    # false discovery rate of pair_fdr.
    pvals = calc_confusion_pvals(n, a, b, ab)
    is_significant = label_significant_pvals(pvals, pair_fdr)
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
         f"Significant at FDR={pair_fdr}": is_significant},
    )
    confusion_matrix.to_csv(table_file)
    return pairs


def _graph_ref_cluster_regions(combined_positions: np.ndarray,
                               pairs: list[tuple[int, int]],
                               pair_fraction: np.ndarray,
                               min_pair_fraction: float,
                               html_file: str | Path):
    # Create a subplot with two rows: top for pair_fraction, bottom for correlated pairs
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.0,
        row_heights=[0.5, 0.5],
    )
    # Indicate the minimum pair fraction.
    min_pos = combined_positions[0]
    max_pos = combined_positions[-1]
    fig.add_trace(
        go.Scatter(
            x=[min_pos, max_pos],
            y=[min_pair_fraction, min_pair_fraction],
            mode='lines',
            line=dict(color="#56B4E9"),
            name="Minimum Pair Fraction",
            showlegend=True,
        ),
        row=1, col=1
    )
    # Graph the pair fraction.
    fig.add_trace(
        go.Bar(
            x=combined_positions,
            y=pair_fraction,
            showlegend=False,
            marker_color="#0072B2"
        ),
        row=1, col=1
    )
    # Graph the pairs as a scatter plot.
    pos_a, pos_b = list(map(np.array, zip(*pairs)))
    x_scatter = (pos_a + pos_b) / 2  # midpoint
    y_scatter = pos_a - pos_b  # distance
    # Above each pair, add a triangle.
    for (a, b), x, y in zip(pairs, x_scatter, y_scatter, strict=True):
        # Vertices of the triangle
        x_triangle = [a, b, x, a]
        y_triangle = [0, 0, y, 0]
        fig.add_trace(
            go.Scatter(
                x=x_triangle,
                y=y_triangle,
                mode='none',
                fill='toself',
                fillcolor=f'rgba(230,159,0,0.02)',
                showlegend=False,
                hoverinfo='skip',
            ),
            row=2, col=1
        )
    # Plot the correlated pairs as points
    fig.add_trace(
        go.Scatter(
            x=x_scatter,
            y=y_scatter,
            mode='markers',
            showlegend=False,
            marker=dict(color="#D55E00"),
            text=list(map(str, pairs)),
            hovertemplate='%{text}<extra></extra>',
        ),
        row=2, col=1
    )
    # Finish the layout.
    fig.update_xaxes(title_text=None,
                     showgrid=True,
                     row=1, col=1)
    fig.update_xaxes(title_text='Position',
                     showgrid=True,
                     row=2, col=1)
    fig.update_yaxes(title_text='Pair Fraction',
                     row=1, col=1,
                     range=[0.0, max(pair_fraction.max(),
                                     min_pair_fraction) * 1.05])
    fig.update_yaxes(title_text='Correlated Pairs',
                     row=2, col=1,
                     range=[y_scatter.min() * 1.05, 0.0])
    # Save the figure.
    fig.write_html(html_file)


def _insert_regions_into_gaps(ends: list[tuple[int, int]],
                              global_end5: int,
                              global_end3: int):
    """ Turn every gap between regions into a new region. """
    assert 1 <= global_end5 <= global_end3
    new_ends = list()
    prev_end3 = global_end5 - 1
    # ends is assumed to be sorted.
    for end5, end3 in ends:
        assert global_end5 <= end5 <= end3 <= global_end3
        # Here, no two regions are allowed to overlap.
        assert end5 > prev_end3
        if end5 > prev_end3 + 1:
            # There is a gap between this region and the previous.
            # Create a new region to fill the gap.
            new_ends.append(((prev_end3 + 1), (end5 - 1)))
        new_ends.append((end5, end3))
        prev_end3 = end3
    if prev_end3 < global_end3:
        # There is a gap between the last region and total_end3.
        # Create a new region to fill the gap.
        new_ends.append(((prev_end3 + 1), global_end3))
    return new_ends


def _expand_regions_into_gaps(ends: list[tuple[int, int]],
                              global_end5: int,
                              global_end3: int):
    """
    Given sorted, non-overlapping inclusive integer intervals (end5 <= end3),
    expand intervals so that they fill all gaps, while:
      - splitting each internal gap as evenly as possible,
      - extending the first interval to start at total_end5,
      - extending the last interval to end at total_end3,

    Example:
      [(1, 5), (10, 20)] with [1, 20] -> [(1, 7), (8, 20)]
    """
    assert 1 <= global_end5 <= global_end3
    if not ends:
        return list()
    # Make mutable end5/end3 arrays.
    end5s, end3s = list(map(list, zip(*ends)))
    assert 1 <= len(ends) == len(end5s) == len(end3s)
    # Expand to cover the global ends.
    end5s[0] = global_end5
    end3s[-1] = global_end3
    # Split each internal gap: left gets floor, right gets ceil
    for i in range(len(ends) - 1):
        left_end3 = end3s[i]
        right_end5 = end5s[i + 1]
        gap = right_end5 - left_end3 - 1
        if gap > 0:
            expand_left = gap // 2
            expand_right = gap - expand_left
            end3s[i] = left_end3 + expand_left
            end5s[i + 1] = right_end5 - expand_right
    return _ends_arrays_to_tuples(end5s, end3s)


def _calc_ref_cluster_regions(datasets: list[MaskMutsDataset],
                              min_pair_fraction: float,
                              min_cluster_length: int,
                              ensembles_gap_fill: str,
                              num_cpus: int,
                              **kwargs):
    """ Calculate the cluster regions for one reference. """
    # Find pairs of bases that correlate significantly in any region.
    pairs = list(chain(*dispatch(_find_correlated_pairs,
                                 num_cpus=num_cpus,
                                 pass_num_cpus=True,
                                 as_list=False,
                                 ordered=False,
                                 raise_on_error=True,
                                 args=as_list_of_tuples(datasets),
                                 kwargs=kwargs)))
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
        num_pos = dataset_positions.size
        counts = np.arange(1, num_pos + 1) * np.arange(num_pos, 0, -1)
        max_pairs.loc[dataset_positions] += counts
    # Determine for which positions the ratio of correlated pairs to
    # maximum possible pairs is at least min_pair_fraction.
    pair_fractions = num_pairs / max_pairs
    adequate_pair_fraction = pair_fractions >= min_pair_fraction
    # Find the regions with adequate fractions of correlated pairs.
    boundaries = np.flatnonzero(np.diff(np.concatenate([[False],
                                                        adequate_pair_fraction,
                                                        [False]])))
    assert boundaries.size % 2 == 0
    end5s = combined_positions[boundaries[::2]]
    end3s = combined_positions[boundaries[1::2] - 1]
    # Remove regions that are too short.
    sufficient_length = end3s - end5s >= min_cluster_length - 1
    end5s = end5s[sufficient_length]
    end3s = end3s[sufficient_length]
    # Determine what to do with gaps between regions.
    ends = _ends_arrays_to_tuples(end5s, end3s)
    if ensembles_gap_fill == ENSEMBLES_GAP_FILL_INSERT:
        ends = _insert_regions_into_gaps(ends, region.end5, region.end3)
    elif ensembles_gap_fill == ENSEMBLES_GAP_FILL_EXPAND:
        ends = _expand_regions_into_gaps(ends, region.end5, region.end3)
    elif ensembles_gap_fill != ENSEMBLES_GAP_FILL_NONE:
        raise ValueError(ensembles_gap_fill)
    # Save the pair fractions as a CSV.
    ref_dir = datasets[0].report_file.parent.parent
    csv_file = ref_dir.joinpath("pair_fractions.csv")
    pair_data = pd.DataFrame.from_dict({
        "Number of Pairs": num_pairs,
        "Potential Pairs": max_pairs,
        "Fraction of Pairs": pair_fractions,
    })
    pair_data.to_csv(csv_file)
    # Graph the results.
    html_file = ref_dir.joinpath("pair_fractions.html")
    try:
        _graph_ref_cluster_regions(combined_positions,
                                   pairs,
                                   pair_fractions,
                                   min_pair_fraction,
                                   html_file)
    except Exception as error:
        logger.error(error)
    return ends


def _calc_cluster_regions(mask_dirs: list[Path],
                          num_cpus: int,
                          **kwargs):
    """ Calculate the cluster regions for all references. """
    # Group the datasets by reference.
    groups = defaultdict(list)
    for dataset in load_mask_dataset.iterate(mask_dirs):
        ref = dataset.ref
        groups[ref].append(dataset)
    # Calculate regions for each reference.
    regions = dispatch(_calc_ref_cluster_regions,
                       num_cpus=num_cpus,
                       pass_num_cpus=True,
                       as_list=False,
                       ordered=True,
                       raise_on_error=True,
                       args=as_list_of_tuples(groups.values()),
                       kwargs=kwargs)
    # Format the regions so that they can be passed into the mask step.
    return [(ref, end5, end3)
            for ref, ref_regions in zip(groups.keys(), regions)
            for end5, end3 in ref_regions]


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
               pair_fdr: float,
               min_pair_fraction: float,
               min_cluster_length: int,
               ensembles_gap_fill: str,
               delete_tiles: bool,
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
    # Divide each reference into overlapping tiles.
    tiled_coords = _calc_tiled_coords(
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
        mask_coords=tiled_coords,
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
        mask_pos_table=False,
        mask_read_table=False,
        brotli_level=brotli_level,
        num_cpus=num_cpus,
        force=force
    )
    # Find regions spanned by correlated base pairs.
    correl_regions = _calc_cluster_regions(
        tiled_dirs,
        pair_fdr=pair_fdr,
        min_pair_fraction=min_pair_fraction,
        min_cluster_length=min_cluster_length,
        ensembles_gap_fill=ensembles_gap_fill,
        num_cpus=num_cpus
    )
    if delete_tiles:
        # Delete the mask reports and batches of the tiles.
        logger.routine("Began deleting tiled reports and batches")
        for file in path.find_files_chain(tiled_dirs,
                                          MaskReport.get_path_seg_types()):
            file.unlink(missing_ok=True)
        for file in path.find_files_chain(tiled_dirs,
                                          MaskBatchIO.get_path_seg_types()):
            file.unlink(missing_ok=True)
        logger.routine("Ended deleting tiled reports and batches")
    # Drop each relate report whose reference contained no correlated
    # base pairs; otherwise, it will default to clustering the full
    # reference, which is not what seismic ensembles is for.
    refs_with_correlated_pairs = {ref for ref, end5, end3 in correl_regions}
    logger.detail(f"Found {len(refs_with_correlated_pairs)} references with "
                  f"correlated pairs: {sorted(refs_with_correlated_pairs)}")
    relate_report_files_with_correlated_pairs = list()
    for file in relate_report_files:
        ref = path.parse(file, RelateReport.get_path_seg_types())[path.REF]
        if ref in refs_with_correlated_pairs:
            relate_report_files_with_correlated_pairs.append(file)
    logger.detail(f"Found {len(relate_report_files_with_correlated_pairs)} "
                  "relate report files with correlated pairs, out of "
                  f"{len(relate_report_files)} total relate report files")
    # Cluster the regions spanned by correlated base pairs.
    correl_dirs = mask_mod.run(
        input_path=relate_report_files_with_correlated_pairs,
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
        # Ensembles options
        joined: str,
        region_length: int,
        region_min_overlap: float,
        max_marcd_join: float):
    """ Infer independent structure ensembles along an entire RNA. """
    # Group reports by sample and branches.
    seg_types = load_relate_dataset.report_path_seg_types
    groups = defaultdict(list)
    for relate_report_file in path.find_files_chain(input_path, seg_types):
        fields = path.parse(relate_report_file, seg_types)
        sample = fields[path.SAMPLE]
        branches = tuple(fields[path.BRANCHES])
        groups[sample, branches].append(relate_report_file)
    # Process each group separately.
    kwargs = dict(branch=branch,
                  tmp_pfx=tmp_pfx,
                  keep_tmp=keep_tmp,
                  brotli_level=brotli_level,
                  force=force,
                  num_cpus=num_cpus,
                  # Ensembles options
                  region_length=region_length,
                  region_min_overlap=region_min_overlap,
                  pair_fdr=0.05,  # TODO: make this an option
                  min_pair_fraction=0.0005,  # TODO: make this an option
                  min_cluster_length=30,  # TODO: make this an option
                  ensembles_gap_fill="none",  # TODO: make this an option
                  delete_tiles=True,  # TODO: make this an option
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
