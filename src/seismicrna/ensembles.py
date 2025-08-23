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

from . import mask as mask_mod, cluster as cluster_mod
from .core import path
from .core.arg import (CMD_ENSEMBLES,
                       GAP_MODE_OMIT,
                       GAP_MODE_INSERT,
                       GAP_MODE_EXPAND,
                       merge_params,
                       extra_defaults,
                       opt_tile_length,
                       opt_tile_min_overlap,
                       opt_erase_tiles,
                       opt_pair_fdr,
                       opt_min_pairs,
                       opt_min_cluster_length,
                       opt_max_cluster_length,
                       opt_gap_mode)
from .core.array import triangular
from .core.batch import (POSITION_A,
                         POSITION_B,
                         accumulate_confusion_matrices,
                         calc_confusion_pvals,
                         calc_confusion_phi,
                         label_significant_pvals)
from .core.dataset import MutsDataset
from .core.error import (IncompatibleValuesError,
                         OutOfBoundsError)
from .core.logs import logger
from .core.run import run_func
from .core.seq import (DNA,
                       DuplicateReferenceNameError,
                       RefRegions,
                       unite)
from .core.task import as_list_of_tuples, dispatch
from .core.validate import require_atleast
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


def _tuples_to_ends_arrays(pairs: list[tuple[int, int]]):
    if not pairs:
        return np.array([], dtype=int), np.array([], dtype=int)
    end5s, end3s = map(np.array, zip(*pairs))
    return end5s, end3s


def _calc_midpoints_distances(end5s: np.ndarray, end3s: np.ndarray):
    assert end5s.shape == end3s.shape
    assert np.all(end5s <= end3s)
    midpoints = (end5s + end3s) / 2
    distances = end3s - end5s
    return midpoints, distances


def _calc_tiles(total_end5: int,
                total_end3: int,
                tile_length: int,
                tile_min_overlap: float):
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
    if not isinstance(tile_length, int):
        raise TypeError(tile_length)
    if tile_length < 1:
        raise OutOfBoundsError(
            f"tile_length must be ≥ 1, but got {tile_length}"
        )
    if tile_length > total_length:
        logger.warning(f"tile_length ({tile_length}) is greater than "
                       f"total length of region ({total_length}): "
                       f"using tile_length of {total_length}")
        return [(total_end5, total_end3)]
    assert 1 <= tile_length <= total_length <= total_end3
    if not isinstance(tile_min_overlap, float):
        raise TypeError(tile_min_overlap)
    if not 0. < tile_min_overlap < 1.:
        raise OutOfBoundsError("tile_min_overlap must be > 0 and < 1, "
                               f"but got {tile_min_overlap}")
    max_step_size = int(tile_length * (1. - tile_min_overlap))
    assert 0 <= max_step_size < tile_length
    if max_step_size == 0:
        raise IncompatibleValuesError(
            f"Cannot have tile_length={tile_length} "
            f"with tile_min_overlap={tile_min_overlap}"
        )
    num_regions = 1 + ceil((total_length - tile_length) / max_step_size)
    region_end5s = np.asarray(
        np.round(np.linspace(total_end5,
                             total_end3 - tile_length + 1,
                             num_regions)),
        dtype=int
    )
    region_end3s = region_end5s + (tile_length - 1)
    return _ends_arrays_to_tuples(region_end5s, region_end3s)


def _calc_ref_tiled_regions(ref: str, *,
                            datasets: dict[str, list[MutsDataset]],
                            total_regions: RefRegions,
                            tile_length: int,
                            tile_min_overlap: float,
                            num_cpus: int):
    """ Calculate the tiles for one reference. """
    logger.routine(f"Began calculating tiles for reference {repr(ref)}")
    ref_total_regions = total_regions.dict[ref]
    # Get the region for this reference.
    assert len(ref_total_regions) > 0
    if len(ref_total_regions) > 1:
        raise DuplicateReferenceNameError(ref)
    ref_total_region = ref_total_regions[0]
    if tile_length > 0:
        # Use a prespecified region length.
        ref_tile_length = tile_length
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
        ref_tile_length = round(2 * median_read_length)
        logger.routine(
            f"Ended calculating optimal tile length: {ref_tile_length}"
        )
    # Calculate the tiles for this reference.
    tiles = _calc_tiles(ref_total_region.end5,
                        ref_total_region.end3,
                        ref_tile_length,
                        tile_min_overlap)
    logger.routine(f"Ended calculating tiles for reference {repr(ref)}")
    return [(ref, end5, end3) for end5, end3 in tiles]


def _calc_tiled_coords(relate_report_files: Iterable[str | Path],
                       coords: Iterable[tuple[str, int, int]],
                       primers: Iterable[tuple[str, DNA, DNA]],
                       primer_gap: int,
                       regions_file: str | None,
                       tile_length: int,
                       tile_min_overlap: float,
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
                    tile_length=tile_length,
                    tile_min_overlap=tile_min_overlap)
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
        dataset.region.unmasked,
        None,
        min_gap=dataset.min_mut_gap,
        num_cpus=num_cpus
    )
    # Determine which pairs of positions correlate significantly at a
    # false discovery rate of pair_fdr.
    pvals = calc_confusion_pvals(n, a, b, ab)
    is_significant = label_significant_pvals(pvals, pair_fdr)
    # Save the confusion matrix as a CSV file.
    csv_file = dataset.report_file.with_name("confusion-matrix.csv")
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
    confusion_matrix.to_csv(csv_file)
    # Return the significantly correlated pairs.
    pairs = n.index
    pairs_significant = pairs[is_significant].to_list()
    logger.detail(f"{dataset} has {len(pairs_significant)} pairs with "
                  f"significant correlations out of {len(pairs)} total pairs")
    return pairs_significant


def _aggregate_pairs(pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """ Aggregate pairs that overlap. """
    if not pairs:
        return list()
    assert len(pairs[0]) == 2
    aggregates = [pairs[0]]
    for pos5, pos3 in pairs[1:]:
        assert 1 <= pos5 <= pos3
        assert pos5 >= aggregates[-1][0]
        if pos5 <= aggregates[-1][1]:
            # This pair overlaps with the previous aggregate.
            if pos3 > aggregates[-1][1]:
                # This pair extends the previous aggregate.
                aggregates[-1] = (aggregates[-1][0], pos3)
        else:
            # This pair does not overlap with the previous aggregate.
            # Create a new aggregate.
            aggregates.append((pos5, pos3))
    return aggregates


def _select_pairs(pairs: list[tuple[int, int]], end5: int, end3: int):
    assert 1 <= end5 <= end3
    return [(pos5, pos3) for pos5, pos3 in pairs
            if end5 <= pos5 and pos3 <= end3]


def _calc_span_per_pos(pairs: list[tuple[int, int]], end5: int, end3: int):
    assert 1 <= end5 <= end3
    spans_per_pos = pd.Series(0, index=range(end5, end3 + 1))
    for pos5, pos3 in pairs:
        assert end5 <= pos5 <= pos3 <= end3
        spans_per_pos.loc[pos5: pos3] += 1
    return spans_per_pos


def _calc_null_span_per_pos_keep_dists(pairs: list[tuple[int, int]],
                                       end5: int,
                                       end3: int):
    assert 1 <= end5 <= end3
    null_spans_per_pos = pd.Series(0., index=range(end5, end3 + 1))
    num_pos = null_spans_per_pos.size
    assert num_pos >= 1
    ramp = np.minimum(np.arange(1, num_pos + 1),
                      np.arange(num_pos, 0, -1))
    fraction_overlap_cache = dict()
    for pos5, pos3 in pairs:
        assert end5 <= pos5 <= pos3 <= end3
        length = pos3 - pos5 + 1
        fraction_overlap = fraction_overlap_cache.get(length)
        if fraction_overlap is None:
            # Number of locations to which the pair could be moved.
            num_loc = num_pos - length + 1
            assert num_loc >= 1
            max_overlap = np.minimum(ramp, min(length, num_loc))
            fraction_overlap = max_overlap / num_loc
            fraction_overlap_cache[length] = fraction_overlap
        null_spans_per_pos += fraction_overlap
    return null_spans_per_pos


def _calc_null_span_per_pos_rand_dists(pairs: list[tuple[int, int]],
                                       end5: int,
                                       end3: int,
                                       min_mut_gap: int):
    assert 1 <= end5 <= end3
    assert min_mut_gap >= 0
    positions = np.arange(end5, end3 + 1)
    num_pos = positions.size
    assert num_pos >= 1
    # Count all intervals (a, b) such that end5 ≤ a ≤ k ≤ b ≤ end3 and
    # b - a > min_mut_gap.
    num_intervals = triangular(num_pos - (min_mut_gap + 1))
    if num_intervals == 0 or len(pairs) == 0:
        return pd.Series(0., index=positions)
    # For each position, count all possible intervals (a, b) such that
    # end5 ≤ a ≤ k ≤ b ≤ end3.
    counts = (positions - (end5 - 1)) * ((end3 + 1) - positions)
    # Remove intervals for which b - a ≤ min_mut_gap.
    ramp = np.minimum(np.arange(1, num_pos + 1),
                      np.arange(num_pos, 0, -1))
    for length in range(1, min(min_mut_gap + 1, num_pos) + 1):
        # Number of possible locations of the interval.
        num_loc = num_pos - length + 1
        counts -= np.minimum(ramp, min(length, num_loc))
    assert np.all(counts >= 0)
    # Calculate the expected number of pairs that overlap each position.
    return pd.Series((len(pairs) / num_intervals) * counts, index=positions)


def _calc_modules_from_pairs(pairs: list[tuple[int, int]],
                             pair_fdr: float,
                             min_mut_gap: int,
                             min_pairs: int = 1,
                             preserve_null_pair_dists: bool = False,
                             max_iter: int = 10000):
    logger.routine("Began calculating modules")
    require_atleast("min_mut_gap", min_mut_gap, 0, classes=int)
    require_atleast("min_pairs", min_pairs, 1, classes=int)
    finished = set()
    # First, naively aggregate all pairs that overlap.
    modules = _aggregate_pairs(pairs)
    for i in range(max_iter):
        logger.detail(f"Modules at iteration {i}: {modules}")
        new_modules = list()
        for module in modules:
            if module in finished:
                new_modules.append(module)
                continue
            end5, end3 = module
            # For each module, count the pairs contained by the module
            # and how many span each position.
            module_pairs = _select_pairs(pairs, end5, end3)
            if len(module_pairs) < min_pairs:
                # Skip modules with insufficient pairs.
                finished.add(module)
                continue
            spans_per_pos = _calc_span_per_pos(module_pairs, end5, end3)
            # Calculate how many pairs expected to span each position
            if preserve_null_pair_dists:
                null_spans_per_pos = _calc_null_span_per_pos_keep_dists(
                    module_pairs,
                    end5,
                    end3,
                )
            else:
                null_spans_per_pos = _calc_null_span_per_pos_rand_dists(
                    module_pairs,
                    end5,
                    end3,
                    min_mut_gap,
                )
            threshold = pair_fdr * null_spans_per_pos
            sufficient_spans_per_pos = spans_per_pos > threshold
            if sufficient_spans_per_pos.all():
                # All positions have enough pairs spanning them: keep
                # the module as is.
                new_modules.append((end5, end3))
            elif sufficient_spans_per_pos.any():
                # Some but not all positions have enough pairs spanning
                # them: split the module where the number of spanning
                # pairs minus the threshold is most negative.
                spans_diff = spans_per_pos - threshold
                split_pos = int(spans_diff.index[np.argmin(spans_diff)])
                if split_pos > end5:
                    new_modules.append((end5, split_pos - 1))
                if split_pos < end3:
                    new_modules.append((split_pos + 1, end3))
            else:
                # No positions have enough pairs spanning them: omit
                # this module.
                pass
            finished.add(module)
        if new_modules == modules:
            break
        modules = new_modules
    else:
        logger.warning(f"Modules did not converge in {max_iter} iterations")
    logger.routine(f"Ended calculating modules: {modules}")
    return modules


def _graph_pairs_and_modules(pairs: list[tuple[int, int]],
                             modules: list[tuple[int, int]],
                             end5: int,
                             end3: int,
                             html_file: str | Path):
    # Create a subplot with two rows: top for pair_fraction,
    # bottom for correlated pairs
    fig = make_subplots(rows=1, cols=1)
    # Graph the modules as triangles.
    end5s, end3s = _tuples_to_ends_arrays(modules)
    modules_midpoints, modules_distances = _calc_midpoints_distances(end5s,
                                                                     end3s)
    for (a, b), x, y in zip(modules,
                            modules_midpoints,
                            modules_distances,
                            strict=True):
        fig.add_trace(go.Scatter(x=[a, b, x, a],
                                 y=[0, 0, y, 0],
                                 mode="none",
                                 fill="toself",
                                 fillcolor=f"rgba(230,159,0,0.5)",
                                 showlegend=False,
                                 name=None,
                                 hoverinfo="text",
                                 text=f"Module {a, b}",
                                 hovertemplate="%{text}<extra></extra>"))
    # Plot the correlated pairs as points.
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    pairs_midpoints, pairs_distances = _calc_midpoints_distances(pos5s,
                                                                 pos3s)
    fig.add_trace(go.Scatter(x=pairs_midpoints,
                             y=pairs_distances,
                             mode="markers",
                             showlegend=False,
                             marker=dict(color="#D55E00"),
                             name=None,
                             hoverinfo="text",
                             text=[f"Pair {pair}" for pair in pairs],
                             hovertemplate="%{text}<extra></extra>"))
    # Finish the layout.
    assert end5 <= end3
    x_range = [end5 - 0.5, end3 + 0.5]
    fig.update_xaxes(title_text="Midpoint Position",
                     showgrid=True,
                     range=x_range)
    if modules_distances.size > 0 or pairs_distances.size > 0:
        y_range = [0.0, 1.05 * np.max(np.concatenate([modules_distances,
                                                      pairs_distances]))]
    else:
        y_range = [0.0, 1.0]
    fig.update_yaxes(title_text="Length of Span",
                     showgrid=True,
                     range=y_range)
    # Save the figure.
    fig.write_html(html_file)


def _insert_modules_into_gaps(modules: list[tuple[int, int]],
                              global_end5: int,
                              global_end3: int):
    """ Turn every gap between modules into a new module. """
    assert 1 <= global_end5 <= global_end3
    new_modules = list()
    prev_pos3 = global_end5 - 1
    # modules is assumed to be sorted.
    for end5, end3 in modules:
        assert global_end5 <= end5 <= end3 <= global_end3
        # Here, no two regions are allowed to overlap.
        assert end5 > prev_pos3
        if end5 > prev_pos3 + 1:
            # There is a gap between this region and the previous.
            # Create a new region to fill the gap.
            new_modules.append(((prev_pos3 + 1), (end5 - 1)))
        new_modules.append((end5, end3))
        prev_pos3 = end3
    if prev_pos3 < global_end3:
        # There is a gap between the last region and total_end3.
        # Create a new region to fill the gap.
        new_modules.append(((prev_pos3 + 1), global_end3))
    return new_modules


def _expand_modules_into_gaps(modules: list[tuple[int, int]],
                              global_end5: int,
                              global_end3: int):
    """ Expand every module to fill gaps on either side. """
    assert 1 <= global_end5 <= global_end3
    if not modules:
        return list()
    # Make mutable end5s/end3s arrays.
    end5s, end3s = map(list, zip(*modules))
    assert 1 <= len(modules) == len(end5s) == len(end3s)
    # Expand to cover the global ends.
    end5s[0] = global_end5
    end3s[-1] = global_end3
    # Split each internal gap: left gets floor, right gets ceil
    for i in range(len(modules) - 1):
        left_end3 = end3s[i]
        right_end5 = end5s[i + 1]
        gap = right_end5 - left_end3 - 1
        if gap > 0:
            expand_left = gap // 2
            expand_right = gap - expand_left
            end3s[i] = left_end3 + expand_left
            end5s[i + 1] = right_end5 - expand_right
    return _ends_arrays_to_tuples(end5s, end3s)


def _filter_modules_length(modules: list[tuple[int, int]],
                           min_length: int | float = 1,
                           max_length: int | float = np.inf):
    """ Remove modules that are too short or too long. """
    return [(end5, end3) for end5, end3 in modules
            if min_length <= (end3 - end5 + 1) <= max_length]


def _write_pairs_to_csv(pairs: list[tuple[int, int]],
                        csv_file: str | Path):
    """ Write the pairs to a CSV file. """
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    df = pd.DataFrame.from_dict({POSITION_A: pos5s, POSITION_B: pos3s},
                                orient="columns")
    df.to_csv(csv_file, index=False)


def _calc_ref_cluster_modules(datasets: list[MaskMutsDataset],
                              pair_fdr: float,
                              min_mut_gap: int,
                              min_pairs: int,
                              min_length: int,
                              max_length: int,
                              gap_mode: str,
                              num_cpus: int):
    """ Calculate the cluster regions for one reference. """
    region = unite([dataset.region for dataset in datasets])
    # Find pairs of bases that correlate significantly in any region.
    # Use sorted(set()) to ensure pairs are unique and sorted.
    pairs = sorted(set(chain(*dispatch(_find_correlated_pairs,
                                       num_cpus=num_cpus,
                                       pass_num_cpus=True,
                                       as_list=False,
                                       ordered=False,
                                       raise_on_error=True,
                                       args=as_list_of_tuples(datasets),
                                       kwargs=dict(pair_fdr=pair_fdr)))))
    # Find modules of correlated pairs.
    modules = _calc_modules_from_pairs(pairs, pair_fdr, min_mut_gap, min_pairs)
    # Determine what to do with gaps between regions.
    if gap_mode == GAP_MODE_INSERT:
        modules = _filter_modules_length(modules, min_length=min_length)
        modules = _insert_modules_into_gaps(modules, region.end5, region.end3)
        modules = _filter_modules_length(modules, max_length=max_length)
    elif gap_mode == GAP_MODE_EXPAND:
        modules = _expand_modules_into_gaps(modules, region.end5, region.end3)
        modules = _filter_modules_length(modules, max_length=max_length)
    elif gap_mode == GAP_MODE_OMIT:
        modules = _filter_modules_length(modules, min_length, max_length)
    else:
        raise ValueError(gap_mode)
    # Write the pairs and modules to CSV files.
    ref_dir = datasets[0].report_file.parent.parent
    _write_pairs_to_csv(pairs, ref_dir.joinpath("pairs.csv"))
    _write_pairs_to_csv(modules, ref_dir.joinpath("modules.csv"))
    # Graph the correlated pairs and modules.
    html_file = ref_dir.joinpath("pairs-and-modules.html")
    try:
        _graph_pairs_and_modules(pairs,
                                 modules,
                                 region.end5,
                                 region.end3,
                                 html_file)
    except Exception as error:
        logger.error(error)
    return modules


def _calc_cluster_modules(mask_dirs: list[Path],
                          num_cpus: int,
                          **kwargs):
    """ Calculate the cluster regions for all references. """
    # Group the datasets by reference.
    groups = defaultdict(list)
    for dataset in load_mask_dataset.iterate(mask_dirs):
        ref = dataset.ref
        groups[ref].append(dataset)
    # Calculate regions for each reference.
    regions = dispatch(_calc_ref_cluster_modules,
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
               tile_length: int,
               tile_min_overlap: float,
               erase_tiles: bool,
               pair_fdr: float,
               min_pairs: int,
               min_cluster_length: int,
               max_cluster_length: int,
               gap_mode: str,
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
        tile_length=tile_length,
        tile_min_overlap=tile_min_overlap,
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
    module_coords = _calc_cluster_modules(
        tiled_dirs,
        pair_fdr=pair_fdr,
        min_mut_gap=min_mut_gap,
        min_pairs=min_pairs,
        min_length=min_cluster_length,
        max_length=max_cluster_length,
        gap_mode=gap_mode,
        num_cpus=num_cpus
    )
    if erase_tiles:
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
    refs_with_correlated_pairs = {ref for ref, end5, end3 in module_coords}
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
    module_dirs = mask_mod.run(
        input_path=relate_report_files_with_correlated_pairs,
        branch=branch,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        mask_coords=tuple(module_coords),
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
        input_path=module_dirs,
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
        # Ensembles options
        tile_length: int,
        tile_min_overlap: float,
        erase_tiles: bool,
        pair_fdr: float,
        min_pairs: int,
        min_cluster_length: int,
        max_cluster_length: int,
        gap_mode: str,
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
                  tile_length=tile_length,
                  tile_min_overlap=tile_min_overlap,
                  erase_tiles=erase_tiles,
                  pair_fdr=pair_fdr,
                  min_pairs=min_pairs,
                  min_cluster_length=min_cluster_length,
                  max_cluster_length=max_cluster_length,
                  gap_mode=gap_mode,
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
                      [opt_tile_length,
                       opt_tile_min_overlap,
                       opt_erase_tiles,
                       opt_pair_fdr,
                       opt_min_pairs,
                       opt_min_cluster_length,
                       opt_max_cluster_length,
                       opt_gap_mode])


@command(CMD_ENSEMBLES, params=params)
def cli(*args, **kwargs):
    """ Infer independent structure ensembles along an entire RNA. """
    return run(*args, **kwargs)
