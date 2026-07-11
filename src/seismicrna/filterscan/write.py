from __future__ import annotations
from datetime import datetime
from itertools import chain
from collections import defaultdict
from math import ceil, inf
from pathlib import Path
from typing import Iterable


from .. import filter as filter_mod
from ..core import path
from ..core.arg.cli import GAP_MODE_OMIT, GAP_MODE_INSERT, GAP_MODE_EXPAND
from ..core.batch.confusion import (
    POSITION_A,
    POSITION_B,
    calc_confusion_pvals,
    calc_confusion_phi,
    calc_bh_adjusted_pvals,
)
from ..core.batch.accum import accumulate_confusion_matrices
from ..core.dataset import MutsDataset
from ..core.error import IncompatibleValuesError, OutOfBoundsError
from ..core.logs import logger
from ..core.seq.xna import DNA
from ..core.seq.region import unite
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write
from ..filter.dataset import FilterMutsDataset, load_filter_dataset
from ..filter.io import FilterBatchIO
from ..filter.main import load_regions, set_mut_gap_params
from ..filter.report import FilterReport
from .report import FilterScanReport

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

PAIRS_CSV = "pairs.csv"
DOMAINS_CSV = "domains.csv"
PAIRS_DOMAINS_HTML = "pairs_and_domains.html"
CONFUSION_MATRIX_CSV = "confusion-matrix.csv"


def _get_batch_read_lengths(batch_num: int, dataset: MutsDataset):
    batch = dataset.get_batch(batch_num)
    return batch.read_lengths


def _ends_arrays_to_tuples(end5s: Iterable, end3s: Iterable):
    # The end5s and end3s arrays must be mapped to Python integers,
    # otherwise they will cause an error during JSON serialization.
    return list(zip(map(int, end5s), map(int, end3s), strict=True))


def _tuples_to_ends_arrays(pairs: list[tuple[int, int]]):
    import numpy as np

    if not pairs:
        return np.array([], dtype=int), np.array([], dtype=int)
    end5s, end3s = map(np.array, zip(*pairs))
    return end5s, end3s


def _calc_midpoints_distances(end5s: np.ndarray, end3s: np.ndarray):
    import numpy as np

    assert end5s.shape == end3s.shape
    assert np.all(end5s <= end3s)
    midpoints = (end5s + end3s) / 2
    distances = end3s - end5s
    return midpoints, distances


def _calc_tiles(
    total_end5: int, total_end3: int, tile_length: int, tile_min_overlap: float
):
    import numpy as np

    if not isinstance(total_end5, int):
        raise TypeError(total_end5)
    if not isinstance(total_end3, int):
        raise TypeError(total_end3)
    if not 1 <= total_end5 <= total_end3:
        raise IncompatibleValuesError(
            "Must have 1 ≤ total_end5 ≤ total_end3, "
            f"but got total_end5={total_end5} "
            f"and total_end3={total_end3}"
        )
    total_length = total_end3 - total_end5 + 1
    assert total_length >= 1
    if not isinstance(tile_length, int):
        raise TypeError(tile_length)
    if tile_length < 1:
        raise OutOfBoundsError(f"tile_length must be ≥ 1, but got {tile_length}")
    if tile_length > total_length:
        logger.warning(
            f"tile_length ({tile_length}) is greater than "
            f"total length of region ({total_length}): "
            f"using tile_length of {total_length}"
        )
        return [(total_end5, total_end3)]
    assert 1 <= tile_length <= total_length <= total_end3
    if not isinstance(tile_min_overlap, float):
        raise TypeError(tile_min_overlap)
    if not 0.0 < tile_min_overlap < 1.0:
        raise OutOfBoundsError(
            f"tile_min_overlap must be > 0 and < 1, but got {tile_min_overlap}"
        )
    max_step_size = int(tile_length * (1.0 - tile_min_overlap))
    assert 0 <= max_step_size < tile_length
    if max_step_size == 0:
        raise IncompatibleValuesError(
            f"Cannot have tile_length={tile_length} "
            f"with tile_min_overlap={tile_min_overlap}"
        )
    num_regions = 1 + ceil((total_length - tile_length) / max_step_size)
    region_end5s = np.asarray(
        np.round(np.linspace(total_end5, total_end3 - tile_length + 1, num_regions)),
        dtype=int,
    )
    region_end3s = region_end5s + (tile_length - 1)
    return _ends_arrays_to_tuples(region_end5s, region_end3s)


def _calc_tile_coords(
    dataset, total_region, tile_length: int, tile_min_overlap: float, num_cpus: int
):
    """Calculate the tiled coordinates for one IDmut dataset."""
    import numpy as np

    ref = dataset.ref
    with logger.debug.single_context("calculating tiles for reference {!r}", ref):
        if tile_length > 0:
            # Use a prespecified region length.
            ref_tile_length = tile_length
        else:
            # Set the region length to twice the median read length.
            logger.debug("Began calculating optimal tile length")
            batches_read_lengths = dispatch(
                _get_batch_read_lengths,
                num_cpus=num_cpus,
                pass_num_cpus=False,
                as_list=True,
                ordered=False,
                raise_on_error=True,
                args=as_list_of_tuples(dataset.batch_nums),
                kwargs=dict(dataset=dataset),
            )
            if sum(a.size for a in batches_read_lengths) == 0:
                raise ValueError(f"{dataset} has 0 reads")
            read_lengths = np.concatenate(batches_read_lengths, axis=0)
            median_read_length = np.median(read_lengths)
            if median_read_length < 1:
                raise ValueError(
                    f"The median read length must be ≥ 1, but got {median_read_length}"
                )
            logger.trace("The median read length is {}", median_read_length)
            ref_tile_length = round(2 * median_read_length)
            logger.debug("Ended calculating optimal tile length: {}", ref_tile_length)
        tiles = _calc_tiles(
            total_region.end5, total_region.end3, ref_tile_length, tile_min_overlap
        )
    return [(ref, end5, end3) for end5, end3 in tiles]


def _find_correlated_pairs(dataset: FilterMutsDataset, pair_fdr: float, num_cpus: int):
    # Calculate the confusion matrix for the region.
    import pandas as pd

    n, a, b, ab = accumulate_confusion_matrices(
        dataset.get_batch,
        dataset.num_batches,
        dataset.pattern,
        dataset.region.unmasked,
        None,
        min_gap=dataset.min_mut_gap,
        num_cpus=num_cpus,
    )
    # Determine which pairs of positions correlate significantly at a
    # false discovery rate of pair_fdr.
    pvals = calc_confusion_pvals(n, a, b, ab)
    pvals_bh_adjusted = calc_bh_adjusted_pvals(pvals)
    is_significant = pvals_bh_adjusted <= pair_fdr
    # Save the confusion matrix as a CSV file.
    csv_file = dataset.report_file.with_name(CONFUSION_MATRIX_CSV)
    phi = calc_confusion_phi(n, a, b, ab)
    confusion_matrix = pd.DataFrame.from_dict(
        {
            "Neither Mutated": n - (a + b - ab),
            "Only A Mutated": a - ab,
            "Only B Mutated": b - ab,
            "Both Mutated": ab,
            "Phi Correlation": phi,
            "Raw P-Value": pvals,
            "Adjusted P-Value": pvals_bh_adjusted,
            f"Significant at FDR={pair_fdr}": is_significant,
        }
    )
    confusion_matrix.to_csv(csv_file)
    # Return the significantly correlated pairs.
    pairs = n.index
    pairs_significant = pairs[is_significant].to_list()
    logger.trace(
        "{} has {} pairs with significant correlations out of {} total pairs",
        dataset,
        len(pairs_significant),
        len(pairs),
    )
    return pairs_significant


def _aggregate_pairs(pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Aggregate pairs that overlap."""
    with logger.debug.single_context(
        "Generating initial domains by aggregating {} correlated pair(s)", len(pairs)
    ):
        if not pairs:
            return list()
        assert len(pairs[0]) == 2
        aggregates = [pairs[0]]
        for pair in pairs[1:]:
            pos5, pos3 = pair
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
                aggregates.append(pair)
            logger.trace(
                "After aggregating pair {}, the last domain is {}", pair, aggregates[-1]
            )
        logger.debug(
            "Generated {} initial domain(s) by aggregating {} correlated pair(s)",
            len(aggregates),
            len(pairs),
        )
        return aggregates


def _select_pairs(pairs: list[tuple[int, int]], end5: int, end3: int):
    assert 1 <= end5 <= end3
    return [(pos5, pos3) for pos5, pos3 in pairs if end5 <= pos5 and pos3 <= end3]


def _compute_keep_dists_null(
    valid_pairs: list[tuple[int, int]],
    total_end5: int,
    total_end3: int,
    unmasked: np.ndarray,
):
    """Expected starts/ends count per position under a null that keeps
    each pair's own length fixed and randomizes only its position.

    A pair may be placed at start ``x`` only if both its 5' end ``x`` and
    its 3' end ``x + L`` fall on **unmasked** positions, since masked
    positions can never be a pair's endpoint.  For a pair of length L,
    the number of valid locations is therefore
    ``num_loc = sum_x u[x] u[x+L]`` over the unmasked indicator ``u``;
    each valid location contributes probability ``1 / num_loc`` of the
    pair landing there.  Summing over all valid pairs (grouped by length,
    for efficiency) gives the Poisson rate at each position.  When no
    positions are masked this reduces to ``num_loc = num_pos_total - L``.

    Returns ``(null_expect_pos5s, null_expect_pos3s)``, each a numpy
    array indexed by position (index 0 == total_end5).
    """
    import numpy as np

    num_pos_total = total_end3 - total_end5 + 1
    # Indicator of unmasked positions over [total_end5, total_end3].
    u = np.zeros(num_pos_total)
    idx = np.asarray(unmasked, dtype=int) - total_end5
    idx = idx[(idx >= 0) & (idx < num_pos_total)]
    u[idx] = 1.0
    # Count the pairs with each length (defined as pos3 - pos5).
    length_counts: defaultdict[int, int] = defaultdict(int)
    for pos5, pos3 in valid_pairs:
        length_counts[pos3 - pos5] += 1
    # Count the pairs expected to have their 5'/3' ends at each position
    # under the null model where pairs are moved randomly while keeping
    # their lengths unchanged, restricted to unmasked positions.
    null_expect_pos5s = np.zeros(num_pos_total)
    null_expect_pos3s = np.zeros(num_pos_total)
    for length, count in length_counts.items():
        if length >= num_pos_total:
            continue
        # A pair of this length can start at x only if x and x + length
        # are both unmasked.
        valid_starts = u[: num_pos_total - length] * u[length:]
        num_loc = valid_starts.sum()
        if num_loc <= 0.0:
            continue
        prob_per_loc = count / num_loc
        null_expect_pos5s[: num_pos_total - length] += prob_per_loc * valid_starts
        null_expect_pos3s[length:] += prob_per_loc * valid_starts
    return null_expect_pos5s, null_expect_pos3s


def _compute_endpoint_significance(
    valid_pairs: list[tuple[int, int]],
    total_end5: int,
    total_end3: int,
    unmasked: np.ndarray,
):
    """Per-position BH-adjusted p-values for endpoint over-representation,
    under the length-preserving ("keep_dists") null restricted to the
    unmasked positions.

    Each position p is tested on its own count of pair 5' (or 3') ends
    against the null expectation at p: a valid pair independently places
    its start (or end) at p with probability ``null_expect[p] / r_valid``,
    so the count follows ``Binomial(r_valid, null_expect[p] / r_valid)``.

    Returns ``(positions, padj_pos5s, padj_pos3s)``, or ``None`` if
    there are no valid pairs.
    """
    import numpy as np
    from scipy.stats import binom

    if not valid_pairs:
        return None

    num_pos = total_end3 - total_end5 + 1
    r_valid = len(valid_pairs)
    positions = np.arange(total_end5, total_end3 + 1)
    null_expect_pos5s, null_expect_pos3s = _compute_keep_dists_null(
        valid_pairs, total_end5, total_end3, unmasked
    )

    pos5s = np.zeros(num_pos, dtype=int)
    pos3s = np.zeros(num_pos, dtype=int)
    for pos5, pos3 in valid_pairs:
        pos5s[pos5 - total_end5] += 1
        pos3s[pos3 - total_end5] += 1

    pvals_pos5s = binom.sf(pos5s - 1, r_valid, null_expect_pos5s / r_valid)
    pvals_pos3s = binom.sf(pos3s - 1, r_valid, null_expect_pos3s / r_valid)
    # Starts and ends are independent hypotheses; each gets its own
    # BH correction rather than one correction over both combined.
    padj_pos5s = calc_bh_adjusted_pvals(pvals_pos5s)
    padj_pos3s = calc_bh_adjusted_pvals(pvals_pos3s)
    return positions, padj_pos5s, padj_pos3s


def _filter_by_l1_distance(
    pairs: list[tuple[int, int]], percentile: float = 95.0, min_nearby_pairs: int = 1
):
    """Keep only pairs that have at least ``min_nearby_pairs`` other
    surviving pairs within the ``percentile``-th percentile of the L1
    (Manhattan) nearest-neighbor distances among ``pairs``.

    Uses the L1-norm ``|Δpos5| + |Δpos3|``, so pairs that are close
    in both coordinates (e.g. a helix with adjacent registers) or
    share one coordinate exactly (e.g. identical start or end) all
    contribute naturally to the minimum distance, with no requirement
    to group by a shared coordinate.

    Setting ``min_nearby_pairs`` above 1 removes coincidental small
    clusters of noise pairs ("buddy noise") at the cost of potentially
    clipping pairs at domain edges.

    Because the endpoint-peak stage already filters out pairs at
    non-significant hub positions, stage-1 survivors are enriched for
    real pairs and depleted for false positives relative to the raw
    input pairs (which start at ~5% FDR).  Calibrating the threshold
    from this already-filtered population therefore reflects a lower
    effective noise rate, making the percentile a self-calibrated and
    dataset-adaptive threshold requiring no manual tuning.

    The filter is applied iteratively until no more pairs are dropped.
    The threshold is computed once from the initial population and held
    fixed; subsequent rounds recompute only the NN-distances within the
    shrinking survivor set, so pairs that lose their dense neighbours
    across rounds are eventually caught.  Convergence is guaranteed
    because ``surviving ⊆ pairs`` after every round.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    if len(pairs) <= min_nearby_pairs:
        return []

    def get_nn_dist(ps: list[tuple[int, int]], n: int):
        # index 0 is the point itself (distance 0); k=(n + 1) selects the
        # n-th nearest neighbor; p=1 gives L1 distance.
        assert len(ps) > 0
        points = np.asarray(ps, dtype=float)
        return cKDTree(points).query(points, k=(n + 1), p=1)[0][:, -1]

    # Calculate the distances to the min_nearby_pairs-th nearest
    # neighbor (n=min_nearby_pairs).
    nn_dist = get_nn_dist(pairs, min_nearby_pairs)
    # Compute threshold once from the initial population.
    threshold = float(np.percentile(nn_dist, percentile))
    # Iterate with the fixed threshold until convergence.
    surviving = pairs
    while True:
        next_surviving = [pair for pair, d in zip(surviving, nn_dist) if d <= threshold]
        if next_surviving == surviving:
            return surviving
        surviving = next_surviving
        if len(surviving) <= min_nearby_pairs:
            return []
        # Calculate the distances to the min_nearby_pairs-th nearest
        # neighbor (n=min_nearby_pairs).
        nn_dist = get_nn_dist(surviving, min_nearby_pairs)


def _calc_domains_from_pairs(
    pairs: list[tuple[int, int]],
    pair_fdr: float,
    min_mut_gap: int,
    min_pairs: int,
    total_end5: int,
    total_end3: int,
    pair_distance_percentile: float,
    min_nearby_pairs: int,
    unmasked: np.ndarray,
) -> list[tuple[int, int]]:
    """Detect domains by finding significant concentrations of pair
    endpoints, rather than testing pairs crossing each interior
    boundary.

    Stage 1 (endpoint peaks): retain a pair if at least one of its
    5' or 3' position is significantly over-represented (``pair_fdr``)
    under the length-preserving null, whose placement space is
    restricted to the ``unmasked`` positions (masked positions can
    never be a pair endpoint).

    Stage 2 (L1 isolation filter): drop any stage-1 survivor whose
    nearest other survivor (L1 distance in (pos5, pos3) space) exceeds
    the ``pair_distance_percentile``-th percentile of the observed L1
    nearest-neighbor distances among stage-1 survivors.  This adapts
    automatically to each dataset.

    Retained pairs are merged into domains by interval overlap.
    """
    with logger.debug.single_context(
        "Identifying domains from {} correlated pair(s)", len(pairs)
    ):
        if not pairs:
            return []
        valid_pairs = [(p5, p3) for p5, p3 in pairs if p3 - p5 > min_mut_gap]
        if not valid_pairs:
            return []

        significance = _compute_endpoint_significance(
            valid_pairs, total_end5, total_end3, unmasked
        )
        if significance is None:
            return []
        positions, padj_starts, padj_ends = significance

        sig_starts = set(positions[padj_starts <= pair_fdr].tolist())
        sig_ends = set(positions[padj_ends <= pair_fdr].tolist())
        surviving = sorted(
            (p5, p3) for p5, p3 in valid_pairs if p5 in sig_starts or p3 in sig_ends
        )
        surviving = _filter_by_l1_distance(
            surviving,
            percentile=pair_distance_percentile,
            min_nearby_pairs=min_nearby_pairs,
        )
        if not surviving:
            return []

        domains = _aggregate_pairs(surviving)
        domains = [
            (e5, e3)
            for e5, e3 in domains
            if len(_select_pairs(surviving, e5, e3)) >= min_pairs
        ]
        logger.debug("Calculated domains: {}", domains)
        return domains


def _graph_pairs_and_domains(
    pairs: list[tuple[int, int]],
    domains: list[tuple[int, int]],
    end5: int,
    end3: int,
    html_file: str | Path,
):
    import numpy as np
    from plotly import graph_objects as go
    from plotly.subplots import make_subplots

    # Create a subplot with two rows: top for pair_fraction,
    # bottom for correlated pairs
    fig = make_subplots(rows=1, cols=1)
    # Graph the domains as triangles.
    end5s, end3s = _tuples_to_ends_arrays(domains)
    domains_midpoints, domains_distances = _calc_midpoints_distances(end5s, end3s)
    for (a, b), x, y in zip(domains, domains_midpoints, domains_distances, strict=True):
        fig.add_trace(
            go.Scatter(
                x=[a, b, x, a],
                y=[0, 0, y, 0],
                mode="none",
                fill="toself",
                fillcolor="rgba(230,159,0,0.5)",
                showlegend=False,
                name=None,
                hoverinfo="text",
                text=f"Domain {a, b}",
                hovertemplate="%{text}<extra></extra>",
            )
        )
    # Plot the correlated pairs as points.
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    pairs_midpoints, pairs_distances = _calc_midpoints_distances(pos5s, pos3s)
    fig.add_trace(
        go.Scatter(
            x=pairs_midpoints,
            y=pairs_distances,
            mode="markers",
            showlegend=False,
            marker=dict(color="#D55E00"),
            name=None,
            hoverinfo="text",
            text=[f"Pair {pair}" for pair in pairs],
            hovertemplate="%{text}<extra></extra>",
        )
    )
    # Finish the layout.
    assert end5 <= end3
    x_range = [end5 - 0.5, end3 + 0.5]
    fig.update_xaxes(title_text="Midpoint Position", showgrid=True, range=x_range)
    if domains_distances.size > 0 or pairs_distances.size > 0:
        y_range = [
            0.0,
            1.05 * np.max(np.concatenate([domains_distances, pairs_distances])),
        ]
    else:
        y_range = [0.0, 1.0]
    fig.update_yaxes(title_text="Length of Span", showgrid=True, range=y_range)
    # Save the figure.
    fig.write_html(html_file)


def _insert_domains_into_gaps(
    domains: list[tuple[int, int]], global_end5: int, global_end3: int
):
    """Turn every gap between domains into a new domain."""
    assert 1 <= global_end5 <= global_end3
    new_domains = list()
    prev_pos3 = global_end5 - 1
    # domains is assumed to be sorted.
    for end5, end3 in domains:
        assert global_end5 <= end5 <= end3 <= global_end3
        # Here, no two regions are allowed to overlap.
        assert end5 > prev_pos3
        if end5 > prev_pos3 + 1:
            # There is a gap between this region and the previous.
            # Create a new region to fill the gap.
            new_domains.append(((prev_pos3 + 1), (end5 - 1)))
        new_domains.append((end5, end3))
        prev_pos3 = end3
    if prev_pos3 < global_end3:
        # There is a gap between the last region and total_end3.
        # Create a new region to fill the gap.
        new_domains.append(((prev_pos3 + 1), global_end3))
    return new_domains


def _expand_domains_into_gaps(
    domains: list[tuple[int, int]], global_end5: int, global_end3: int
):
    """Expand every domain to fill gaps on either side."""
    assert 1 <= global_end5 <= global_end3
    if not domains:
        return list()
    # Make mutable end5s/end3s arrays.
    end5s, end3s = map(list, zip(*domains))
    assert 1 <= len(domains) == len(end5s) == len(end3s)
    # Expand to cover the global ends.
    end5s[0] = global_end5
    end3s[-1] = global_end3
    # Split each internal gap: left gets floor, right gets ceil
    for i in range(len(domains) - 1):
        left_end3 = end3s[i]
        right_end5 = end5s[i + 1]
        gap = right_end5 - left_end3 - 1
        if gap > 0:
            expand_left = gap // 2
            expand_right = gap - expand_left
            end3s[i] = left_end3 + expand_left
            end5s[i + 1] = right_end5 - expand_right
    return _ends_arrays_to_tuples(end5s, end3s)


def _filter_domains_length(
    domains: list[tuple[int, int]],
    min_length: int | float = 1,
    max_length: int | float = inf,
):
    """Remove domains that are too short or too long."""
    return [
        (end5, end3)
        for end5, end3 in domains
        if min_length <= (end3 - end5 + 1) <= max_length
    ]


def _write_pairs_to_csv(pairs: list[tuple[int, int]], csv_file: str | Path):
    """Write the pairs to a CSV file."""
    import pandas as pd

    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    df = pd.DataFrame.from_dict(
        {POSITION_A: pos5s, POSITION_B: pos3s}, orient="columns"
    )
    df.to_csv(csv_file, index=False)


def _calc_cluster_domains(
    filter_dirs: list[Path],
    report_dir: Path,
    num_cpus: int,
    pair_fdr: float,
    probe: str,
    min_mut_gap: int | None,
    min_pairs: int,
    pair_distance_percentile: float,
    min_nearby_pairs: int,
    min_length: int,
    gap_mode: str,
):
    """Calculate the cluster regions for all tiles of one reference."""
    # Each dataset corresponds to one tile.
    datasets = list(load_filter_dataset.iterate(filter_dirs))
    if not datasets:
        return list(), 0
    # Find the common ref, refseq, top, and branches (must be
    # identical among all datasets).
    ref = datasets[0].ref
    refseq = datasets[0].refseq
    branches = datasets[0].branches
    for dataset in datasets[1:]:
        if dataset.ref != ref:
            raise ValueError(
                f"Expected all tile datasets to have reference "
                f"{repr(ref)}, but got {repr(dataset.ref)}"
            )
        if dataset.refseq != refseq:
            raise ValueError(
                f"Expected all tile datasets to have the same "
                f"reference sequence as {repr(ref)}, but got a "
                f"different sequence for {repr(dataset.ref)}"
            )
        if dataset.branches != branches:
            raise ValueError(
                f"Expected all tile datasets to have branches "
                f"{branches}, but got {dataset.branches}"
            )
    # The region is the union of all tiles' regions.
    region = unite([dataset.region for dataset in datasets], refseq=refseq)
    # Find pairs of bases that correlate significantly in any region.
    # Use sorted(set()) to ensure pairs are unique and sorted.
    pairs = sorted(
        set(
            chain(
                *dispatch(
                    _find_correlated_pairs,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=False,
                    ordered=False,
                    raise_on_error=True,
                    args=as_list_of_tuples(datasets),
                    kwargs=dict(pair_fdr=pair_fdr),
                )
            )
        )
    )
    # Find domains of correlated pairs.
    min_mut_gap, _ = set_mut_gap_params(probe, min_mut_gap)
    domains = _calc_domains_from_pairs(
        pairs,
        pair_fdr,
        min_mut_gap,
        min_pairs,
        total_end5=region.end5,
        total_end3=region.end3,
        pair_distance_percentile=pair_distance_percentile,
        min_nearby_pairs=min_nearby_pairs,
        unmasked=region.unmasked_int,
    )
    n_domains_before_filter = len(domains)
    # Determine what to do with gaps between regions.
    if gap_mode == GAP_MODE_INSERT:
        domains = _filter_domains_length(domains, min_length=min_length)
        domains = _insert_domains_into_gaps(domains, region.end5, region.end3)
    elif gap_mode == GAP_MODE_EXPAND:
        domains = _expand_domains_into_gaps(domains, region.end5, region.end3)
    elif gap_mode == GAP_MODE_OMIT:
        domains = _filter_domains_length(domains, min_length=min_length)
    else:
        raise ValueError(gap_mode)
    logger.debug("Ended adjusting domains: {}", domains)
    if n_domains_before_filter > 0 and not domains:
        logger.warning(
            "All {} domain(s) were removed by the minimum-length filter (min={})",
            n_domains_before_filter,
            min_length,
        )
    # Write the pairs and domains to CSV files.
    report_dir.mkdir(parents=True, exist_ok=True)
    _write_pairs_to_csv(pairs, report_dir.joinpath(PAIRS_CSV))
    _write_pairs_to_csv(domains, report_dir.joinpath(DOMAINS_CSV))
    # Graph the correlated pairs and domains.
    html_file = report_dir.joinpath(PAIRS_DOMAINS_HTML)
    try:
        _graph_pairs_and_domains(pairs, domains, region.end5, region.end3, html_file)
    except Exception as error:
        logger.error(error)
    return [(ref, end5, end3) for end5, end3 in domains], len(pairs)


def filterscan(
    idmut_report_file: Path,
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
    pair_fdr: float,
    min_pairs: int,
    pair_distance_percentile: float,
    min_nearby_pairs: int,
    min_cluster_length: int,
    gap_mode: str,
    # Filter options
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    primer_gap: int,
    regions_file: Path | None,
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
    """Scan one IDmut dataset for domains of correlated base pairs.

    Run the filter step over overlapping tiles spanning the RNA, detect
    domains of correlated base pairs, filter the reads over each domain,
    and write a FilterScanReport recording the domain coordinates. The
    tiles are then deleted, leaving the domain filter results for
    clusterscan to cluster.
    """
    # Load region info cheaply (reads JSON only, no batch I/O).
    datasets, total_regions = load_regions(
        [idmut_report_file], region_coords, region_primers, primer_gap, regions_file
    )
    refs = total_regions.refs
    assert len(refs) == 1
    ref = refs[0]
    ref_total_regions = total_regions.dict[ref]
    assert len(ref_total_regions) == 1
    total_region = ref_total_regions[0]
    assert list(datasets.keys()) == [ref]
    assert len(datasets[ref]) == 1
    idmut_dataset = datasets[ref][0]
    # Check if the FilterScanReport already exists.
    report_branches = path.add_branch(
        path.FILTERSCAN_STEP, branch, idmut_dataset.branches
    )
    report_file = FilterScanReport.build_path(
        {
            path.TOP: idmut_dataset.top,
            path.SAMPLE: idmut_dataset.sample,
            path.BRANCHES: report_branches,
            path.REF: idmut_dataset.ref,
            path.REG: total_region.name,
        }
    )
    if need_write(report_file, force):
        began = datetime.now()
        # Compute tile coordinates (potentially expensive when tile_length <= 0).
        tile_coords = _calc_tile_coords(
            idmut_dataset, total_region, tile_length, tile_min_overlap, num_cpus
        )
        tiled_dirs = filter_mod.run(
            input_path=[idmut_report_file],
            branch=branch,
            tmp_pfx=tmp_pfx,
            keep_tmp=keep_tmp,
            region_coords=tile_coords,
            region_primers=(),
            primer_gap=0,
            regions_file=None,
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
            filter_pos_table=False,
            filter_read_table=False,
            self_contained=self_contained,
            brotli_level=brotli_level,
            num_cpus=num_cpus,
            force=force,
        )
        # Find regions spanned by correlated base pairs.
        report_dir = report_file.parent
        domain_coords, n_signif_pairs = _calc_cluster_domains(
            tiled_dirs,
            report_dir=report_dir,
            pair_fdr=pair_fdr,
            probe=probe,
            min_mut_gap=min_mut_gap,
            min_pairs=min_pairs,
            pair_distance_percentile=pair_distance_percentile,
            min_nearby_pairs=min_nearby_pairs,
            min_length=min_cluster_length,
            gap_mode=gap_mode,
            num_cpus=num_cpus,
        )
        if domain_coords:
            # Filter the reads over each domain so that clusterscan can
            # cluster them without re-running the filter step.
            filter_mod.run(
                input_path=[idmut_report_file],
                branch=branch,
                tmp_pfx=tmp_pfx,
                keep_tmp=keep_tmp,
                region_coords=tuple(domain_coords),
                region_primers=(),
                primer_gap=0,
                regions_file=None,
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
                brotli_level=brotli_level,
                num_cpus=num_cpus,
                force=force,
            )
        else:
            logger.warning(
                "No domains of correlated pairs found in {}", idmut_report_file
            )
        FilterScanReport(
            sample=idmut_dataset.sample,
            ref=idmut_dataset.ref,
            reg=total_region.name,
            branches=report_branches,
            # Domain-detection parameters.
            tile_length=tile_length,
            tile_min_overlap=tile_min_overlap,
            erase_tiles=erase_tiles,
            pair_fdr=pair_fdr,
            min_pairs=min_pairs,
            pair_distance_percentile=pair_distance_percentile,
            min_nearby_pairs=min_nearby_pairs,
            min_cluster_length=min_cluster_length,
            gap_mode=gap_mode,
            # Results (store coordinates without the reference, which is
            # already recorded in the report).
            tile_coords=[(end5, end3) for _, end5, end3 in tile_coords],
            n_signif_pairs=n_signif_pairs,
            n_domains=len(domain_coords),
            domain_coords=[(end5, end3) for _, end5, end3 in domain_coords],
            began=began,
            ended=datetime.now(),
        ).save(idmut_dataset.top, force=force)
        if erase_tiles:
            # Delete the filter reports and batches of the tiles.
            with logger.debug.single_context("Erasing tiles"):
                for file in path.find_files_chain(
                    tiled_dirs, FilterReport.get_path_seg_types()
                ):
                    file.unlink(missing_ok=True)
                    logger.trace("Erased {}", file)
                for file in path.find_files_chain(
                    tiled_dirs, FilterBatchIO.get_path_seg_types()
                ):
                    file.unlink(missing_ok=True)
                    logger.trace("Erased {}", file)
    return report_file.parent
