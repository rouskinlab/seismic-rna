from __future__ import annotations
from datetime import datetime
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
    calc_confusion_chi_square,
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
from ..filter.main import load_regions
from ..filter.report import FilterReport
from .report import FilterScanReport

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

PAIRS_CSV = "pairs.csv"
DOMAINS_CSV = "domains.csv"
PAIRS_DOMAINS_HTML = "pairs_and_domains.html"
CONFUSION_MATRIX_CSV = "confusion-matrix.csv"

N_COL = "N"
PHI_COL = "Phi"
CHI_SQUARE_COL = "Chi-Square"
SIGNIFICANT_COL = "Significant"


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


def _find_correlated_pairs(
    dataset: FilterMutsDataset, pair_fdr: float, band_width: int, num_cpus: int
):
    """Calculate the confusion matrix for the region, restricted to
    pairs of positions no more than ``band_width`` apart (0 means no
    additional restriction beyond the tile itself, which already
    bounds every pair's separation to less than the tile length).

    Returns a table indexed by every in-band pair of positions, with
    columns for its coverage (N), phi correlation, chi-square score
    (``N * phi ** 2``), and whether it is significant at the given
    false discovery rate (by the same hypergeometric + BH test as
    before; the chi-square score is used only downstream, to find
    domains, not to flag individual pairs).
    """
    import pandas as pd

    n, a, b, ab = accumulate_confusion_matrices(
        dataset.get_batch,
        dataset.num_batches,
        dataset.pattern,
        dataset.region.unmasked,
        None,
        min_gap=dataset.min_mut_gap,
        max_gap=(band_width if band_width > 0 else None),
        num_cpus=num_cpus,
    )
    # Determine which pairs of positions correlate significantly at a
    # false discovery rate of pair_fdr.
    pvals = calc_confusion_pvals(n, a, b, ab)
    pvals_bh_adjusted = calc_bh_adjusted_pvals(pvals)
    is_significant = pvals_bh_adjusted <= pair_fdr
    phi = calc_confusion_phi(n, a, b, ab)
    chi_square = calc_confusion_chi_square(n, phi)
    # Save the confusion matrix as a CSV file.
    csv_file = dataset.report_file.with_name(CONFUSION_MATRIX_CSV)
    confusion_matrix = pd.DataFrame.from_dict(
        {
            "Neither Mutated": n - (a + b - ab),
            "Only A Mutated": a - ab,
            "Only B Mutated": b - ab,
            "Both Mutated": ab,
            "Phi Correlation": phi,
            CHI_SQUARE_COL: chi_square,
            "Raw P-Value": pvals,
            "Adjusted P-Value": pvals_bh_adjusted,
            f"Significant at FDR={pair_fdr}": is_significant,
        }
    )
    confusion_matrix.to_csv(csv_file)
    table = pd.DataFrame(
        {
            N_COL: n,
            PHI_COL: phi,
            CHI_SQUARE_COL: chi_square,
            SIGNIFICANT_COL: is_significant,
        }
    )
    logger.trace(
        "{} has {} pair(s) within the band, of which {} are significant",
        dataset,
        len(table.index),
        int(is_significant.sum()),
    )
    return table


def _build_banded_table(
    per_tile_tables: list[pd.DataFrame], band_width: int
) -> pd.DataFrame:
    """Merge per-tile pair tables into one table with one row per
    unique in-band pair of positions, keeping the tile observation
    with the greatest coverage (N) for any pair seen in more than one
    overlapping tile (the highest-power estimate, and avoids double-
    counting the same reads under two different tilings).

    ``band_width`` caps the separation ``j - i`` of every row (0 means
    no additional cap beyond what the tiles already impose, since
    every pair's separation is already less than the tile length).
    """
    import pandas as pd

    if not per_tile_tables:
        return pd.DataFrame(
            {col: pd.Series(dtype=float) for col in (N_COL, PHI_COL, CHI_SQUARE_COL)}
            | {SIGNIFICANT_COL: pd.Series(dtype=bool)},
            index=pd.MultiIndex.from_arrays([[], []], names=[POSITION_A, POSITION_B]),
        )
    combined = pd.concat(per_tile_tables, axis=0)
    combined = combined[combined[CHI_SQUARE_COL].notna()]
    # For pairs seen in more than one tile, keep only the observation
    # with the greatest coverage (N).
    combined = combined.sort_values(N_COL, ascending=False)
    combined = combined[~combined.index.duplicated(keep="first")]
    if band_width > 0:
        gaps = combined.index.get_level_values(
            POSITION_B
        ) - combined.index.get_level_values(POSITION_A)
        combined = combined[gaps <= band_width]
    return combined.sort_index()


def _pair_band_row_cumsum(
    table: pd.DataFrame, total_end5: int, total_end3: int
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Build banded row-cumulative sums of the observable-pair and
    significant-pair grids over positions [total_end5, total_end3],
    0-indexed, so that the number of observable/significant pairs
    whose triangle lies within any interval [s, e] (both endpoints in
    [s, e]) can be read off via ``_triangle_sum_banded`` -- see
    ``_calc_domains_by_dp_segmentation`` for how this drives the exact
    dynamic-program segmentation.

    Every pair (a, b) that reaches this function has already passed
    through ``_build_banded_table``'s ``band_width`` cap (and, before
    that, ``_calc_domains_by_dp_segmentation``'s own coverage filter),
    so ``b - a`` is bounded by some ``max_gap`` read directly off the
    table's own data. Representing the grid banded, as an
    ``(n_positions, max_gap + 1)`` array of gaps rather than a dense
    ``(n_positions, n_positions)`` array of positions, cuts memory
    from O(L^2) to O(L * max_gap) -- for a reference tens of thousands
    of positions long with a band a few hundred wide, the difference
    between tens of gigabytes and tens of megabytes.

    ``table`` is assumed already restricted to observable pairs (see
    ``_calc_domains_by_dp_segmentation``): a pair filtered out there
    (including every pair touching a fully masked position) never
    reaches this function, so it is scattered into neither returned
    array and cannot contribute to any downstream count, rate, or
    likelihood.

    Returns ``(obs_row_cum, sig_row_cum, n_positions, max_gap)``, where
    ``row_cum[a, d]`` (for ``d`` in ``0..max_gap+1``) is the number of
    observable/significant pairs ``(a, a+d')`` for ``d' = 0..d-1``
    (``row_cum[a, 0] == 0``; ``row_cum[a, max_gap+1]`` is position
    ``a``'s full row total within the band).
    """
    import numpy as np

    n_positions = total_end3 - total_end5 + 1
    pos_a = table.index.get_level_values(POSITION_A).to_numpy() - total_end5
    pos_b = table.index.get_level_values(POSITION_B).to_numpy() - total_end5
    sig = table[SIGNIFICANT_COL].to_numpy()
    gaps = pos_b - pos_a
    max_gap = int(gaps.max()) if len(gaps) else 0

    obs_grid = np.zeros((n_positions, max_gap + 1), dtype=np.int64)
    sig_grid = np.zeros((n_positions, max_gap + 1), dtype=np.int64)
    np.add.at(obs_grid, (pos_a, gaps), 1)
    np.add.at(sig_grid, (pos_a[sig], gaps[sig]), 1)

    obs_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
    sig_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
    obs_row_cum[:, 1:] = np.cumsum(obs_grid, axis=1)
    sig_row_cum[:, 1:] = np.cumsum(sig_grid, axis=1)
    return obs_row_cum, sig_row_cum, n_positions, max_gap


def _triangle_sum_banded(row_cum: np.ndarray, e: int, max_gap: int) -> np.ndarray:
    """Return an array of length e+1: the triangle sum n(s, e) (or
    k(s, e), depending on which row-cumulative array is passed) for
    every candidate block start s = 0..e (0-indexed), given a banded
    row-cumulative array from ``_pair_band_row_cumsum``.

    Splits each candidate start's contribution into two parts to
    avoid ever touching an O(L) x O(L) array: positions
    ``a <= e - max_gap`` contribute their *full* row (every pair they
    have is within [a, e] automatically, since gap <= max_gap <=
    e - a) -- a single O(1) range-sum via a precomputed 1-D prefix of
    row totals; positions within ``max_gap`` of ``e`` contribute a
    *partial* row up to distance ``e - a`` -- at most ``max_gap + 1``
    positions, gathered directly. Both parts are then broadcast across
    all s in one pass, so building the length-(e+1) output array
    (needed regardless, for the DP step at this e) is the only O(e)
    work; the rest is O(max_gap).
    """
    import numpy as np

    boundary = max(
        0, e - max_gap + 1
    )  # first position in the "near" (partial-row) window
    full_row = row_cum[:, max_gap + 1]
    cum_full = np.concatenate(([0], np.cumsum(full_row)))

    near_a = np.arange(boundary, e + 1)
    near_d = e - near_a
    near_vals = row_cum[near_a, near_d + 1]
    near_suffix = np.cumsum(near_vals[::-1])[
        ::-1
    ]  # near_suffix[i] = sum(near_vals[i:])
    near_total = near_suffix[0] if len(near_suffix) else 0

    result = np.zeros(e + 1, dtype=np.int64)
    if boundary > 0:
        s_far = np.arange(0, boundary)
        result[s_far] = (cum_full[boundary] - cum_full[s_far]) + near_total
    if len(near_a):
        result[boundary : e + 1] = near_suffix
    return result


def _block_gain(n, k, theta0: float):
    """Log-likelihood-ratio gain of modeling a block's ``n``
    observable pairs (``k`` of them significant) at their own rate
    ``theta_hat = k / n`` instead of the shared background rate
    ``theta0``: ``n * KL(theta_hat || theta0)``, i.e.
    ``k * log(theta_hat / theta0) + (n - k) * log((1 - theta_hat) /
    (1 - theta0))``. Zero unless ``theta_hat > theta0`` (a block only
    "wins" by being *denser* than background, never sparser) and zero
    wherever ``n == 0`` (an empty candidate carries no evidence)."""
    import numpy as np

    with np.errstate(divide="ignore", invalid="ignore"):
        theta_hat = np.where(n > 0, k / np.maximum(n, 1), 0.0)
        denser = (n > 0) & (theta_hat > theta0)
        term1 = np.where(
            k > 0,
            k * np.log(np.where(theta_hat > 0, theta_hat, 1.0) / max(theta0, 1e-300)),
            0.0,
        )
        n_minus_k = n - k
        complement = np.maximum(1.0 - theta_hat, 1e-300)
        term2 = np.where(
            n_minus_k > 0,
            n_minus_k * np.log(complement / max(1.0 - theta0, 1e-300)),
            0.0,
        )
    return np.where(denser, term1 + term2, 0.0)


def _dp_segment_blocks(
    obs_row_cum: np.ndarray,
    sig_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    theta0: float,
    penalty: float,
) -> list[tuple[int, int]]:
    """Find the partition of [0, n_positions) into background and
    domain blocks that maximizes the total block-model log-likelihood
    ratio (sum of ``_block_gain`` over the chosen blocks, minus
    ``penalty`` per block), via exact dynamic programming:
    ``dp[t] = max(dp[t-1], max_s dp[s] + gain(s, t-1) - penalty)``,
    where ``dp[t]`` is the best score using positions [0, t) and
    ``gain(s, t-1)`` scores a candidate block [s, t-1]. Every interval
    length at every position is considered (a block may be far longer
    than ``max_gap`` even though the band restricts any single pair's
    span, so this loop is still O(n_positions^2), unchanged from a
    dense grid -- only the O(n_positions * max_gap) memory improves),
    so no window scale is ever chosen; splitting a real domain costs
    the positive evidence of every pair that would have to move to the
    background (see ``_calc_domains_by_dp_segmentation`` for why this
    discourages over-splitting a domain bridged by many crossing
    pairs). Returns the chosen blocks as sorted, non-overlapping
    0-indexed (start, end) pairs, inclusive.
    """
    import numpy as np

    dp = np.zeros(n_positions + 1)
    backtrack = np.full(n_positions + 1, -1, dtype=np.int64)
    for t in range(1, n_positions + 1):
        e = t - 1
        s_vals = np.arange(t)
        n_se = _triangle_sum_banded(obs_row_cum, e, max_gap)
        k_se = _triangle_sum_banded(sig_row_cum, e, max_gap)
        gain = _block_gain(n_se, k_se, theta0)
        candidate_vals = np.where(n_se > 0, dp[s_vals] + gain - penalty, -np.inf)
        best_idx = int(np.argmax(candidate_vals))
        best_val = candidate_vals[best_idx]
        if best_val > dp[e]:
            dp[t] = best_val
            backtrack[t] = best_idx
        else:
            dp[t] = dp[e]
            backtrack[t] = -1
    blocks = []
    t = n_positions
    while t > 0:
        s = int(backtrack[t])
        if s == -1:
            t -= 1
        else:
            blocks.append((s, t - 1))
            t = s
    blocks.reverse()
    return blocks


def _estimate_background_rate(
    obs_row_cum: np.ndarray,
    sig_row_cum: np.ndarray,
    max_gap: int,
    blocks: list[tuple[int, int]],
    fallback: float,
) -> float:
    """Estimate the background significant-pair rate ``theta0`` from
    every observable pair *not* inside any of ``blocks`` -- i.e. every
    background and inter-block pair, which is inter-domain by
    construction and so uncontaminated by any called domain's own
    interior density (unlike a global average over the whole scanned
    region, which any real domain inflates). Falls back to
    ``fallback`` if there are no such pairs (e.g. one block covers the
    entire region)."""
    total_obs = int(obs_row_cum[:, max_gap + 1].sum())
    total_sig = int(sig_row_cum[:, max_gap + 1].sum())
    block_obs = sum(
        int(_triangle_sum_banded(obs_row_cum, e, max_gap)[s]) for s, e in blocks
    )
    block_sig = sum(
        int(_triangle_sum_banded(sig_row_cum, e, max_gap)[s]) for s, e in blocks
    )
    background_obs = total_obs - block_obs
    background_sig = total_sig - block_sig
    if background_obs <= 0:
        return fallback
    return background_sig / background_obs


def _merge_adjacent_blocks(
    obs_row_cum: np.ndarray,
    sig_row_cum: np.ndarray,
    max_gap: int,
    blocks: list[tuple[int, int]],
    theta0: float,
) -> list[tuple[int, int]]:
    """Greedily merge adjacent DP blocks, strongest bridge first,
    whenever the *bridge* between them -- everything in their union
    except their own two triangles, i.e. any gap between them plus
    every pair crossing from one into the other -- is itself denser
    than background (``_block_gain(bridge_n, bridge_k, theta0) > 0``),
    regardless of whether the DP's exact optimum happened to split
    them there.

    The DP in ``_dp_segment_blocks`` finds the globally optimal
    partition, so it never leaves a *beneficial* merge on the table:
    if joining two of its blocks scored higher, the DP would already
    have chosen that partition (it considers every candidate block,
    including their union, at every step). So a plain "does the union
    still look like a domain on its own" check is useless here -- by
    construction it can never fire when the DP declined a merge, and
    proved too permissive besides: it stayed true across a genuinely
    sparse bridge between two otherwise very dense sub-domains, purely
    because both sub-domains alone were dense enough to carry the
    union (see ``TestDpSegmentationAntiOverSplit``'s sparse-bridge
    case; that test still requires this pass to leave the split
    alone).

    What the DP's "sum of block gains" optimum discards, that a human
    eye doesn't, is a bridge that is real but *weaker* than either
    neighbor: a region with genuine internal density variation (one
    stretch at 3%, an adjacent one at 7%, both well above a background
    of 0.3%) is optimally split into purer pieces, because two own-
    rate fits beat one diluted shared rate -- even with zero gap
    between them. Testing the bridge alone, instead of the whole
    union, targets exactly that case without the sparse-bridge test's
    false positive: a truly silent gap or crossing region (bridge rate
    <= background) is left split, while a bridge that is merely
    weaker than its neighbors but still above background is merged.

    Merges are applied strongest-bridge-first (by ``_block_gain``,
    the same n-weighted evidence statistic used to decide whether a
    bridge qualifies at all -- not raw density, so a small, noisy
    high-ratio bridge can't jump ahead of a larger, more reliable
    one) rather than left to right, via a max-heap of candidate
    adjacent pairs. Every merge can change the bridges on either side
    of the absorbing block (its end position moves), so both
    neighboring pairs are recomputed and re-queued after each merge;
    stale heap entries (referring to a block that has since been
    absorbed, or whose end has since moved) are detected by an
    end-position version counter per block and discarded lazily
    rather than eagerly purged from the heap. This still lets 3+
    consecutive blocks with real bridges cascade into one, but the
    order in which ambiguous, weaker bridges get resolved no longer
    depends on which end of the region they happen to sit at.

    Uses the same ``obs_row_cum``/``sig_row_cum`` arrays as the DP, so
    a bridge's ``n``/``k`` are built only from pairs that already
    survived ``_calc_domains_by_dp_segmentation``'s coverage filter --
    a masked/excluded pair cannot contribute to a merge decision any
    more than it could to the DP itself.
    """
    import heapq

    if not blocks:
        return blocks

    n_blocks = len(blocks)
    starts = [s for s, _ in blocks]
    ends = [e for _, e in blocks]
    prev_idx: list[int | None] = [i - 1 if i > 0 else None for i in range(n_blocks)]
    next_idx: list[int | None] = [
        i + 1 if i < n_blocks - 1 else None for i in range(n_blocks)
    ]
    alive = [True] * n_blocks
    end_version = [0] * n_blocks

    def bridge_gain(left: int, right: int) -> float:
        prev_s, prev_e = starts[left], ends[left]
        s, e = starts[right], ends[right]
        union_obs_row = _triangle_sum_banded(obs_row_cum, e, max_gap)
        union_sig_row = _triangle_sum_banded(sig_row_cum, e, max_gap)
        union_n = int(union_obs_row[prev_s])
        union_k = int(union_sig_row[prev_s])
        next_n = int(union_obs_row[s])
        next_k = int(union_sig_row[s])
        prev_n = int(_triangle_sum_banded(obs_row_cum, prev_e, max_gap)[prev_s])
        prev_k = int(_triangle_sum_banded(sig_row_cum, prev_e, max_gap)[prev_s])
        bridge_n = union_n - prev_n - next_n
        bridge_k = union_k - prev_k - next_k
        return float(_block_gain(bridge_n, bridge_k, theta0))

    heap: list[tuple[float, int, int, int, int]] = []

    def push(left: int, right: int) -> None:
        gain = bridge_gain(left, right)
        if gain > 0:
            heapq.heappush(
                heap, (-gain, left, right, end_version[left], end_version[right])
            )

    for i in range(n_blocks - 1):
        push(i, i + 1)

    while heap:
        _, left, right, v_left, v_right = heapq.heappop(heap)
        if not (alive[left] and alive[right]):
            continue
        if next_idx[left] != right:
            continue
        if end_version[left] != v_left or end_version[right] != v_right:
            continue
        ends[left] = ends[right]
        end_version[left] += 1
        alive[right] = False
        next_idx[left] = next_idx[right]
        if next_idx[left] is not None:
            prev_idx[next_idx[left]] = left
            push(left, next_idx[left])
        if prev_idx[left] is not None:
            push(prev_idx[left], left)

    return [(starts[i], ends[i]) for i in range(n_blocks) if alive[i]]


N_BACKGROUND_RATE_ITERS = 3


def _calc_domains_by_dp_segmentation(
    table: pd.DataFrame,
    total_end5: int,
    total_end3: int,
    bic_multiplier: float,
    min_pair_coverage: int = 1000,
) -> list[tuple[int, int]]:
    """Call domains by exact dynamic-program segmentation of the
    significant-pair contact map into background/domain blocks (see
    ``_dp_segment_blocks``): partition [total_end5, total_end3] into
    consecutive intervals, each either background (rate ``theta0``) or
    a domain with its own rate, maximizing the total block-model
    log-likelihood ratio against an all-background null, penalized by
    ``bic_multiplier`` per block (Bayesian information criterion:
    ``0.5 * log(N observable pairs)`` per extra block, scaled by the
    multiplier).

    This replaces an earlier insulation-score boundary caller, which
    projected the 2-D contact map to a 1-D crossing-density track and
    so needed one fixed window scale that could not fit domains of
    very different sizes in the same reference (over-fragmenting
    some, under-resolving others). Scoring each candidate block
    directly by its own 2-D triangle density -- exactly what a person
    sees by eye in the pairs-vs-domains plot -- removes the window
    entirely: the DP considers every possible interval length at every
    position.

    Splitting a real domain into pieces is discouraged even before the
    per-block penalty: every pair crossing the split point is a
    significant pair that would have to be re-classified as
    background (rate ``theta0``, much lower than a domain's own rate),
    forfeiting its positive evidence. A split therefore only wins when
    the crossing pairs are genuinely sparse and/or the two halves have
    sharply different densities -- i.e. it really is two domains, not
    a single domain bridged by many crossing pairs.

    ``theta0`` is estimated from the pairs outside all called blocks
    (see ``_estimate_background_rate``), re-estimated for a few
    passes (EM-style) after each DP run so it stays uncontaminated by
    domain interiors: initialized to the global observable rate, then
    re-fit from the off-block pairs after each pass.

    Only pairs whose total coverage (``N``, the sum of the 2x2
    confusion matrix) is at least ``min_pair_coverage`` are used. This
    is the *only* place masking is applied, and it is applied before
    anything else: a pair below the threshold -- including every pair
    touching a fully masked position, which has none above it -- is
    dropped from ``obs_table`` here and so never reaches
    ``_pair_band_row_cumsum``. It is therefore scattered into neither
    the observable nor the significant grid, and cannot contribute to
    any triangle count, background rate, or block gain computed
    downstream; it is not merely down-weighted or diluted, it simply
    does not exist for the rest of this function.

    Represents the grid banded (see ``_pair_band_row_cumsum``), using
    O(L * max_gap) memory where L is the scanned region's length and
    max_gap is the largest separation of any observable pair -- for a
    reference tens of thousands of positions long with a band a few
    hundred wide, tens of megabytes instead of the tens of gigabytes a
    dense O(L^2) grid would need. The DP itself is still O(L^2) time
    (a domain may be far longer than the band), unchanged from a dense
    grid.

    After the DP converges, adjacent blocks are greedily merged,
    strongest bridge first (see ``_merge_adjacent_blocks``), whenever
    the material bridging them is itself denser than background, even
    though the DP's exact optimum kept them split: real regions with
    substantial internal density variation are otherwise reported as
    several small, tightly-adjacent (often zero-gap) domains, which is
    a true statistical result but a finer partition than most
    downstream, flat (non-hierarchical) uses want.
    """
    import numpy as np

    if min_pair_coverage < 1:
        raise OutOfBoundsError(
            f"min_pair_coverage must be >= 1, but got {min_pair_coverage}"
        )
    if bic_multiplier < 0:
        raise OutOfBoundsError(f"bic_multiplier must be >= 0, but got {bic_multiplier}")
    obs_table = table[table[N_COL] >= min_pair_coverage]
    with logger.debug.single_context(
        "Identifying domains from {} observable pair(s) (of {} in-band) "
        "by DP block-diagonal segmentation",
        len(obs_table.index),
        len(table.index),
    ):
        obs_row_cum, sig_row_cum, n_positions, max_gap = _pair_band_row_cumsum(
            obs_table, total_end5, total_end3
        )
        total_obs = int(obs_row_cum[:, max_gap + 1].sum())
        total_sig = int(sig_row_cum[:, max_gap + 1].sum())
        theta0 = total_sig / total_obs if total_obs > 0 else 0.5
        penalty = 0.5 * float(np.log(max(total_obs, 2))) * bic_multiplier
        blocks: list[tuple[int, int]] = []
        for _ in range(N_BACKGROUND_RATE_ITERS):
            blocks = _dp_segment_blocks(
                obs_row_cum, sig_row_cum, n_positions, max_gap, theta0, penalty
            )
            theta0 = _estimate_background_rate(
                obs_row_cum, sig_row_cum, max_gap, blocks, theta0
            )
        blocks = _merge_adjacent_blocks(
            obs_row_cum, sig_row_cum, max_gap, blocks, theta0
        )
        domains = sorted((total_end5 + s, total_end5 + e) for s, e in blocks)
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

    fig = go.Figure()
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
    fig.update_xaxes(title_text="Position", showgrid=True, range=x_range)
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
    band_width: int,
    bic_multiplier: float,
    min_length: int,
    gap_mode: str,
    min_pair_coverage: int = 1000,
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
    # Build the banded per-pair chi-square table across all tiles.
    per_tile_tables = dispatch(
        _find_correlated_pairs,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=True,
        args=as_list_of_tuples(datasets),
        kwargs=dict(pair_fdr=pair_fdr, band_width=band_width),
    )
    table = _build_banded_table(per_tile_tables, band_width)
    pairs = table.index[table[SIGNIFICANT_COL]].to_list()
    # Find domains via DP block-diagonal segmentation.
    domains = _calc_domains_by_dp_segmentation(
        table,
        total_end5=region.end5,
        total_end3=region.end3,
        bic_multiplier=bic_multiplier,
        min_pair_coverage=min_pair_coverage,
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
    band_width: int,
    bic_multiplier: float,
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
            band_width=band_width,
            bic_multiplier=bic_multiplier,
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
            band_width=band_width,
            bic_multiplier=bic_multiplier,
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
