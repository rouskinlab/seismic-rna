from __future__ import annotations
from datetime import datetime
from math import ceil
from pathlib import Path
from typing import Iterable


from .. import filter as filter_mod
from ..core import path
from ..core.batch.confusion import (
    POSITION_A,
    POSITION_B,
    calc_bh_adjusted_pvals,
    calc_confusion_matrix,
)
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
ACTION_COL = "Action"

# How a domain in the final, reported set was produced: called as-is
# (ACTION_ORIGINAL), grown into a gap by --widen (ACTION_WIDENED), or
# inserted into a gap by --fill (ACTION_FILLED).
ACTION_ORIGINAL = "original"
ACTION_WIDENED = "widened"
ACTION_FILLED = "filled"

# Per-pair 0/1 indicator of whether an eligible pair is a "bridge" (the DP's
# discrete currency; see _bridge_mask). Stored as a float column on the
# eligible-pairs table so _pair_band_row_cumsum can sum it into per-block
# bridge counts.
BRIDGE_COL = "Bridge"

# Floor keeping a block's fitted bridge rate strictly inside (0, 1) so the
# binomial log-likelihood ratio stays finite (mirrors the old chi-square
# score's _SCORE_MIN_SHAPE role).
_MIN_RATE = 1e-9

# The four cells of each pair's 2x2 confusion matrix. These are the raw
# counts from which every derived quantity (N and, e.g., fold change) can be
# recomputed; they are written to pairs.csv along with the p-value and fold
# change for convenience.

NEITHER_COL = "Neither Mutated"
ONLY_A_COL = "Only A Mutated"
ONLY_B_COL = "Only B Mutated"
BOTH_COL = "Both Mutated"
CONFUSION_COLS = (NEITHER_COL, ONLY_A_COL, ONLY_B_COL, BOTH_COL)
P_VALUE_COL = "P-value"
Q_VALUE_COL = "Q-value"
FOLD_CHANGE_COL = "Fold Change"


def _expected_both_mutated(table: pd.DataFrame) -> pd.Series:
    """Expected number of jointly-mutated reads under the assumption that
    the two positions mutate independently: with ``a = P(A mutated) * N``
    (``ONLY_A_COL + BOTH_COL``) and ``b`` defined analogously for B,
    ``E[AB] = a * b / N``. This is *not* the observed both-mutated count
    (``BOTH_COL``); it is what independence would predict, which is the
    quantity a chi-square test's minimum-expected-count rule of thumb
    (conventionally >= 5) is about."""
    a = table[ONLY_A_COL] + table[BOTH_COL]
    b = table[ONLY_B_COL] + table[BOTH_COL]
    return a * b / table[N_COL]


def _analyzed_pairs_mask(
    table: pd.DataFrame,
    min_pair_coverage: int,
    min_expect_both: float,
) -> pd.Series:
    """Pairs usable for statistical analysis: enough joint coverage
    (``min_pair_coverage``) and a large enough independence-expected
    both-mutated count (``min_expect_both``) for the exact hypergeometric
    null to be well-defined. This says nothing about the *direction* of
    correlation: both anti-correlated and co-occurring pairs are eligible,
    which is exactly what keeps the one-sided p-values in ``_bridge_mask``
    honestly uniform under the null (dropping the co-occurring half would
    select the small-p tail and inflate the false-discovery rate). Whether a
    pair is anti-correlated enough to be a *bridge* is decided in
    ``_bridge_mask`` by its p-value and effect size."""
    return (table[N_COL] >= min_pair_coverage) & (
        _expected_both_mutated(table) >= min_expect_both
    )


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


def _gather_pooled_reads(dataset: FilterMutsDataset):
    """Read all of one tile's batches exactly once, pooling each
    position's covered and mutated read-index sets across every batch.

    Read numbers are batch-local (a read's identity is the (batch, read)
    pair), so each batch's read numbers are offset by a running stride
    (``max_read + 1``) to make them unique -- and, since batches are
    concatenated in order with a monotonically increasing offset, still
    sorted -- across the pooled tile before the set operations in
    ``calc_confusion_matrix``.

    Pooling first and computing the confusion matrix once over the pooled
    sets gives the identical *observed* matrix as summing one matrix per
    batch (see ``accumulate_confusion_matrices``): every pair's
    intersection count is computed over read numbers that are disjoint
    across batches, so a cross-batch pair contributes nothing and the
    pooled count equals the per-batch sum.
    """
    import numpy as np

    pooled_covered: dict[int, list] = dict()
    pooled_mutated: dict[int, list] = dict()
    offset = 0
    with logger.debug.single_context(
        "Reading {} batch(es) of {}", dataset.num_batches, dataset
    ):
        for batch_num in range(dataset.num_batches):
            batch = dataset.get_batch(batch_num)
            if batch.read_weights is not None:
                raise NotImplementedError(
                    "The confusion-matrix null does not support weighted reads"
                )
            covered = batch.covered_reads_per_pos
            mutated = batch.reads_per_pos(dataset.pattern)
            # Cast to int64 before offsetting: batch read numbers are a small
            # uint type (e.g. uint16) that the running offset would overflow.
            for pos, reads in covered.items():
                pooled_covered.setdefault(pos, []).append(
                    np.asarray(reads, dtype=np.int64) + offset
                )
            for pos, reads in mutated.items():
                pooled_mutated.setdefault(pos, []).append(
                    np.asarray(reads, dtype=np.int64) + offset
                )
            # Offset by max_read + 1 (>= any read number in the batch) so the
            # pooled read numbers stay unique and globally sorted across
            # batches.
            offset += int(batch.max_read) + 1
            logger.trace(
                "{} read batch {}/{} ({} read(s))",
                dataset,
                batch_num + 1,
                dataset.num_batches,
                int(batch.max_read) + 1,
            )
    covered_pooled = {pos: np.concatenate(arrs) for pos, arrs in pooled_covered.items()}
    mutated_pooled = {
        pos: (
            np.concatenate(pooled_mutated[pos])
            if pos in pooled_mutated
            else np.array([], dtype=int)
        )
        for pos in covered_pooled
    }
    logger.debug(
        "Pooled {} read(s) across {} position(s) for {}",
        offset,
        len(covered_pooled),
        dataset,
    )
    return covered_pooled, mutated_pooled


def _confusion_to_table(
    n: pd.Series,
    a: pd.Series,
    b: pd.Series,
    ab: pd.Series,
    *,
    dataset: FilterMutsDataset,
    write_csv: bool,
):
    """Turn one tile's raw 2x2 confusion-matrix components into the pair
    table used downstream: ``N`` and the four raw cells (kept so
    ``pairs.csv`` can record them per bridge pair; every derived quantity,
    e.g. fold change, is recomputable from the cells).

    ``write_csv`` saves ``confusion-matrix.csv`` for the observed data;
    it must be ``False`` for a null replicate, to avoid clobbering the
    real one.
    """
    import pandas as pd

    neither = n - (a + b - ab)
    only_a = a - ab
    only_b = b - ab
    both = ab
    if write_csv:
        csv_file = dataset.report_file.with_name(CONFUSION_MATRIX_CSV)
        confusion_matrix = pd.DataFrame.from_dict(
            {
                NEITHER_COL: neither,
                ONLY_A_COL: only_a,
                ONLY_B_COL: only_b,
                BOTH_COL: both,
            }
        )
        confusion_matrix.to_csv(csv_file)
    table = pd.DataFrame(
        {
            N_COL: n,
            NEITHER_COL: neither,
            ONLY_A_COL: only_a,
            ONLY_B_COL: only_b,
            BOTH_COL: both,
        }
    )
    logger.trace("{} has {} pair(s) within the band", dataset, len(table.index))
    return table


def _find_correlated_pairs(
    dataset: FilterMutsDataset, tile_index: int, *, band_width: int, n_tiles: int
) -> pd.DataFrame:
    """Read one tile's batches (``_gather_pooled_reads``) and compute its
    observed pair table."""
    with logger.debug.single_context(
        "Finding correlated pairs for tile {}/{} ({})", tile_index + 1, n_tiles, dataset
    ):
        pos_index = dataset.region.unmasked
        max_gap = band_width if band_width > 0 else None
        # Widen the exclusion zone to 2 * min_mut_gap (rather than just
        # min_mut_gap) for correlation purposes only: if positions within
        # min_mut_gap of each other are mutually exclusive (by construction
        # of the mutation-collision filter), positions from min_mut_gap + 1
        # to 2 * min_mut_gap show weak positive correlation as an indirect
        # artifact of that exclusivity, not a real domain signal. This does
        # not affect min_mut_gap as used by the actual filter step.
        min_gap = dataset.min_mut_gap * 2
        covered_pooled, mutated_pooled = _gather_pooled_reads(dataset)
        n, a, b, ab = calc_confusion_matrix(
            pos_index, covered_pooled, mutated_pooled, min_gap=min_gap, max_gap=max_gap
        )
        return _confusion_to_table(n, a, b, ab, dataset=dataset, write_csv=True)


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
            {N_COL: pd.Series(dtype=float)}
            | {col: pd.Series(dtype=int) for col in CONFUSION_COLS},
            index=pd.MultiIndex.from_arrays([[], []], names=[POSITION_A, POSITION_B]),
        )
    combined = pd.concat(per_tile_tables, axis=0)
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


def _calc_median_read_length(dataset: MutsDataset, num_cpus: int) -> float:
    """Median length of the reads in one dataset, from a single read of each
    batch's read lengths (see ``_get_batch_read_lengths``). Used to derive the
    default tile length and the default domain-length thresholds."""
    import numpy as np

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
    median_read_length = float(np.median(read_lengths))
    if median_read_length < 1:
        raise ValueError(
            f"The median read length must be ≥ 1, but got {median_read_length}"
        )
    return median_read_length


def _bridge_mask(
    table: pd.DataFrame,
    min_pair_coverage: int,
    min_expect_both: float,
    pair_fdr: float,
    min_fold_change: float,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """``(eligible, bridge, pvalue, qvalue)`` over the pairs of ``table``:
    two boolean masks, and the raw and Benjamini-Hochberg-adjusted p-values
    (NaN where not eligible, since ineligible pairs are excluded from the
    family).

    An **eligible** pair is analyzable (``_analyzed_pairs_mask``: enough joint
    coverage, ``min_pair_coverage``, and expected-both count,
    ``min_expect_both``), regardless of the direction of its correlation.
    Both anti-correlated and co-occurring pairs are eligible on purpose: the
    left-tail p-value below is one-sided, so under the null its values are
    uniform only when *every* analyzable pair is in the family. Restricting
    the family to anti-correlated pairs would select the small-p half and
    inflate the realized false-discovery rate roughly twofold.

    A **bridge** is an eligible pair that is both statistically significant
    and large in effect:

    - *significance*: the exact hypergeometric (Fisher) left-tail p-value
      ``P(X ≤ ab)`` for its 2×2 margins ``(n, a, b)`` -- the probability that
      if ``a`` of the ``n`` jointly-covered reads are mutated at position A
      and ``b`` at position B independently, at most ``ab`` are mutated at
      both. This is the coverage-preserving permutation null evaluated
      exactly, needing no chi-square approximation. Benjamini-Hochberg-
      adjusted over the eligible pairs, kept when ``q < pair_fdr``. Because it
      is a left-tail test, only depleted (anti-correlated) pairs can reach
      small p-values, so significance alone already enforces the direction
      that matters: a per-molecule modification-loading factor scales every
      position together and so can create only positive correlation (its
      contribution to the covariance is ``μ_i·μ_j·Var(loading) ≥ 0``), while
      mutual exclusivity (negative correlation) is the signature of genuine
      alternative structure, immune to that confound.
    - *effect size*: the observed joint-mutation count is at most the
      independence expectation reduced by ``min_fold_change``, i.e.
      ``ab ≤ (a·b)/n / min_fold_change`` (``ab = 0`` always passes). This is
      the fold change ``expected/observed = (a·b)/(n·ab) ≥ min_fold_change``
      rearranged to avoid the division by ``ab``. The fold change reflects the
      cluster mixing (a 50/50 split depletes joint mutations differently than
      a 90/10 split), so a floor on it removes the trivially-significant,
      biologically-negligible pairs that a raw p-value admits at high read
      depth.
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom

    if min_fold_change < 1:
        raise OutOfBoundsError(
            f"min_fold_change must be ≥ 1, but got {min_fold_change}"
        )
    eligible = _analyzed_pairs_mask(
        table, min_pair_coverage, min_expect_both
    )
    bridge = pd.Series(False, index=table.index)
    pvalue = pd.Series(np.nan, index=table.index)
    qvalue = pd.Series(np.nan, index=table.index)
    if not eligible.any():
        return eligible, bridge, pvalue, qvalue
    sub = table[eligible]
    n = sub[N_COL].to_numpy(dtype=np.int64)
    ab = sub[BOTH_COL].to_numpy(dtype=np.int64)
    a = (sub[ONLY_A_COL] + sub[BOTH_COL]).to_numpy(dtype=np.int64)
    b = (sub[ONLY_B_COL] + sub[BOTH_COL]).to_numpy(dtype=np.int64)
    pvals = hypergeom.cdf(ab, n, a, b)
    qvals = calc_bh_adjusted_pvals(pvals)
    pvalue.loc[sub.index] = pvals
    qvalue.loc[sub.index] = qvals
    # Effect size: observed at most the independence expectation reduced by
    # min_fold_change, i.e. ab <= (a*b)/n / min_fold_change. Compare the
    # observed count against this threshold (in float, to avoid integer
    # overflow) rather than forming the fold change, so ab == 0 needs no
    # special-casing -- it is <= any non-negative threshold and always passes.
    expected = (a.astype(float) * b) / n
    effect = ab <= expected / min_fold_change
    is_bridge = (qvals < pair_fdr) & effect
    bridge.loc[sub.index[is_bridge]] = True
    return eligible, bridge, pvalue, qvalue


def _block_score(n_elig: np.ndarray, n_bridge: np.ndarray, pi0: float) -> np.ndarray:
    """Score a candidate block (or array of candidate blocks) holding
    ``n_elig`` eligible pairs of which ``n_bridge`` are bridges, by its own
    locally-fit binomial log-likelihood ratio against the null bridge rate
    ``pi0``, BIC-corrected for the one free parameter it fits.

    Each candidate block fits its *own* bridge rate ``pi_hat = n_bridge /
    n_elig`` from just its own pairs, so a bridge-dense domain elsewhere in the
    scanned region cannot inflate the score of an unrelated block. The rate is
    clamped to ``[pi0, 1 - _MIN_RATE]`` so the score rewards only *enrichment*
    above ``pi0`` (a block at or below the null rate fits ``pi_hat = pi0``,
    giving ``llr = 0`` and a negative gain, exactly as the old chi-square
    score floored its inflation shape at zero)::

        pi_hat = clip(n_bridge / n_elig, pi0, 1 - _MIN_RATE)
        llr    = n_bridge*log(pi_hat/pi0)
                 + (n_elig - n_bridge)*log((1 - pi_hat)/(1 - pi0))

    ``llr`` is the block's own maximized binomial log-likelihood ratio (the MLE
    ``pi_hat`` fits at least as well as ``pi0``), which would let a small subset
    of pairs with an above-average bridge count by pure chance look significant
    against itself. Charging a per-parameter BIC cost (``0.5 * log(n_elig)``,
    for the one free parameter ``pi_hat``) corrects that: a truly null block
    (``pi_hat`` at ``pi0``) scores negative for any ``n_elig > 1``, so extending
    a block into bridge-free space still costs.

    Returns ``gain = llr - log(n_elig) / 2``, the raw (uncharged beyond BIC)
    per-block gain that ``_dp_segment_blocks`` sums directly, with no further
    per-block penalty (see ``_dp_segment_blocks`` for why: gating each block
    against its own admission bar, rather than subtracting a flat fee, keeps
    two adjacent real domains from ever being cheaper to merge than to keep
    separate). An empty block (``n_elig == 0``) scores ``-inf``."""
    import numpy as np

    n_elig = np.asarray(n_elig, dtype=np.float64)
    n_bridge = np.asarray(n_bridge, dtype=np.float64)
    scalar = n_elig.ndim == 0
    n_elig = np.atleast_1d(n_elig)
    n_bridge = np.atleast_1d(n_bridge)
    gain = np.full(n_elig.shape, -np.inf)
    ok = n_elig > 0
    if ok.any():
        ne = n_elig[ok]
        nb = n_bridge[ok]
        with np.errstate(divide="ignore", invalid="ignore"):
            pi_hat = np.clip(nb / ne, pi0, 1.0 - _MIN_RATE)
        llr = nb * np.log(pi_hat / pi0) + (ne - nb) * np.log(
            (1.0 - pi_hat) / (1.0 - pi0)
        )
        bic = np.log(np.maximum(ne, 1.0)) / 2.0
        gain[ok] = llr - bic
    return gain[0] if scalar else gain


def _pair_band_row_cumsum(
    table: pd.DataFrame, total_end5: int, total_end3: int, value_col: str
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Build banded row-cumulative sums of the pair count and of ``value_col``
    (e.g. the per-pair bridge indicator ``BRIDGE_COL``, for domain calling)
    over positions [total_end5, total_end3], 0-indexed, so that the pair count
    and the summed value whose triangle lies within any interval [s, e] (both
    endpoints in [s, e]) can be read off via ``_triangle_sum_banded`` -- see
    ``_calc_domains_by_dp_segmentation`` for how this drives the exact dynamic-
    program segmentation.

    Every pair (a, b) that reaches this function has already passed through
    ``_build_banded_table``'s ``band_width`` cap, so ``b - a`` is bounded by
    some ``max_gap`` read directly off the table's own data. Representing the
    grid banded, as an ``(n_positions, max_gap + 1)`` array of gaps rather than
    a dense ``(n_positions, n_positions)`` array of positions, cuts memory from
    O(L^2) to O(L * max_gap).

    ``table`` is assumed already restricted to the eligible pairs (see
    ``_calc_domains_by_dp_segmentation``) and to carry ``value_col``.

    Returns ``(count_row_cum, value_row_cum, n_positions, max_gap)``, where
    ``count_row_cum[a, d]`` (for ``d`` in ``0..max_gap+1``) is the number of
    pairs ``(a, a+d')`` for ``d' = 0..d-1`` and ``value_row_cum[a, d]`` is the
    sum of their ``value_col`` values (``row_cum[a, 0] == 0``;
    ``row_cum[a, max_gap+1]`` is position ``a``'s full row total within the
    band).
    """
    import numpy as np

    n_positions = total_end3 - total_end5 + 1
    pos_a = table.index.get_level_values(POSITION_A).to_numpy() - total_end5
    pos_b = table.index.get_level_values(POSITION_B).to_numpy() - total_end5
    if len(pos_a) and (
        pos_a.min() < 0 or pos_b.max() >= n_positions or (pos_b < pos_a).any()
    ):
        raise OutOfBoundsError(
            "All pairs must satisfy total_end5 ≤ a ≤ b ≤ total_end3 "
            f"(0 ≤ a ≤ b < {n_positions} after offsetting by "
            f"total_end5={total_end5})"
        )
    value = table[value_col].to_numpy(dtype=float)
    gaps = pos_b - pos_a
    max_gap = int(gaps.max()) if len(gaps) else 0

    count_grid = np.zeros((n_positions, max_gap + 1), dtype=np.int64)
    value_grid = np.zeros((n_positions, max_gap + 1), dtype=float)
    np.add.at(count_grid, (pos_a, gaps), 1)
    np.add.at(value_grid, (pos_a, gaps), value)

    count_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
    value_row_cum = np.zeros((n_positions, max_gap + 2), dtype=float)
    count_row_cum[:, 1:] = np.cumsum(count_grid, axis=1)
    value_row_cum[:, 1:] = np.cumsum(value_grid, axis=1)
    return count_row_cum, value_row_cum, n_positions, max_gap


def _triangle_sum_banded(
    row_cum: np.ndarray, e: int, max_gap: int, s_min: int = 0
) -> np.ndarray:
    """Return an array of length e+1: the triangle sum over every candidate
    block start ``s`` in ``[s_min, e]`` (0-indexed) whose block is ``[s, e]``,
    given a banded row-cumulative array from ``_pair_band_row_cumsum``. Entries
    below ``s_min`` are left zero (the caller slices ``[s_min : e+1]``).

    Bounding the computed range to ``[s_min, e]`` (rather than ``[0, e]``) is
    what lets ``_dp_segment_blocks`` cap each ``e``'s work at
    O(max_domain_length + max_gap) instead of O(e): with ``s_min = e -
    max_domain_length + 1`` the far (full-row) prefix is only summed over the
    admissible starts, not the whole region. With the default ``s_min = 0`` the
    result is identical to summing over all starts.

    Positions ``s <= e - max_gap`` contribute their *full* row (every pair they
    have is within [s, e], since gap <= max_gap <= e - s); positions within
    ``max_gap`` of ``e`` contribute a *partial* row up to distance ``e - s``.
    """
    import numpy as np

    s_min = max(0, s_min)
    true_boundary = max(0, e - max_gap + 1)  # first start with a partial row
    boundary = max(s_min, true_boundary)  # first "near" start we compute
    result = np.zeros(e + 1, dtype=row_cum.dtype)
    # Near/partial region: s in [boundary, e], each contributing row_cum up to
    # distance e - s.
    near_a = np.arange(boundary, e + 1)
    near_d = e - near_a
    near_vals = row_cum[near_a, near_d + 1]
    near_suffix = np.cumsum(near_vals[::-1])[::-1]  # near_suffix[i]=sum(vals[i:])
    if len(near_a):
        result[boundary : e + 1] = near_suffix
    near_total = near_suffix[0] if len(near_suffix) else row_cum.dtype.type(0)
    # Far/full region: s in [s_min, boundary), each contributing its full row
    # plus everything in the near region.
    if boundary > s_min:
        full_row = row_cum[:, max_gap + 1]
        far_vals = full_row[s_min:boundary]
        far_suffix = np.cumsum(far_vals[::-1])[::-1]
        result[s_min:boundary] = far_suffix + near_total
    return result


def _calc_block_pvalue_cutoff(
    elig_row_cum: np.ndarray,
    bridge_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    pi0: float,
    detect_fdr: float,
    max_domain_length: int,
) -> float:
    """Benjamini-Hochberg p-value cutoff over every candidate block the DP
    could ever consider (every ``[s, e]`` with length ``<= max_domain_length``),
    computed analytically with no null replicates.

    Each candidate block's own exact null p-value is ``binom.sf(n_bridge - 1,
    n_elig, pi0)`` -- the probability that at least ``n_bridge`` of its
    ``n_elig`` eligible pairs would be bridges if each were an independent
    bridge with probability ``pi0`` (the null bridge rate). Validated this
    session: under a per-position-independent null the bridge rate is
    essentially zero and flat across coverage, so a single scalar ``pi0`` and
    an independent-Bernoulli block model are justified.

    This complements ``_block_score``'s admission bar: a small "lucky" block
    can score positively under the BIC-corrected LLR, but is unlikely to
    survive Benjamini-Hochberg correction against the full candidate pool.

    Returns the largest raw p-value whose BH-adjusted q-value is
    ``<= detect_fdr``, or ``-1.0`` if none survive (so no p-value, always
    ``>= 0``, is ever admitted -- i.e. "admit nothing")."""
    import numpy as np
    from scipy.stats import binom

    all_ne = []
    all_nb = []
    for t in range(1, n_positions + 1):
        e = t - 1
        s_lo = max(0, t - max_domain_length)
        ne = _triangle_sum_banded(elig_row_cum, e, max_gap, s_min=s_lo)[s_lo:t]
        nb = np.rint(
            _triangle_sum_banded(bridge_row_cum, e, max_gap, s_min=s_lo)[s_lo:t]
        ).astype(np.int64)
        mask = ne > 0
        all_ne.append(ne[mask])
        all_nb.append(nb[mask])
    ne_all = np.concatenate(all_ne) if all_ne else np.zeros(0, dtype=np.int64)
    nb_all = np.concatenate(all_nb) if all_nb else np.zeros(0, dtype=np.int64)
    if ne_all.size == 0:
        return -1.0
    pvals = binom.sf(nb_all - 1, ne_all, pi0)
    qvals = calc_bh_adjusted_pvals(pvals)
    significant = pvals[qvals <= detect_fdr]
    if significant.size == 0:
        return -1.0
    return float(significant.max())


def _dp_segment_blocks(
    elig_row_cum: np.ndarray,
    bridge_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    pi0: float,
    p_cutoff: float,
    max_domain_length: int,
) -> tuple[list[tuple[int, int]], float]:
    """Find the partition of [0, n_positions) into background and domain blocks
    that maximizes the total admitted block score, via exact dynamic
    programming: ``dp[t] = max(dp[t-1], max_s dp[s] + gain(s, t-1))`` over every
    candidate block ``(s, t-1)`` of length ``<= max_domain_length`` whose own
    ``gain`` (``_block_score``, fit locally from that block's own eligible and
    bridge counts) is ``>= 0`` *and* whose own exact binomial p-value clears the
    Benjamini-Hochberg cutoff ``p_cutoff`` (``_calc_block_pvalue_cutoff``);
    ``dp[t]`` is the best score using positions [0, t), and ``gain(s, t-1)`` is
    summed *raw*, with no per-block subtraction.

    ``gain >= 0`` needs no calibration: ``_block_score``'s BIC charge already
    makes a truly null block score negative for any ``n_elig > 1``. Gating on a
    fixed bar instead of subtracting a penalty means picking more blocks is
    never penalized -- each just has to clear the bar -- so two strong domains
    keep their own high scores when kept separate, while a diluted merged block
    can only lose out.

    Bounding candidate starts to ``s >= t - max_domain_length`` enforces the
    hard maximum domain length and cuts the cost from O(n_positions^2) to
    O(n_positions * max_domain_length).

    Returns ``(blocks, objective)``: the chosen blocks as sorted,
    non-overlapping 0-indexed (start, end) pairs, inclusive, and the maximized
    objective ``dp[n_positions]``."""
    import numpy as np
    from scipy.stats import binom

    dp = np.zeros(n_positions + 1)
    backtrack = np.full(n_positions + 1, -1, dtype=np.int64)
    for t in range(1, n_positions + 1):
        e = t - 1
        s_lo = max(0, t - max_domain_length)
        ne = _triangle_sum_banded(elig_row_cum, e, max_gap, s_min=s_lo)[s_lo:t]
        nb = np.rint(
            _triangle_sum_banded(bridge_row_cum, e, max_gap, s_min=s_lo)[s_lo:t]
        ).astype(np.int64)
        s_vals = np.arange(s_lo, t)
        gain = _block_score(ne, nb, pi0)
        with np.errstate(divide="ignore", invalid="ignore"):
            pvals = binom.sf(nb - 1, np.maximum(ne, 1), pi0)
        eligible = (ne > 0) & (gain >= 0.0) & (pvals <= p_cutoff)
        candidate_vals = np.where(eligible, dp[s_vals] + gain, -np.inf)
        best_local = int(np.argmax(candidate_vals))
        best_val = candidate_vals[best_local]
        if best_val > dp[e]:
            dp[t] = best_val
            backtrack[t] = s_vals[best_local]
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
    return blocks, float(dp[n_positions])


def _cut_crossing_scores(
    elig_table: pd.DataFrame, total_end5: int, total_end3: int
) -> tuple[np.ndarray, np.ndarray]:
    """For every cut -- the boundary between position ``m`` and ``m + 1``, for
    ``m`` in ``[total_end5, total_end3)`` -- return the raw ``(n_elig_cross,
    n_bridge_cross)``: the number of eligible pairs ``(a, b)`` with ``a <= m <
    b`` (straddling that cut) and how many of them are bridges.

    Computed as an ``O(n_pairs + L)`` difference-array sweep: each pair ``(a,
    b)`` straddles every cut ``m = a .. b - 1``, so it contributes ``+1`` (and
    ``+bridge``) at cut ``a`` and ``-1`` (and ``-bridge``) at cut ``b``, before
    a single cumulative sum -- not via inclusion-exclusion on triangle sums, so
    a cut with truly zero crossing pairs reads back as exactly zero.

    This is what ``_merge_connected_blocks`` uses to decide whether two
    DP-chosen blocks separated by a gap are directly connected: every cut sees
    every pair that spans it, including one reaching from deep inside one block,
    past an intervening block, into another -- the long-range evidence a block
    (which sees only its own contained pairs) cannot.

    Returns ``(n_elig_cross, n_bridge_cross)``, each of length ``total_end3 -
    total_end5`` (one entry per cut, 0-indexed by ``m - total_end5``)."""
    import numpy as np

    n_positions = total_end3 - total_end5 + 1
    n_cuts = max(n_positions - 1, 0)
    pos_a = (
        elig_table.index.get_level_values(POSITION_A).to_numpy(dtype=np.int64)
        - total_end5
    )
    pos_b = (
        elig_table.index.get_level_values(POSITION_B).to_numpy(dtype=np.int64)
        - total_end5
    )
    bridge_val = elig_table[BRIDGE_COL].to_numpy(dtype=float)
    diff_n = np.zeros(n_positions + 1, dtype=np.int64)
    diff_b = np.zeros(n_positions + 1, dtype=float)
    np.add.at(diff_n, pos_a, 1)
    np.add.at(diff_n, pos_b, -1)
    np.add.at(diff_b, pos_a, bridge_val)
    np.add.at(diff_b, pos_b, -bridge_val)
    n_elig_cross = np.cumsum(diff_n)[:n_cuts]
    n_bridge_cross = np.cumsum(diff_b)[:n_cuts]
    return n_elig_cross, n_bridge_cross


def _merge_connected_blocks(
    blocks: list[tuple[int, int]],
    n_elig_cross: np.ndarray,
    n_bridge_cross: np.ndarray,
    pi0: float,
    merge_fdr: float,
    max_domain_length: int,
) -> list[tuple[int, int]]:
    """Merge adjacent DP-chosen blocks whose intervening gap is, at *every* cut
    within it, directly connected by real crossing bridges -- a targeted test
    for a real gap that nonetheless has real correlations crossing it (e.g. a
    helix with an unpaired internal loop), which the DP's own contained-block
    objective cannot see.

    A cut is connected iff its crossing bridges (``_cut_crossing_scores``) clear
    the same two-gate admission standard used to call a domain:
    ``_block_score(n_elig_cross, n_bridge_cross, pi0) >= 0`` (effect size) *and*
    an exact binomial p-value (``binom.sf(n_bridge_cross - 1, n_elig_cross,
    pi0)``) that survives Benjamini-Hochberg correction across every cut in the
    scanned region, at ``merge_fdr``. Both gates are required: a cut's crossing
    eligible-pair count is band-sized, so BH-significance alone would flag a
    biologically-trivial rate bump as "connected" and over-merge; the effect-
    size floor (crossing bridge *rate* significantly above ``pi0``, not merely
    nonzero) holds true boundaries.

    Two adjacent blocks merge iff every cut in their gap is connected **and**
    the combined block would not exceed ``max_domain_length`` (a hard cap absent
    from the original design: a genuinely-connected span longer than the cap is
    left as multiple domains, since it must fit the clustering window). The
    length guard makes the left-to-right sweep the defined merge order."""
    import numpy as np
    from scipy.stats import binom

    if len(blocks) < 2:
        return blocks
    valid = n_elig_cross > 0
    ne = np.maximum(n_elig_cross, 1)
    nb = np.rint(n_bridge_cross).astype(np.int64)
    gains = np.where(valid, _block_score(ne, nb, pi0), -1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pvals = np.where(valid, binom.sf(nb - 1, ne, pi0), 1.0)
    qvals = np.ones_like(pvals)
    if valid.any():
        qvals[valid] = calc_bh_adjusted_pvals(pvals[valid])
    connected_cut = valid & (gains >= 0.0) & (qvals <= merge_fdr)
    merged: list[tuple[int, int]] = [blocks[0]]
    n_connected = 0
    for i in range(len(blocks) - 1):
        gap_connected = bool(connected_cut[blocks[i][1] : blocks[i + 1][0]].all())
        combined_length = blocks[i + 1][1] - merged[-1][0] + 1
        if gap_connected and combined_length <= max_domain_length:
            merged[-1] = (merged[-1][0], blocks[i + 1][1])
            n_connected += 1
        else:
            merged.append(blocks[i + 1])
    logger.debug(
        "Cut-crossing merge: joined {} of {} adjacent block pair(s) (merge FDR {})",
        n_connected,
        len(blocks) - 1,
        merge_fdr,
    )
    return merged


def _extend_domains_by_bridges(
    domains: list[tuple[int, int]],
    elig_table: pd.DataFrame,
    pi0: float,
    extend_fdr: float,
    max_domain_length: int,
    total_end5: int,
    total_end3: int,
) -> list[tuple[int, int]]:
    """Grow each domain outward to absorb the "thin line" of bridge pairs that
    reach into it from just outside its edge -- long-range contacts that the
    density-optimizing DP excludes because extending over the sparse space
    between them and the dense core would dilute the block's mean bridge rate.

    A bridge ``(a, b)`` links two anti-correlated positions, so if one endpoint
    is inside a domain the other belongs to the same structural module. For each
    edge, extend to the *farthest* outer endpoint ``p`` of the domain's crossing
    bridges such that the crossing bridges over the whole candidate extension --
    pairs from ``[p, s)`` into ``[s, e]`` on the left, symmetric on the right --
    are *collectively* significantly enriched over ``pi0`` (the same two-gate
    admission standard as the merge: ``_block_score >= 0`` and an exact binomial
    p-value ``< extend_fdr``). Testing the aggregate, not each cut, is essential:
    a thin line's per-cut crossing count sits right at the BIC gain margin (it
    flickers between admit/reject cut to cut), but the whole line's bridges are
    unambiguously significant as a group. Extending to the farthest *significant*
    endpoint (not merely the farthest endpoint) ignores a lone distant stray
    bridge, whose singleton aggregate is not significant.

    Extension never crosses into a neighbouring domain and never produces a
    domain longer than ``max_domain_length``. Iterated to a fixpoint, since
    absorbing one region can expose further crossing bridges beyond it."""
    import numpy as np
    from scipy.stats import binom

    if not domains:
        return domains
    pos_a = elig_table.index.get_level_values(POSITION_A).to_numpy()
    pos_b = elig_table.index.get_level_values(POSITION_B).to_numpy()
    is_bridge = elig_table[BRIDGE_COL].to_numpy() > 0.5

    def _significant(mask: np.ndarray) -> bool:
        n_elig = int(mask.sum())
        if n_elig == 0:
            return False
        n_bridge = int(is_bridge[mask].sum())
        if _block_score(np.array([n_elig]), np.array([n_bridge]), pi0)[0] < 0.0:
            return False
        return bool(binom.sf(n_bridge - 1, n_elig, pi0) < extend_fdr)

    doms = [list(d) for d in sorted(domains)]
    changed = True
    while changed:
        changed = False
        for i, (s, e) in enumerate(doms):
            left_bound = doms[i - 1][1] + 1 if i > 0 else total_end5
            right_bound = doms[i + 1][0] - 1 if i < len(doms) - 1 else total_end3
            # Left: outer endpoints a < s of bridges reaching into [s, e],
            # tried farthest-first; extend to the farthest whose whole
            # [p, s) -> [s, e] crossing is significant.
            left_ends = pos_a[is_bridge & (pos_a < s) & (pos_b >= s) & (pos_b <= e)]
            for p in sorted(set(int(x) for x in left_ends)):
                if p < left_bound or (e - p + 1) > max_domain_length:
                    continue
                conn = (pos_a >= p) & (pos_a < s) & (pos_b >= s) & (pos_b <= e)
                if _significant(conn):
                    doms[i][0] = p
                    changed = True
                    break
            s = doms[i][0]
            # Right: symmetric, outer endpoints b > e tried farthest-first.
            right_ends = pos_b[is_bridge & (pos_b > e) & (pos_a >= s) & (pos_a <= e)]
            for p in sorted((int(x) for x in right_ends), reverse=True):
                if p > right_bound or (p - s + 1) > max_domain_length:
                    continue
                conn = (pos_b > e) & (pos_b <= p) & (pos_a >= s) & (pos_a <= e)
                if _significant(conn):
                    doms[i][1] = p
                    changed = True
                    break
    return [(d[0], d[1]) for d in doms]


def _estimate_null_bridge_rate(
    elig_table: pd.DataFrame, domains: list[tuple[int, int]], floor: float
) -> float:
    """Estimate the null bridge rate ``pi0`` as the bridge rate among the
    eligible pairs lying *outside* the (preliminary) domains: with no real
    domain there, those bridges are false positives, so their rate is the
    background rate a real domain must beat. Captures real background and bias
    that a permutation or pure-FDR estimate misses, and excludes the domains so
    a dense domain cannot inflate it (immune to the global-rate contamination
    of the previous model).

    ``domains`` are absolute ``(end5, end3)`` intervals; a pair ``(a, b)`` is
    "inside" a domain iff ``end5 <= a`` and ``b <= end3``. Floored at ``floor``
    (the pass-1 seed ``pair_fdr * global_rate``) so an all-clean background
    never yields ``pi0 = 0``."""
    import numpy as np

    pos_a = elig_table.index.get_level_values(POSITION_A).to_numpy()
    pos_b = elig_table.index.get_level_values(POSITION_B).to_numpy()
    bridge_val = elig_table[BRIDGE_COL].to_numpy(dtype=float)
    inside = np.zeros(len(pos_a), dtype=bool)
    for lo, hi in domains:
        inside |= (pos_a >= lo) & (pos_b <= hi)
    outside = ~inside
    n_elig_out = int(outside.sum())
    if n_elig_out == 0:
        return floor
    rate = float(bridge_val[outside].sum()) / n_elig_out
    return max(rate, floor)


def _filter_domains_length(domains: list[tuple[int, int]], min_length: int = 1):
    """Remove domains shorter than ``min_length`` positions."""
    return [(end5, end3) for end5, end3 in domains if (end3 - end5 + 1) >= min_length]


def _filter_domains_min_pairs(
    domains: list[tuple[int, int]],
    pos_a: np.ndarray,
    pos_b: np.ndarray,
    min_pairs: int,
) -> list[tuple[int, int]]:
    """Remove domains containing fewer than ``min_pairs`` bridge pairs
    (``pos_a``, ``pos_b``: the endpoints of every bridge pair), i.e. pairs
    whose both endpoints lie within the domain."""
    import numpy as np

    return [
        (end5, end3)
        for end5, end3 in domains
        if int(((pos_a >= end5) & (pos_b <= end3)).sum()) >= min_pairs
    ]


def _split_gap_evenly(
    end5: int, end3: int, max_domain_length: int
) -> list[tuple[int, int]]:
    """Split [end5, end3] into the minimum number of domains of as equal
    length as possible, none exceeding max_domain_length."""
    length = end3 - end5 + 1
    n = -(-length // max_domain_length)  # ceil division
    base, extra = divmod(length, n)
    domains = []
    pos = end5
    for i in range(n):
        size = base + (1 if i < extra else 0)
        domains.append((pos, pos + size - 1))
        pos += size
    return domains


def _fill_domains_into_gaps(
    domains: list[tuple[int, int]],
    global_end5: int,
    global_end3: int,
    max_domain_length: int,
) -> list[tuple[int, int]]:
    """Insert domains into every gap -- leading (before the first domain),
    interior (between two consecutive domains), and trailing (after the
    last domain) -- so that afterward every position in
    [global_end5, global_end3] is covered by exactly one domain. Splits any
    gap longer than max_domain_length via _split_gap_evenly."""

    def _fill_one_gap(gap_end5: int, gap_end3: int) -> list[tuple[int, int]]:
        if gap_end3 < gap_end5:
            return []
        if gap_end3 - gap_end5 + 1 > max_domain_length:
            return _split_gap_evenly(gap_end5, gap_end3, max_domain_length)
        return [(gap_end5, gap_end3)]

    if not domains:
        return _fill_one_gap(global_end5, global_end3)
    filled = _fill_one_gap(global_end5, domains[0][0] - 1)
    filled.append(domains[0])
    for (_, prev_end3), (next_end5, next_end3) in zip(domains, domains[1:]):
        filled.extend(_fill_one_gap(prev_end3 + 1, next_end5 - 1))
        filled.append((next_end5, next_end3))
    filled.extend(_fill_one_gap(domains[-1][1] + 1, global_end3))
    return filled


def _widen_domains_into_gaps(
    domains: list[tuple[int, int]],
    global_end5: int,
    global_end3: int,
    max_domain_length: int,
) -> list[tuple[int, int]]:
    """Grow each domain into adjacent gap space (leading edge, each interior
    gap, trailing edge, processed in that order), capping every domain's
    total length at max_domain_length. An interior gap splits its space
    between its two neighbors in half, each capped by that neighbor's own
    remaining budget; any space neither neighbor can absorb is left as a
    (possibly shorter) residual gap."""
    if not domains:
        return []
    end5s = [d[0] for d in domains]
    end3s = [d[1] for d in domains]
    n = len(domains)

    def budget(i: int) -> int:
        return max_domain_length - (end3s[i] - end5s[i] + 1)

    lead_gap = end5s[0] - global_end5
    if lead_gap > 0:
        end5s[0] -= min(lead_gap, budget(0))
    for i in range(n - 1):
        gap = end5s[i + 1] - end3s[i] - 1
        if gap > 0:
            desired_left = gap // 2
            desired_right = gap - desired_left
            b_left, b_right = budget(i), budget(i + 1)
            actual_left = min(desired_left, b_left)
            actual_right = min(desired_right, b_right)
            leftover = gap - actual_left - actual_right
            if leftover > 0:
                extra = min(leftover, b_left - actual_left)
                actual_left += extra
                leftover -= extra
            if leftover > 0:
                extra = min(leftover, b_right - actual_right)
                actual_right += extra
            end3s[i] += actual_left
            end5s[i + 1] -= actual_right
    trail_gap = global_end3 - end3s[-1]
    if trail_gap > 0:
        end3s[-1] += min(trail_gap, budget(n - 1))
    return _ends_arrays_to_tuples(end5s, end3s)


def _label_widened(
    raw_domains: list[tuple[int, int]], widened_domains: list[tuple[int, int]]
) -> dict[tuple[int, int], str]:
    """Pair raw (pre-widen) domains with their widened counterparts (same
    count/order, since widening only grows boundaries) and label each
    ACTION_ORIGINAL (unchanged) or ACTION_WIDENED (grown)."""
    return {
        wide: (ACTION_ORIGINAL if wide == orig else ACTION_WIDENED)
        for orig, wide in zip(raw_domains, widened_domains)
    }


def _calc_domains_by_dp_segmentation(
    table: pd.DataFrame,
    eligible: pd.Series,
    bridge: pd.Series,
    total_end5: int,
    total_end3: int,
    pair_fdr: float,
    detect_fdr: float,
    merge_fdr: float,
    max_domain_length: int,
) -> tuple[list[tuple[int, int]], float]:
    """Call domains by exact dynamic-program segmentation of the per-pair
    *bridge indicator* into background/domain blocks (see ``_dp_segment_blocks``):
    partition [total_end5, total_end3] into consecutive intervals and keep as a
    domain each interval whose own locally-fit binomial log-likelihood ratio
    (``_block_score``: its own bridge rate ``pi_hat`` vs the null rate ``pi0``,
    BIC-charged) is ``>= 0`` *and* whose own exact binomial p-value clears a
    Benjamini-Hochberg cutoff over every candidate this scan considers
    (``_calc_block_pvalue_cutoff``), then merge adjacent blocks connected across
    a gap by crossing bridges (``_cut_crossing_scores`` /
    ``_merge_connected_blocks``). Every candidate block fits its *own* bridge
    rate, so a bridge-dense domain elsewhere cannot inflate an unrelated block.

    ``pi0`` (the null bridge rate) is estimated from the data in one refinement
    pass rather than assumed: seed ``pi0 = pair_fdr * global_bridge_rate`` (the
    FDR guarantee's expected false-bridge rate), run the DP + merge, then
    re-estimate ``pi0`` as the bridge rate among eligible pairs *outside* the
    preliminary domains (``_estimate_null_bridge_rate``) and run once more.

    Returns ``(domains, pi0)``: the sorted absolute ``(end5, end3)`` domain
    intervals and the estimated null bridge rate used for the final pass."""
    if not 1 <= total_end5 <= total_end3:
        raise IncompatibleValuesError(
            "Must have 1 ≤ total_end5 ≤ total_end3, "
            f"but got total_end5={total_end5} and total_end3={total_end3}"
        )
    if max_domain_length < 1:
        raise OutOfBoundsError(
            f"max_domain_length must be ≥ 1, but got {max_domain_length}"
        )
    elig_table = table.loc[eligible].copy()
    elig_table[BRIDGE_COL] = bridge.loc[eligible].to_numpy(dtype=float)
    n_elig_total = int(eligible.sum())
    n_bridge_total = int(bridge.sum())
    with logger.debug.single_context(
        "Identifying domains from {} eligible pair(s) ({} bridge(s)) "
        "by DP block segmentation",
        n_elig_total,
        n_bridge_total,
    ):
        # Grids and cut-crossing scores depend only on the pairs, not on pi0,
        # so compute them once and reuse across both DP passes.
        elig_row_cum, bridge_row_cum, n_positions, max_gap = _pair_band_row_cumsum(
            elig_table, total_end5, total_end3, value_col=BRIDGE_COL
        )
        n_elig_cross, n_bridge_cross = _cut_crossing_scores(
            elig_table, total_end5, total_end3
        )
        global_rate = n_bridge_total / n_elig_total if n_elig_total else 0.0
        pi0_seed = max(pair_fdr * global_rate, _MIN_RATE)

        def _run(pi0: float) -> list[tuple[int, int]]:
            p_cutoff = _calc_block_pvalue_cutoff(
                elig_row_cum,
                bridge_row_cum,
                n_positions,
                max_gap,
                pi0,
                detect_fdr,
                max_domain_length,
            )
            blocks, _ = _dp_segment_blocks(
                elig_row_cum,
                bridge_row_cum,
                n_positions,
                max_gap,
                pi0,
                p_cutoff,
                max_domain_length,
            )
            blocks = _merge_connected_blocks(
                blocks, n_elig_cross, n_bridge_cross, pi0, merge_fdr, max_domain_length
            )
            domains = sorted((total_end5 + s, total_end5 + e) for s, e in blocks)
            # Grow domains outward to absorb thin-line edge bridges the density
            # DP excludes (uses merge_fdr, the same crossing-bridge standard).
            return _extend_domains_by_bridges(
                domains,
                elig_table,
                pi0,
                merge_fdr,
                max_domain_length,
                total_end5,
                total_end3,
            )

        # Pass 1: seed pi0 from the FDR guarantee, find preliminary domains.
        prelim_domains = _run(pi0_seed)
        # Pass 2: re-estimate pi0 as the out-of-domain background bridge rate.
        pi0 = _estimate_null_bridge_rate(elig_table, prelim_domains, floor=pi0_seed)
        domains = sorted(_run(pi0))
        logger.debug(
            "Calculated {} domain(s) at pi0={} (seed {}): {}",
            len(domains),
            pi0,
            pi0_seed,
            domains,
        )
        return domains, pi0


def _graph_pairs_and_domains(
    pairs: list[tuple[int, int]],
    fold_changes: list[float],
    domains: list[tuple[int, int]],
    end5: int,
    end3: int,
    html_file: str | Path,
):
    import numpy as np
    from plotly import graph_objects as go

    fig = go.Figure()
    # Graph the domains as triangles. Only the first trace shows a legend
    # entry, so all domains appear under a single "Domains" label.
    end5s, end3s = _tuples_to_ends_arrays(domains)
    domains_midpoints, domains_distances = _calc_midpoints_distances(end5s, end3s)
    for i, ((a, b), x, y) in enumerate(
        zip(domains, domains_midpoints, domains_distances, strict=True)
    ):
        fig.add_trace(
            go.Scatter(
                x=[a, b, x, a],
                y=[0, 0, y, 0],
                mode="none",
                fill="toself",
                fillcolor="rgba(230,159,0,0.5)",
                showlegend=(i == 0),
                name="Domains",
                hoverinfo="text",
                text=f"Domain {a, b}",
                hovertemplate="%{text}<extra></extra>",
            )
        )
    # Plot the bridge pairs as points, uniformly colored the fully-opaque
    # shade that used to sit at the top of the fold-change color scale
    # (#D55E00).
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    pairs_midpoints, pairs_distances = _calc_midpoints_distances(pos5s, pos3s)
    pairs_color = "rgb(213,94,0)"  # #D55E00
    marker = dict(color=pairs_color)
    pairs_text = [
        f"Pair {pair}, fold change={fc:.2f}"
        for pair, fc in zip(pairs, fold_changes, strict=True)
    ]
    fig.add_trace(
        go.Scatter(
            x=pairs_midpoints,
            y=pairs_distances,
            mode="markers",
            showlegend=True,
            marker=marker,
            name="Bridge pairs",
            hoverinfo="text",
            text=pairs_text,
            hovertemplate="%{text}<extra></extra>",
        )
    )
    # A static legend key (not an interactive toggle): placed above the plot.
    fig.update_layout(
        legend=dict(
            itemclick=False,
            itemdoubleclick=False,
            orientation="h",
            x=0.0,
            y=1.08,
            xanchor="left",
            yanchor="bottom",
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


def _write_pairs_to_csv(pairs: list[tuple[int, int]], csv_file: str | Path):
    """Write the pairs to a CSV file."""
    import pandas as pd

    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    df = pd.DataFrame.from_dict(
        {POSITION_A: pos5s, POSITION_B: pos3s}, orient="columns"
    )
    df.to_csv(csv_file, index=False)


def _write_domains_to_csv(
    domain_actions: dict[tuple[int, int], str], csv_file: str | Path
):
    """Write the final domains, and how each was produced, to a CSV file."""
    import pandas as pd

    end5s, end3s = _tuples_to_ends_arrays(list(domain_actions.keys()))
    df = pd.DataFrame.from_dict(
        {
            POSITION_A: end5s,
            POSITION_B: end3s,
            ACTION_COL: list(domain_actions.values()),
        }
    )
    df.to_csv(csv_file, index=False)


def _write_pairs_with_confusion(
    pos_table: pd.DataFrame,
    pvalue: pd.Series,
    qvalue: pd.Series,
    csv_file: str | Path,
):
    """Write the given pairs' 2x2 confusion-matrix counts to a CSV file,
    along with the exact hypergeometric (Fisher) left-tail p-value
    (``pvalue``) and its Benjamini-Hochberg-adjusted q-value (``qvalue``),
    both aligned to ``pos_table``'s index and passed in rather than
    recomputed here -- ``_bridge_mask`` already computed them once, and the
    BH adjustment in particular depends on the full eligible family, not
    just the written pairs -- plus the fold change (independence
    expectation over the observed both-mutated count)."""
    import numpy as np
    import pandas as pd

    with np.errstate(divide="ignore", invalid="ignore"):
        fold_changes = (
            _expected_both_mutated(pos_table) / pos_table[BOTH_COL]
        ).to_numpy()
    df = pd.DataFrame(
        {
            POSITION_A: pos_table.index.get_level_values(POSITION_A).to_numpy(),
            POSITION_B: pos_table.index.get_level_values(POSITION_B).to_numpy(),
        }
        | {col: pos_table[col].to_numpy() for col in CONFUSION_COLS}
        | {
            P_VALUE_COL: pvalue.loc[pos_table.index].to_numpy(),
            Q_VALUE_COL: qvalue.loc[pos_table.index].to_numpy(),
            FOLD_CHANGE_COL: fold_changes,
        }
    )
    df.to_csv(csv_file, index=False)


def _calc_cluster_domains(
    filter_dirs: list[Path],
    report_dir: Path,
    num_cpus: int,
    band_width: int,
    min_pair_coverage: int,
    min_expect_both: float,
    pair_fdr: float,
    min_fold_change: float,
    detect_fdr: float,
    merge_fdr: float,
    widen: bool,
    fill: bool,
    max_domain_length: int,
    min_domain_length: int,
    min_pairs: int,
):
    """Calculate the cluster regions for all tiles of one reference.

    Domains are called by exact dynamic-program segmentation of the per-pair
    bridge indicator into background/domain blocks
    (``_calc_domains_by_dp_segmentation``): each candidate block is scored by
    its own locally-fit bridge rate against an estimated null rate ``pi0``, the
    globally-optimal partition is found by DP, and adjacent blocks connected by
    crossing bridges are merged. ``widen`` grows domains into their neighboring
    gaps (capped at ``max_domain_length``) and ``fill`` then inserts domains
    into whatever gaps remain (splitting any that exceed ``max_domain_length``).
    """
    import numpy as np

    # Each dataset corresponds to one tile.
    datasets = list(load_filter_dataset.iterate(filter_dirs))
    if not datasets:
        return list(), {}, 0, 0.0
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
    # Build the banded per-pair confusion table across all tiles, from a
    # single read of each tile's batches (see _find_correlated_pairs).
    n_tiles = len(datasets)
    with logger.debug.single_context(
        "Finding correlated pairs across {} tile(s) of {}", n_tiles, ref
    ):
        per_tile_tables = dispatch(
            _find_correlated_pairs,
            num_cpus=num_cpus,
            pass_num_cpus=False,
            as_list=True,
            ordered=True,
            raise_on_error=True,
            args=list(zip(datasets, range(n_tiles))),
            kwargs=dict(band_width=band_width, n_tiles=n_tiles),
        )
    table = _build_banded_table(per_tile_tables, band_width)
    # Compute the eligible pairs and bridge indicator once, and reuse it for
    # both domain calling and reporting. A pair is "positive"
    # (displayed/saved/counted in n_positive_pairs) when it is a bridge -- the
    # same criterion the domain caller itself uses -- so the reported pairs
    # match what actually drove the called domains.
    eligible, bridge, pvalue, qvalue = _bridge_mask(
        table, min_pair_coverage, min_expect_both, pair_fdr, min_fold_change
    )
    pos_table = table[bridge]
    pairs = pos_table.index.to_list()
    with np.errstate(divide="ignore", invalid="ignore"):
        pair_fold_changes = (
            _expected_both_mutated(pos_table) / pos_table[BOTH_COL]
        ).to_list()
    # Call domains by DP block segmentation of the per-pair bridge indicator.
    raw_domains, null_bridge_rate = _calc_domains_by_dp_segmentation(
        table,
        eligible,
        bridge,
        total_end5=region.end5,
        total_end3=region.end3,
        pair_fdr=pair_fdr,
        detect_fdr=detect_fdr,
        merge_fdr=merge_fdr,
        max_domain_length=max_domain_length,
    )
    raw_domains = _filter_domains_min_pairs(
        raw_domains,
        pos_table.index.get_level_values(POSITION_A).to_numpy(),
        pos_table.index.get_level_values(POSITION_B).to_numpy(),
        min_pairs,
    )
    # Widen domains into their neighboring gaps (capped at max_domain_length),
    # then fill whatever gaps remain (splitting any that are still too long).
    if widen:
        working_domains = _widen_domains_into_gaps(
            raw_domains, region.end5, region.end3, max_domain_length
        )
        pre_fill_actions = _label_widened(raw_domains, working_domains)
    else:
        working_domains = _filter_domains_length(
            raw_domains, min_length=min_domain_length
        )
        pre_fill_actions = {d: ACTION_ORIGINAL for d in working_domains}
    final_domains = (
        _fill_domains_into_gaps(
            working_domains, region.end5, region.end3, max_domain_length
        )
        if fill
        else working_domains
    )
    domain_actions = {d: pre_fill_actions.get(d, ACTION_FILLED) for d in final_domains}
    # Write the pairs and domains to CSV files.
    report_dir.mkdir(parents=True, exist_ok=True)
    _write_pairs_with_confusion(
        pos_table, pvalue, qvalue, report_dir.joinpath(PAIRS_CSV)
    )
    _write_domains_to_csv(domain_actions, report_dir.joinpath(DOMAINS_CSV))
    logger.debug(
        "Wrote {} pair(s) and {} domain(s) to {}",
        len(pairs),
        len(final_domains),
        report_dir,
    )
    # Graph the bridge pairs and domains.
    html_file = report_dir.joinpath(PAIRS_DOMAINS_HTML)
    try:
        _graph_pairs_and_domains(
            pairs, pair_fold_changes, final_domains, region.end5, region.end3, html_file
        )
    except Exception as error:
        logger.error(error)
    return (
        [(ref, end5, end3) for end5, end3 in final_domains],
        domain_actions,
        len(pairs),
        null_bridge_rate,
    )


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
    band_width: int,
    min_pair_coverage: int,
    min_expect_both: float,
    pair_fdr: float,
    min_fold_change: float,
    detect_fdr: float,
    merge_fdr: float,
    widen: bool,
    fill: bool,
    max_domain_length: int,
    min_domain_length: int,
    min_pairs: int,
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
        if tile_length <= 0 or max_domain_length <= 0:
            # Median read length determines both the default tile length and
            # the maximum domain length.
            median_read_length = _calc_median_read_length(idmut_dataset, num_cpus)
            logger.trace("The median read length is {}", median_read_length)
            if tile_length <= 0:
                tile_length = round(median_read_length * 2)
                logger.trace("Using tile_length={}", tile_length)
            if max_domain_length <= 0:
                max_domain_length = round(median_read_length * 2)
                logger.trace("Using max_domain_length={}", max_domain_length)
        # Compute tile coordinates.
        with logger.debug.single_context("Calculating tiles for reference {!r}", ref):
            tile_coords = [
                (ref, end5, end3)
                for end5, end3 in _calc_tiles(
                    total_region.end5, total_region.end3, tile_length, tile_min_overlap
                )
            ]
        logger.debug(
            "Began filtering {} tile(s) of {}", len(tile_coords), idmut_report_file
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
        logger.debug(
            "Ended filtering {} tile(s) of {}", len(tile_coords), idmut_report_file
        )
        # Find regions spanned by correlated base pairs.
        report_dir = report_file.parent
        domain_coords, domain_actions, n_positive_pairs, null_bridge_rate = (
            _calc_cluster_domains(
                tiled_dirs,
                report_dir=report_dir,
                band_width=band_width,
                min_pair_coverage=min_pair_coverage,
                min_expect_both=min_expect_both,
                pair_fdr=pair_fdr,
                min_fold_change=min_fold_change,
                detect_fdr=detect_fdr,
                merge_fdr=merge_fdr,
                widen=widen,
                fill=fill,
                max_domain_length=max_domain_length,
                min_domain_length=min_domain_length,
                min_pairs=min_pairs,
                num_cpus=num_cpus,
            )
        )
        logger.debug(
            "Found {} domain(s) and {} positive-score pair(s) in {}",
            len(domain_coords),
            n_positive_pairs,
            idmut_report_file,
        )
        if domain_coords:
            # Filter the reads over each domain so that clusterscan can
            # cluster them without re-running the filter step.
            logger.debug(
                "Began filtering {} domain(s) of {}",
                len(domain_coords),
                idmut_report_file,
            )
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
            logger.debug(
                "Ended filtering {} domain(s) of {}",
                len(domain_coords),
                idmut_report_file,
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
            band_width=band_width,
            min_pair_coverage=min_pair_coverage,
            min_expect_both=min_expect_both,
            pair_fdr=pair_fdr,
            min_fold_change=min_fold_change,
            detect_fdr=detect_fdr,
            merge_fdr=merge_fdr,
            widen=widen,
            fill=fill,
            max_domain_length=max_domain_length,
            min_domain_length=min_domain_length,
            min_pairs=min_pairs,
            # Results (store coordinates without the reference, which is
            # already recorded in the report).
            tile_coords=[(end5, end3) for _, end5, end3 in tile_coords],
            n_positive_pairs=n_positive_pairs,
            null_bridge_rate=null_bridge_rate,
            n_domains=len(domain_coords),
            domain_coords=domain_actions,
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
