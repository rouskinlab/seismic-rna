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
    calc_confusion_matrix,
    calc_confusion_chi_square,
    resample_mutated_reads,
)
from ..core.dataset import MutsDataset
from ..core.error import IncompatibleValuesError, OutOfBoundsError
from ..core.logs import logger
from ..core.random import get_random_integer_generator
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
CHI_SQUARE_COL = "Chi-Square"
# Per-pair block-scoring currency (see _calc_pair_scores): the null-
# referenced, coverage-normalized correlation excess g_ij that the DP
# accumulates over a block's triangle in place of the old significant-pair
# count.
SCORE_COL = "Score"

# The four cells of each pair's 2x2 confusion matrix. These are the raw
# counts from which every derived quantity (N and chi-square) can be
# recomputed, so only these are written to pairs.csv.

# Floor for the estimated domain variance-inflation shape r (see
# _calc_pair_scores): keeps r strictly positive so the score stays defined
# even when the observed data carries essentially no correlation signal
# (mean chi-square at or below the null's 1), in which case the penalty/FDR
# calibration is what declines to call any domain.
_SCORE_MIN_SHAPE = 1e-9

NEITHER_COL = "Neither Mutated"
ONLY_A_COL = "Only A Mutated"
ONLY_B_COL = "Only B Mutated"
BOTH_COL = "Both Mutated"
CONFUSION_COLS = (NEITHER_COL, ONLY_A_COL, ONLY_B_COL, BOTH_COL)


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
    pooled count equals the per-batch sum. This lets one read of the
    tile's batches serve both the observed matrix and every null
    replicate (see ``_find_correlated_pairs_with_nulls``), instead of
    re-reading the batches once per replicate.
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
    covered_pooled = {
        pos: np.concatenate(arrs) for pos, arrs in pooled_covered.items()
    }
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
    table used downstream: chi-square and the four raw cells (kept so
    ``pairs.csv`` can record them per positive-score pair; chi-square is
    used only by the domain caller, not written to ``pairs.csv``, being
    recomputable from the cells).

    ``write_csv`` saves ``confusion-matrix.csv`` for the observed data;
    it must be ``False`` for a null replicate, to avoid clobbering the
    real one.
    """
    import pandas as pd

    chi_square = calc_confusion_chi_square(n, a, b, ab)
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
                CHI_SQUARE_COL: chi_square,
            }
        )
        confusion_matrix.to_csv(csv_file)
    table = pd.DataFrame(
        {
            N_COL: n,
            CHI_SQUARE_COL: chi_square,
            NEITHER_COL: neither,
            ONLY_A_COL: only_a,
            ONLY_B_COL: only_b,
            BOTH_COL: both,
        }
    )
    logger.trace(
        "{} has {} pair(s) within the band",
        dataset,
        len(table.index),
    )
    return table


def _compute_null_confusion_table(
    replicate_seed: int,
    *,
    pos_index,
    covered_pooled: dict,
    mutated_pooled: dict,
    dataset: FilterMutsDataset,
    min_gap: int,
    max_gap: int | None,
) -> pd.DataFrame:
    """Compute one null-replicate confusion-matrix table for a tile from
    its already-pooled reads (see ``_find_correlated_pairs_with_nulls``).

    Factored out of ``_find_correlated_pairs_with_nulls`` so the null
    replicates can be dispatched in parallel. ``replicate_seed`` fully
    determines this replicate's pseudo-random draws; the caller draws one
    such seed per replicate from a per-tile stream so the seeds are
    distinct and reproducible.
    """
    mutated_null = resample_mutated_reads(
        covered_pooled, mutated_pooled, replicate_seed
    )
    n_null, a_null, b_null, ab_null = calc_confusion_matrix(
        pos_index,
        covered_pooled,
        mutated_null,
        None,
        min_gap=min_gap,
        max_gap=max_gap,
    )
    null_table = _confusion_to_table(
        n_null, a_null, b_null, ab_null, dataset=dataset, write_csv=False
    )
    logger.debug(
        "Computed a null-replicate confusion matrix (seed {}) for {}",
        replicate_seed,
        dataset,
    )
    return null_table


def _find_correlated_pairs_with_nulls(
    dataset: FilterMutsDataset,
    tile_index: int,
    tile_seed: int,
    *,
    band_width: int,
    n_null_replicates: int,
    n_tiles: int,
    num_cpus: int,
):
    """Read one tile's batches exactly once (``_gather_pooled_reads``)
    and compute both its observed pair table and every null-replicate
    pair table from that single pooled read.

    Returns ``(real_table, null_tables)``, where ``null_tables`` is a
    list of length ``n_null_replicates`` (empty if 0). The null
    replicates are computed in parallel over ``num_cpus`` processors
    (``_compute_null_confusion_table`` via ``dispatch``); ``ordered=True``
    keeps ``null_tables`` in replicate order. Each replicate's seed is
    drawn from a random-integer stream initialized with this tile's own
    seed ``tile_seed`` (itself drawn by the caller from a stream seeded
    with the run's seed), so every replicate gets its own reproducible
    seed without ever computing one by offsetting another.
    """
    with logger.debug.single_context(
        "Finding correlated pairs for tile {}/{} ({})",
        tile_index + 1,
        n_tiles,
        dataset,
    ):
        pos_index = dataset.region.unmasked
        max_gap = band_width if band_width > 0 else None
        min_gap = dataset.min_mut_gap
        covered_pooled, mutated_pooled = _gather_pooled_reads(dataset)
        n, a, b, ab = calc_confusion_matrix(
            pos_index,
            covered_pooled,
            mutated_pooled,
            None,
            min_gap=min_gap,
            max_gap=max_gap,
        )
        real_table = _confusion_to_table(
            n, a, b, ab, dataset=dataset, write_csv=True
        )
        logger.debug(
            "Computed the observed confusion matrix for {}", dataset
        )
        if n_null_replicates > 0:
            # Draw one seed per replicate from a stream initialized with
            # this tile's own seed, and pass each as the positional argument
            # so every replicate gets its own distinct, reproducible seed.
            seed_stream = get_random_integer_generator(tile_seed)
            replicate_seeds = [
                next(seed_stream) for _ in range(n_null_replicates)
            ]
            null_tables = dispatch(
                _compute_null_confusion_table,
                num_cpus=num_cpus,
                pass_num_cpus=False,
                as_list=True,
                ordered=True,
                raise_on_error=True,
                args=as_list_of_tuples(replicate_seeds),
                kwargs=dict(
                    pos_index=pos_index,
                    covered_pooled=covered_pooled,
                    mutated_pooled=mutated_pooled,
                    dataset=dataset,
                    min_gap=min_gap,
                    max_gap=max_gap,
                ),
            )
        else:
            null_tables = []
        return real_table, null_tables


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
            {col: pd.Series(dtype=float) for col in (N_COL, CHI_SQUARE_COL)}
            | {col: pd.Series(dtype=int) for col in CONFUSION_COLS},
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


def _calc_pair_scores(
    table: pd.DataFrame,
    null_tables_per_replicate: list[pd.DataFrame],
    min_pair_coverage: int,
) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    """Attach the per-pair block-scoring currency ``SCORE_COL`` (``g_ij``)
    to the observed ``table`` and to each null replicate table.

    The score tests, per block, whether the *distribution* of the pairs'
    correlations is drawn from a wider-than-null distribution. Work in
    chi-square space, where the null is coverage-invariant: under the
    per-position-independence null a pair's ``chi^2 = N * phi^2`` is
    ``chi^2(1)`` regardless of coverage. Model a domain as inflating that
    variance by a factor ``1 + r`` (the pair's phi is drawn from a wider
    spread), and score each pair by the log-likelihood ratio of the inflated
    model against the null::

        g_ij = log f_alt(chi2_ij) - log f_null(chi2_ij)
             = -0.5 * log(1 + r) + 0.5 * chi2_ij * r / (1 + r)

    which the DP sums over a block's triangle (see ``_dp_segment_blocks``).
    ``g`` is *signed*: a null-like pair (``chi^2`` near 1) is dominated by
    the ``-0.5 * log(1 + r)`` intercept and contributes a small *negative*
    value, so widening a block across dead space costs -- the pressure that
    keeps distinct domains apart without any tuned threshold -- while a
    genuinely correlated pair (large ``chi^2``) contributes a large positive
    value. A strongly-connected hub is absorbed by the *accumulation* of its
    many individually-modest edges, not by up-weighting any single one.

    ``r`` (how much wider a domain's correlation distribution is than the
    null) is estimated from the data, not tuned: under the model
    ``E[chi^2] = 1 + r``, so ``r = mean(chi^2) - 1`` over the analyzed
    (``N >= min_pair_coverage``) pairs, floored above 0 (``_SCORE_MIN_SHAPE``).
    The null replicates are scored with the *same* ``r`` and transform (each
    is a chi-square draw under the null, so it scores ~0 or negative), and
    they calibrate the block-level penalty in ``_calc_null_penalty`` -- which
    is where the joint, cross-pair (transitivity) structure of the null
    enters, since each replicate is one full joint draw.

    Returns ``(scored_table, scored_null_tables)`` (copies with
    ``SCORE_COL`` added). Works with no null replicates too (``r`` needs only
    the observed data; the penalty then falls back to the analytical BIC).
    """
    with logger.debug.single_context("Scoring {} pairs", table.index.size):
        import numpy as np

        chi2 = table[CHI_SQUARE_COL].astype(float)
        analyzed = chi2[table[N_COL] >= min_pair_coverage]
        r = (float(np.nanmean(analyzed)) - 1.0) if len(analyzed.index) else 0.0
        r = max(r, _SCORE_MIN_SHAPE)
        logger.debug(
            "Estimated domain variance-inflation shape r={} from {} analyzed pair(s)",
            r,
            len(analyzed.index),
        )
        intercept = -np.log1p(r) / 2.0
        slope = r / (2.0 * (1.0 + r))

        def score(col: pd.Series) -> np.ndarray:
            c = col.to_numpy(dtype=float)
            return np.where(np.isfinite(c), intercept + slope * c, 0.0)
        
        scored_table = table.copy()
        scored_table[SCORE_COL] = score(chi2)
        scored_null_tables = []
        for nt in null_tables_per_replicate:
            scored_nt = nt.copy()
            scored_nt[SCORE_COL] = score(nt[CHI_SQUARE_COL])
            scored_null_tables.append(scored_nt)
        return scored_table, scored_null_tables


def _pair_band_row_cumsum(
    table: pd.DataFrame, total_end5: int, total_end3: int
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Build banded row-cumulative sums of the observable-pair count and
    the per-pair score (``SCORE_COL``) grids over positions
    [total_end5, total_end3], 0-indexed, so that the pair count and the
    summed score whose triangle lies within any interval [s, e] (both
    endpoints in [s, e]) can be read off via ``_triangle_sum_banded`` --
    see ``_calc_domains_by_dp_segmentation`` for how this drives the exact
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
    ``_calc_domains_by_dp_segmentation``) and to carry ``SCORE_COL`` (see
    ``_calc_pair_scores``): a pair filtered out there (including every pair
    touching a fully masked position) never reaches this function, so it is
    scattered into neither returned array and cannot contribute to any
    downstream count or block score.

    Returns ``(obs_row_cum, score_row_cum, n_positions, max_gap)``, where
    ``obs_row_cum[a, d]`` (for ``d`` in ``0..max_gap+1``) is the number of
    observable pairs ``(a, a+d')`` for ``d' = 0..d-1`` and
    ``score_row_cum[a, d]`` is the sum of their ``SCORE_COL`` values
    (``row_cum[a, 0] == 0``; ``row_cum[a, max_gap+1]`` is position ``a``'s
    full row total within the band).
    """
    import numpy as np

    n_positions = total_end3 - total_end5 + 1
    pos_a = table.index.get_level_values(POSITION_A).to_numpy() - total_end5
    pos_b = table.index.get_level_values(POSITION_B).to_numpy() - total_end5
    # Every pair must lie within [total_end5, total_end3]; otherwise the
    # scatter below would silently wrap a negative index or raise an
    # IndexError, corrupting the counts. Fail explicitly instead.
    if len(pos_a) and (
        pos_a.min() < 0 or pos_b.max() >= n_positions or (pos_b < pos_a).any()
    ):
        raise OutOfBoundsError(
            "All pairs must satisfy total_end5 ≤ a ≤ b ≤ total_end3 "
            f"(0 ≤ a ≤ b < {n_positions} after offsetting by "
            f"total_end5={total_end5})"
        )
    score = table[SCORE_COL].to_numpy(dtype=float)
    gaps = pos_b - pos_a
    max_gap = int(gaps.max()) if len(gaps) else 0

    obs_grid = np.zeros((n_positions, max_gap + 1), dtype=np.int64)
    score_grid = np.zeros((n_positions, max_gap + 1), dtype=float)
    np.add.at(obs_grid, (pos_a, gaps), 1)
    np.add.at(score_grid, (pos_a, gaps), score)

    obs_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
    score_row_cum = np.zeros((n_positions, max_gap + 2), dtype=float)
    obs_row_cum[:, 1:] = np.cumsum(obs_grid, axis=1)
    score_row_cum[:, 1:] = np.cumsum(score_grid, axis=1)
    return obs_row_cum, score_row_cum, n_positions, max_gap


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

    # Preserve the row-cumulative array's dtype: this helper sums the
    # integer observable-pair count grid (int64) *and* the float per-pair
    # score grid (see _pair_band_row_cumsum), so a hardcoded int64 output
    # would silently truncate the float block scores.
    result = np.zeros(e + 1, dtype=row_cum.dtype)
    if boundary > 0:
        s_far = np.arange(0, boundary)
        result[s_far] = (cum_full[boundary] - cum_full[s_far]) + near_total
    if len(near_a):
        result[boundary : e + 1] = near_suffix
    return result


def _dp_segment_blocks(
    obs_row_cum: np.ndarray,
    score_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    penalty: float,
) -> tuple[list[tuple[int, int]], float]:
    """Find the partition of [0, n_positions) into background and
    domain blocks that maximizes the total block score (the sum of each
    chosen block's summed per-pair score ``Sum s_ij`` over its triangle,
    minus ``penalty`` per block), via exact dynamic programming:
    ``dp[t] = max(dp[t-1], max_s dp[s] + gain(s, t-1) - penalty)``,
    where ``dp[t]`` is the best score using positions [0, t) and
    ``gain(s, t-1)`` is the summed score of a candidate block [s, t-1].
    Every interval length at every position is considered (a block may be
    far longer than ``max_gap`` even though the band restricts any single
    pair's span, so this loop is still O(n_positions^2), unchanged from a
    dense grid -- only the O(n_positions * max_gap) memory improves), so no
    window scale is ever chosen.

    The per-pair score ``s_ij`` (see ``_calc_pair_scores``) is signed: a
    genuinely correlated pair contributes a positive value and a null-like
    pair a small negative one, so extending a block across dead space costs
    -- the pressure that keeps distinct domains apart and that a strongly-
    connected hub's positive edges must overcome to be pulled in. Returns
    ``(blocks, objective)``: the chosen blocks as sorted, non-overlapping
    0-indexed (start, end) pairs, inclusive, and the maximized objective
    ``dp[n_positions]`` (the total summed score, already net of ``penalty``
    per block).
    """
    import numpy as np

    dp = np.zeros(n_positions + 1)
    backtrack = np.full(n_positions + 1, -1, dtype=np.int64)
    for t in range(1, n_positions + 1):
        e = t - 1
        s_vals = np.arange(t)
        n_se = _triangle_sum_banded(obs_row_cum, e, max_gap)
        gain = _triangle_sum_banded(score_row_cum, e, max_gap)
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
    return blocks, float(dp[n_positions])


def _scan_block_gains(
    table: pd.DataFrame,
    total_end5: int,
    total_end3: int,
    min_pair_coverage: int,
) -> np.ndarray:
    """Scan a scored ``table`` (carrying ``SCORE_COL``; see
    ``_calc_pair_scores``) for the summed scores of the positive-score
    blocks a penalty-free DP carves out, returning their gains.

    A penalty-free segmentation carves out every block whose pairs' scores
    sum to a positive value -- i.e. every stretch whose correlations exceed
    the null on balance. On the observed table these are the candidate
    domains; on a null replicate (scored leave-one-out, so its pairs are ~0
    on average) they are the chance blocks used to calibrate the per-block
    penalty (see ``_calc_null_penalty``)."""
    import numpy as np

    obs_table = table[table[N_COL] >= min_pair_coverage]
    obs_row_cum, score_row_cum, n_positions, max_gap = _pair_band_row_cumsum(
        obs_table, total_end5, total_end3
    )
    if int(obs_row_cum[:, max_gap + 1].sum()) == 0:
        return np.zeros(0, dtype=float)
    blocks, _ = _dp_segment_blocks(
        obs_row_cum, score_row_cum, n_positions, max_gap, 0.0
    )
    gains = np.empty(len(blocks), dtype=float)
    for i, (s, e) in enumerate(blocks):
        gains[i] = float(_triangle_sum_banded(score_row_cum, e, max_gap)[s])
    return gains


def _calibrate_penalty_from_null(
    real_gains: np.ndarray,
    null_gains: np.ndarray,
    n_replicates: int,
    domain_fdr: float,
) -> float:
    """Per-block penalty from a Benjamini-Hochberg-style FDR comparison of
    the real candidate-block gains to the pooled null block gains.

    For a threshold ``T``, the expected number of *false* domains is
    ``V(T) = (# null blocks with gain > T) / n_replicates`` and the number
    of *called* domains is ``R(T) = # real blocks with gain > T``, so the
    false-discovery proportion is ``V(T) / R(T)``. Scanning ``T`` over the
    real gains (each admits one more real block), take the lowest ``T``
    (most domains) whose FDR is still ``<= domain_fdr``, and return a
    penalty just below it so the DP admits exactly those blocks. Returns
    ``+inf`` if no threshold satisfies the FDR (call nothing)."""
    import numpy as np

    if real_gains.size == 0:
        return inf
    real_descending = np.sort(real_gains)[::-1]
    null_ascending = np.sort(null_gains)
    n_rep = max(int(n_replicates), 1)
    best_i = 0
    for i in range(1, real_descending.size + 1):
        threshold = real_descending[i - 1]
        n_null_above = null_gains.size - int(
            np.searchsorted(null_ascending, threshold, side="right")
        )
        fdr = (n_null_above / n_rep) / i
        if fdr <= domain_fdr:
            best_i = i  # largest i satisfying the FDR (BH step-up)
    if best_i == 0:
        # Even the single strongest real block fails the FDR: call none.
        return inf
    # Place the penalty between the weakest admitted block and the next
    # (rejected) one, so the DP keeps exactly the top best_i blocks.
    admitted = float(real_descending[best_i - 1])
    below = float(real_descending[best_i]) if best_i < real_descending.size else 0.0
    return (admitted + below) / 2.0


def _calc_null_penalty(
    table: pd.DataFrame,
    null_tables_per_replicate: list[pd.DataFrame],
    region,
    min_pair_coverage: int,
    domain_fdr: float,
) -> float:
    """Calibrate the DP's per-block penalty from a simulated null.

    Scan the observed ``table`` for its real candidate-domain gains (summed
    per-pair scores), then scan each of ``null_tables_per_replicate`` (one
    already-banded table per null replicate, scored with the *same*
    variance-inflation transform as the observed data and built by the caller
    from tables gathered while reading each tile's batches -- see
    ``_find_correlated_pairs_with_nulls`` and ``_calc_pair_scores``) for its
    chance block gains (``_scan_block_gains``). Each replicate is one joint
    draw under the null, so this is where the cross-pair (transitivity)
    dependence of the block-level null enters; a replicate's pairs are
    chi-square draws under the null, so its blocks score near 0. Return the
    FDR-controlled penalty comparing the two at ``domain_fdr``
    (``_calibrate_penalty_from_null``). Requires at least one null replicate,
    since the block-level null distribution (which captures the cross-pair
    transitivity structure) has no analytical substitute."""
    import numpy as np

    n_null_replicates = len(null_tables_per_replicate)
    if n_null_replicates < 1:
        raise OutOfBoundsError(
            f"n_null_replicates must be >= 1, but got {n_null_replicates}"
        )
    real_gains = _scan_block_gains(
        table, region.end5, region.end3, min_pair_coverage
    )
    pooled = []
    for r, null_table in enumerate(null_tables_per_replicate):
        gains = _scan_block_gains(
            null_table, region.end5, region.end3, min_pair_coverage
        )
        pooled.append(gains)
        logger.trace(
            "Null replicate {}/{} has {} spurious block(s)",
            r + 1,
            n_null_replicates,
            gains.size,
        )
    pooled_gains = np.concatenate(pooled) if pooled else np.zeros(0)
    penalty = _calibrate_penalty_from_null(
        real_gains, pooled_gains, n_null_replicates, domain_fdr
    )
    logger.debug(
        "Null-calibrated DP penalty {} from {} real candidate block(s) vs "
        "{} pooled null block(s) over {} replicate(s) (target FDR {})",
        penalty,
        real_gains.size,
        pooled_gains.size,
        n_null_replicates,
        domain_fdr,
    )
    return penalty


def _crossing_gain(
    score_row_cum: np.ndarray,
    max_gap: int,
    block_left: tuple[int, int],
    block_right: tuple[int, int],
) -> float:
    """Sum of ``SCORE_COL`` for exactly the pairs with one endpoint in
    ``block_left`` and the other in ``block_right`` -- excluding every
    pair that touches the gap between them (gap-to-gap or block-to-gap).

    Summed directly, one row at a time, from the banded row-cumulative
    score grid (``score_row_cum``; see ``_pair_band_row_cumsum``): for each
    position ``a`` in ``block_left``, the pairs ``(a, b)`` with ``b`` in
    ``block_right`` occupy a contiguous gap range in ``a``'s row, read off
    with two lookups. This is deliberately *not* computed as a difference
    of larger triangle sums (e.g. via inclusion-exclusion on
    ``_triangle_sum_banded``): when the true crossing gain is exactly zero
    (no pair of ``block_left`` and ``block_right`` is even within
    ``max_gap`` of each other), subtracting several ``O(block size)``
    triangle sums that are each individually nonzero can leave a
    floating-point residue of the wrong sign, spuriously triggering a join
    in ``_join_bridged_blocks`` on numerical noise alone. Summing only the
    pairs that actually exist avoids that cancellation entirely: with no
    such pairs, every row's contribution is exactly ``0.0``.

    This is what actually decides whether two DP-chosen blocks separated
    by dead space should be reported as one domain (see
    ``_join_bridged_blocks``): the crossing pairs are the only direct
    evidence connecting the two blocks, and testing them alone -- rather
    than summing in every pair that merely touches the gap -- avoids the
    gap's own (expectedly null) pairs, which can vastly outnumber the
    crossing pairs, diluting real evidence into a net negative."""
    import numpy as np

    s1, e1 = block_left
    s2, e2 = block_right
    a_vals = np.arange(s1, e1 + 1)
    gap_lo = s2 - a_vals
    gap_hi = np.minimum(e2 - a_vals, max_gap)
    in_band = gap_lo <= max_gap
    gap_lo_idx = np.clip(gap_lo, 0, max_gap + 1)
    gap_hi_idx = np.clip(gap_hi + 1, 0, max_gap + 1)
    row_hi = score_row_cum[a_vals, gap_hi_idx]
    row_lo = score_row_cum[a_vals, gap_lo_idx]
    return float(np.where(in_band, row_hi - row_lo, 0.0).sum())


def _join_bridged_blocks(
    table: pd.DataFrame,
    null_tables_per_replicate: list[pd.DataFrame],
    blocks: list[tuple[int, int]],
    total_end5: int,
    total_end3: int,
    min_pair_coverage: int,
    domain_fdr: float,
) -> list[tuple[int, int]]:
    """Join adjacent DP-chosen blocks whose *direct crossing pairs*
    (``_crossing_gain``) -- not the whole bridge between them, which is
    dominated by however much (expectedly null) dead space happens to sit
    in between -- clear the same null-calibrated FDR standard used to call
    a domain in the first place (``_calc_null_penalty``).

    The DP's own optimum is exact under the additive per-pair score (see
    ``_dp_segment_blocks``): whenever it leaves two blocks split, the
    *whole* bridge between them (their union's score minus each block's
    own) is necessarily <= -penalty, by construction of the DP recurrence
    -- so a post-hoc pass re-testing that same whole-bridge quantity could
    never find anything to join, and is not what this function does. A
    real gap can still sit between two blocks that are nonetheless
    directly connected by real long-range correlations (e.g. a helix with
    an unpaired linker or internal loop in between): the gap's own pairs
    are correctly null, yet summing them in with the real crossing pairs
    swamps the latter, masking a genuine relationship that a block -- a
    contiguous position range with no way to isolate "just the crossing
    pairs" from what merely surrounds them -- has no way to see on its
    own. This function is the direct, targeted test for that case.

    Requires at least one null replicate: the real crossing gains for
    every adjacent pair of DP blocks are compared against the crossing
    gains of the *same* block boundaries reapplied to each null
    replicate's scored table (one joint draw of the null each), and
    calibrated via the identical BH-style FDR logic used for the main
    domain-calling penalty (``_calibrate_penalty_from_null``) -- so
    joining two blocks is held to the same ``domain_fdr`` standard as
    calling one, with no separate tuned threshold. A single static pass
    over the original DP block boundaries suffices (no re-scoring after
    each join): each edge's crossing gain depends only on the two fixed
    original blocks it connects, so joins that fuse 3+ consecutive
    blocks fall out of chaining consecutive admitted edges together."""
    import numpy as np

    if len(blocks) < 2:
        return blocks
    n_null_replicates = len(null_tables_per_replicate)
    if n_null_replicates < 1:
        raise OutOfBoundsError(
            f"n_null_replicates must be >= 1, but got {n_null_replicates}"
        )

    def crossing_gains(scored_table: pd.DataFrame) -> np.ndarray:
        obs = scored_table[scored_table[N_COL] >= min_pair_coverage]
        _, score_row_cum, _, max_gap = _pair_band_row_cumsum(
            obs, total_end5, total_end3
        )
        return np.array(
            [
                _crossing_gain(score_row_cum, max_gap, blocks[i], blocks[i + 1])
                for i in range(len(blocks) - 1)
            ]
        )

    real_gains = crossing_gains(table)
    pooled = [crossing_gains(nt) for nt in null_tables_per_replicate]
    pooled_gains = np.concatenate(pooled) if pooled else np.zeros(0)
    threshold = _calibrate_penalty_from_null(
        real_gains, pooled_gains, n_null_replicates, domain_fdr
    )
    logger.debug(
        "Crossing-pair join threshold {} from {} adjacent block pair(s) vs "
        "{} pooled null crossing gain(s) over {} replicate(s) (target FDR {})",
        threshold,
        real_gains.size,
        pooled_gains.size,
        n_null_replicates,
        domain_fdr,
    )
    joined: list[tuple[int, int]] = [blocks[0]]
    for i in range(1, len(blocks)):
        if real_gains[i - 1] > threshold:
            s0, _ = joined[-1]
            _, e1 = blocks[i]
            joined[-1] = (s0, e1)
        else:
            joined.append(blocks[i])
    return joined


def _calc_domains_by_dp_segmentation(
    table: pd.DataFrame,
    null_tables_per_replicate: list[pd.DataFrame],
    total_end5: int,
    total_end3: int,
    penalty: float,
    domain_fdr: float,
    min_pair_coverage: int = 1000,
) -> list[tuple[int, int]]:
    """Call domains by exact dynamic-program segmentation of the
    per-pair correlation-score contact map into background/domain blocks
    (see ``_dp_segment_blocks``): partition [total_end5, total_end3] into
    consecutive intervals and keep as a domain each interval whose pairs'
    null-referenced scores (``SCORE_COL``; see ``_calc_pair_scores``) sum
    to more than ``penalty``, maximizing the total summed score net of one
    ``penalty`` per block.

    This scores each candidate block by the *summed null-referenced
    correlation excess* over its 2-D triangle, replacing an earlier
    significant-pair *density* (a rate ``k / n``). Density collapsed every
    pair to one significant/not bit and, being a ratio, was diluted by the
    empty span around a strongly-correlated but isolated hub position, so
    such a hub was excluded from its domain. The summed signed score fixes
    both: it keeps each pair's coverage-normalized correlation magnitude,
    and because a null-like pair contributes a small *negative* value while
    a hub's genuine edges contribute large positive ones, widening a block
    to absorb a hub across a short dead span is favorable, while widening it
    across a long genuine gap is not (that is what keeps distinct domains
    apart). No window scale is ever chosen: the DP considers every possible
    interval length at every position.

    Only pairs whose total coverage (``N``, the sum of the 2x2 confusion
    matrix) is at least ``min_pair_coverage`` are used. This is the *only*
    place masking is applied, and it is applied before anything else: a pair
    below the threshold -- including every pair touching a fully masked
    position, which has none above it -- is dropped from ``obs_table`` here
    and so never reaches ``_pair_band_row_cumsum``. It cannot contribute to
    any triangle count or block score downstream; it is not merely
    down-weighted, it simply does not exist for the rest of this function.

    ``penalty`` is calibrated from a simulated null (``_calc_null_penalty``).

    Represents the grid banded (see ``_pair_band_row_cumsum``), using
    O(L * max_gap) memory where L is the scanned region's length and
    max_gap is the largest separation of any observable pair -- tens of
    megabytes instead of the tens of gigabytes a dense O(L^2) grid would
    need. The DP itself is still O(L^2) time (a domain may be far longer
    than the band), unchanged from a dense grid.

    After the DP, adjacent blocks are joined when their *direct crossing
    pairs* alone clear the null-calibrated FDR bar (``_join_bridged_blocks``,
    at the same ``domain_fdr``) -- not a post-hoc re-test of the DP's own
    whole-bridge objective, which the DP's exactness already makes
    unwinnable, but a targeted test for a real gap that nonetheless has
    real long-range correlations crossing it.
    """
    if min_pair_coverage < 1:
        raise OutOfBoundsError(
            f"min_pair_coverage must be >= 1, but got {min_pair_coverage}"
        )
    if not 1 <= total_end5 <= total_end3:
        raise IncompatibleValuesError(
            "Must have 1 ≤ total_end5 ≤ total_end3, "
            f"but got total_end5={total_end5} and total_end3={total_end3}"
        )
    obs_table = table[table[N_COL] >= min_pair_coverage]
    with logger.debug.single_context(
        "Identifying domains from {} observable pair(s) (of {} in-band) "
        "by DP block-diagonal segmentation",
        len(obs_table.index),
        len(table.index),
    ):
        obs_row_cum, score_row_cum, n_positions, max_gap = _pair_band_row_cumsum(
            obs_table, total_end5, total_end3
        )
        blocks, _ = _dp_segment_blocks(
            obs_row_cum, score_row_cum, n_positions, max_gap, penalty
        )
        blocks = _join_bridged_blocks(
            table,
            null_tables_per_replicate,
            blocks,
            total_end5,
            total_end3,
            min_pair_coverage,
            domain_fdr,
        )
        domains = sorted((total_end5 + s, total_end5 + e) for s, e in blocks)
        logger.debug("Calculated domains: {}", domains)
        return domains


def _graph_pairs_and_domains(
    pairs: list[tuple[int, int]],
    scores: list[float],
    domains: list[tuple[int, int]],
    end5: int,
    end3: int,
    html_file: str | Path,
):
    import numpy as np
    from plotly import graph_objects as go

    fig = go.Figure()
    # Graph the domains as triangles. All domain traces share a legend
    # group so that clicking the single "Domains" legend entry toggles
    # every domain triangle on/off together (groupclick="togglegroup",
    # set below), acting as a switch for the whole layer.
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
                legendgroup="domains",
                showlegend=(i == 0),
                name="Domains",
                hoverinfo="text",
                text=f"Domain {a, b}",
                hovertemplate="%{text}<extra></extra>",
            )
        )
    # Plot the correlated pairs as points in the original pairs color
    # (#D55E00), with opacity set by their score on a square-root scale
    # anchored at zero (0% opacity at score 0, regardless of the minimum
    # score actually present, 100% opacity at the maximum score in the
    # data). A square-root scale sits between linear (too many low
    # scores crowded near 0% opacity) and log (too many points pulled up
    # toward 100%). Plotly ties a colorbar to marker.color + colorscale,
    # not marker.opacity, so the opacity ramp is expressed as a two-stop
    # colorscale in that same color (0% to 100% alpha) applied to
    # sqrt(score), which renders identically to varying opacity while
    # keeping the colorbar legend; the colorbar ticks are relabeled back
    # to score units since the underlying values are sqrt-scaled.
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    pairs_midpoints, pairs_distances = _calc_midpoints_distances(pos5s, pos3s)
    scores_arr = np.asarray(scores, dtype=float)
    finite_scores = scores_arr[np.isfinite(scores_arr)]
    pairs_color = "rgb(213,94,0)"  # #D55E00
    if finite_scores.size and finite_scores.max() > 0:
        max_score = float(finite_scores.max())
        sqrt_scores = np.sqrt(np.clip(scores_arr, 0.0, None))
        sqrt_max = float(np.sqrt(max_score))
        tick_vals = np.linspace(0.0, sqrt_max, num=6)
        marker = dict(
            color=sqrt_scores.tolist(),
            colorscale=[
                [0.0, "rgba(213,94,0,0)"],
                [1.0, "rgba(213,94,0,1)"],
            ],
            cmin=0.0,
            cmax=sqrt_max,
            showscale=True,
            colorbar=dict(
                title="Score",
                tickvals=tick_vals.tolist(),
                ticktext=[f"{v:.2g}" for v in (tick_vals**2)],
            ),
        )
    else:
        marker = dict(color=pairs_color)
    pairs_trace_index = len(domains)
    pairs_text = [
        f"Pair {pair}, score={score:.2f}"
        for pair, score in zip(pairs, scores, strict=True)
    ]
    fig.add_trace(
        go.Scatter(
            x=pairs_midpoints,
            y=pairs_distances,
            mode="markers",
            legendgroup="pairs",
            showlegend=True,
            marker=marker,
            name="Correlated pairs",
            hoverinfo="text",
            text=pairs_text,
            hovertemplate="%{text}<extra></extra>",
        )
    )
    # Clicking a legend entry toggles every trace in its group, turning
    # the legend into an on/off switch for the domains and the pairs.
    # Placed above the plot (horizontal orientation) so it doesn't
    # overlap the score colorbar, which sits to the right of the plot.
    fig.update_layout(
        legend=dict(
            groupclick="togglegroup",
            orientation="h",
            x=0.0,
            y=1.08,
            xanchor="left",
            yanchor="bottom",
        )
    )
    # Add a slider that hides pairs scoring below a chosen minimum, by
    # restyling the pairs trace's x/y/text/marker.color to the subset at
    # or above each threshold (cmin/cmax stay fixed at the full-data
    # range so the color/opacity scale doesn't shift as points drop out).
    if finite_scores.size:
        thresholds = np.square(np.linspace(0.0, np.sqrt(finite_scores.max()), 31))
        steps = []
        for t in thresholds:
            keep = scores_arr >= t
            step_args = {
                "x": [pairs_midpoints[keep].tolist()],
                "y": [pairs_distances[keep].tolist()],
                "text": [[txt for txt, k in zip(pairs_text, keep) if k]],
            }
            if "color" in marker and isinstance(marker["color"], list):
                step_args["marker.color"] = [
                    [c for c, k in zip(marker["color"], keep) if k]
                ]
            steps.append(
                dict(
                    method="restyle",
                    args=[step_args, [pairs_trace_index]],
                    label=f"{t:.2g}",
                )
            )
        fig.update_layout(
            sliders=[
                dict(
                    active=0,
                    currentvalue=dict(prefix="Min score: "),
                    pad=dict(t=50),
                    steps=steps,
                )
            ]
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


def _write_pairs_with_confusion(pos_table: pd.DataFrame, csv_file: str | Path):
    """Write the positive-score pairs (``SCORE_COL > 0``) to a CSV file,
    with their 2x2 confusion-matrix counts and null-referenced score
    (``SCORE_COL``; see ``_calc_pair_scores``). The four raw counts are
    kept because chi-square is recomputable from them (writing it too
    would be redundant); the score is kept anyway, since it also depends
    on the table-wide variance-inflation shape ``r`` and so cannot be
    recomputed from a single pair's counts alone."""
    import pandas as pd

    df = pd.DataFrame(
        {
            POSITION_A: pos_table.index.get_level_values(POSITION_A).to_numpy(),
            POSITION_B: pos_table.index.get_level_values(POSITION_B).to_numpy(),
        }
        | {col: pos_table[col].to_numpy() for col in CONFUSION_COLS}
        | {SCORE_COL: pos_table[SCORE_COL].to_numpy()}
    )
    df.to_csv(csv_file, index=False)


def _calc_cluster_domains(
    filter_dirs: list[Path],
    report_dir: Path,
    num_cpus: int,
    band_width: int,
    min_length: int,
    gap_mode: str,
    min_pair_coverage: int = 1000,
    n_null_replicates: int = 10,
    domain_fdr: float = 0.1,
    seed: int | None = None,
):
    """Calculate the cluster regions for all tiles of one reference.

    Domains are called by calibrating the DP's per-block penalty from
    ``n_null_replicates`` simulated per-position independence-null
    replicates at false-discovery rate ``domain_fdr`` (see
    ``_calc_null_penalty``); the block-level null distribution captures
    the cross-pair transitivity structure that has no analytical
    substitute, so at least one replicate is required."""
    # The FDR estimate has resolution 1 / n_null_replicates for the
    # strongest candidate domain (where the number of admitted real
    # domains is 1), so domain_fdr cannot be resolved below that: with
    # too few replicates, calling the strongest domain degenerates into
    # requiring zero null exceedances (much stricter than the target FDR).
    if domain_fdr > 0 and n_null_replicates < 1 / domain_fdr:
        logger.warning(
            "n_null_replicates ({}) is less than 1 / domain_fdr ({}); the "
            "target FDR {} cannot be resolved for the strongest domain, so "
            "domain calling will be overly conservative. Increase "
            "--n-null-replicates to at least {} or raise --domain-fdr.",
            n_null_replicates,
            ceil(1 / domain_fdr),
            domain_fdr,
            ceil(1 / domain_fdr),
        )
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
    # Build the banded per-pair chi-square table across all tiles, and (if
    # n_null_replicates > 0) every null replicate's table, from a single
    # read of each tile's batches (see _find_correlated_pairs_with_nulls;
    # this is what previously required re-reading each tile's batches once
    # per null replicate in addition to once for the observed data).
    n_tiles = len(datasets)
    # Draw one seed per tile from a stream initialized with the run's seed,
    # so each tile gets its own reproducible seed without offsetting the
    # base seed by an index (if seed is None, the stream is nondeterministic
    # but every drawn tile seed is still a concrete integer).
    tile_seed_stream = get_random_integer_generator(seed)
    tile_seeds = [next(tile_seed_stream) for _ in range(n_tiles)]
    with logger.debug.single_context(
        "Finding correlated pairs across {} tile(s) of {} "
        "(with {} null replicate(s) each)",
        n_tiles,
        ref,
        n_null_replicates,
    ):
        per_tile_results = dispatch(
            _find_correlated_pairs_with_nulls,
            num_cpus=num_cpus,
            pass_num_cpus=True,
            as_list=True,
            ordered=True,
            raise_on_error=True,
            args=list(zip(datasets, range(n_tiles), tile_seeds)),
            kwargs=dict(
                band_width=band_width,
                n_null_replicates=n_null_replicates,
                n_tiles=n_tiles,
            ),
        )
    per_tile_tables = [real_table for real_table, _ in per_tile_results]
    table = _build_banded_table(per_tile_tables, band_width)
    # Rebuild each null replicate's banded table across all tiles (ordered
    # dispatch above guarantees per_tile_results[t][1][r] is tile t's r-th
    # replicate).
    null_tables_per_replicate = [
        _build_banded_table(
            [null_tables[r] for _, null_tables in per_tile_results], band_width
        )
        for r in range(n_null_replicates)
    ]
    # Attach the per-pair block-scoring currency (SCORE_COL) to the observed
    # table and to each null replicate: the variance-inflation log-likelihood
    # ratio of each pair's chi-square under a wider-than-null (domain) model,
    # with the shape r estimated from the data (see _calc_pair_scores).
    table, null_tables_per_replicate = _calc_pair_scores(
        table, null_tables_per_replicate, min_pair_coverage
    )
    # Restrict to the same observable pairs the domain caller uses, so the
    # displayed/saved pairs and n_positive_pairs match what was actually
    # analyzed (a pair below min_pair_coverage is invisible to the DP and can
    # never lie inside a called domain). A pair is "positive" when its score
    # (the same null-referenced currency the DP itself sums) is > 0, i.e. it
    # looks more like the domain model than the null -- the one criterion for
    # what counts as a correlated pair anywhere in this pipeline.
    obs_table = table[table[N_COL] >= min_pair_coverage]
    pos_table = obs_table[obs_table[SCORE_COL] > 0]
    pairs = pos_table.index.to_list()
    pair_scores = pos_table[SCORE_COL].to_list()
    logger.debug("Began calibrating the domain-calling penalty for {}", ref)
    penalty = _calc_null_penalty(
        table,
        null_tables_per_replicate,
        region,
        min_pair_coverage=min_pair_coverage,
        domain_fdr=domain_fdr,
    )
    logger.debug("Ended calibrating the domain-calling penalty for {}", ref)
    # Find domains via DP block-diagonal segmentation of the per-pair score
    # contact map, using the null-calibrated penalty.
    domains = _calc_domains_by_dp_segmentation(
        table,
        null_tables_per_replicate,
        total_end5=region.end5,
        total_end3=region.end3,
        min_pair_coverage=min_pair_coverage,
        penalty=penalty,
        domain_fdr=domain_fdr,
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
    _write_pairs_with_confusion(pos_table, report_dir.joinpath(PAIRS_CSV))
    _write_pairs_to_csv(domains, report_dir.joinpath(DOMAINS_CSV))
    logger.debug(
        "Wrote {} pair(s) and {} domain(s) to {}", len(pairs), len(domains), report_dir
    )
    # Graph the correlated pairs and domains.
    html_file = report_dir.joinpath(PAIRS_DOMAINS_HTML)
    try:
        _graph_pairs_and_domains(
            pairs, pair_scores, domains, region.end5, region.end3, html_file
        )
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
    band_width: int,
    domain_fdr: float,
    n_null_replicates: int,
    min_pair_coverage: int,
    seed: int | None,
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
        domain_coords, n_positive_pairs = _calc_cluster_domains(
            tiled_dirs,
            report_dir=report_dir,
            band_width=band_width,
            domain_fdr=domain_fdr,
            n_null_replicates=n_null_replicates,
            min_pair_coverage=min_pair_coverage,
            seed=seed,
            min_length=min_cluster_length,
            gap_mode=gap_mode,
            num_cpus=num_cpus,
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
            domain_fdr=domain_fdr,
            n_null_replicates=n_null_replicates,
            seed=seed,
            min_cluster_length=min_cluster_length,
            gap_mode=gap_mode,
            # Results (store coordinates without the reference, which is
            # already recorded in the report).
            tile_coords=[(end5, end3) for _, end5, end3 in tile_coords],
            n_positive_pairs=n_positive_pairs,
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
