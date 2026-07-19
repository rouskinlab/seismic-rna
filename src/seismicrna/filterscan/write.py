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
    calc_bh_adjusted_pvals,
    calc_confusion_matrix,
    calc_confusion_chi_square,
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
CHI_SQUARE_COL = "Chi-Square"

# The four cells of each pair's 2x2 confusion matrix. These are the raw
# counts from which every derived quantity (N and chi-square) can be
# recomputed, so only these are written to pairs.csv.

# Floor for a candidate block's locally-fit variance-inflation shape r (see
# _block_score): keeps r strictly positive so the score stays defined even
# when a block carries essentially no correlation signal (mean chi-square at
# or below the null's 1), in which case the admission threshold declines to
# call any domain.
_SCORE_MIN_SHAPE = 1e-9

NEITHER_COL = "Neither Mutated"
ONLY_A_COL = "Only A Mutated"
ONLY_B_COL = "Only B Mutated"
BOTH_COL = "Both Mutated"
CONFUSION_COLS = (NEITHER_COL, ONLY_A_COL, ONLY_B_COL, BOTH_COL)


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
    table: pd.DataFrame, min_pair_coverage: int, min_expect_both: float
) -> pd.Series:
    """Pairs usable for chi-square-based analysis: enough joint coverage
    (``min_pair_coverage``) and a large enough independence-expected
    both-mutated count (``min_expect_both``) for the chi-square
    approximation to be reliable."""
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
            pos_index,
            covered_pooled,
            mutated_pooled,
            None,
            min_gap=min_gap,
            max_gap=max_gap,
        )
        real_table = _confusion_to_table(n, a, b, ab, dataset=dataset, write_csv=True)
        logger.debug("Computed the observed confusion matrix for {}", dataset)
        return real_table


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


def _block_score(n: np.ndarray, chi2_sum: np.ndarray) -> np.ndarray:
    """Score a candidate block (or array of candidate blocks) of ``n``
    pairs whose chi-squares sum to ``chi2_sum``, by its own locally-fit
    log-likelihood ratio against the null, BIC-corrected for the one free
    parameter it fits.

    Each candidate block fits its *own* variance-inflation shape from just
    its own pairs, so that a strong domain elsewhere in the scanned region
    can't inflate the score of an unrelated bridge::

        r_hat     = max(chi2_sum / n - 1, _SCORE_MIN_SHAPE)
        intercept = -0.5 * log1p(r_hat)
        slope     = r_hat / (2 * (1 + r_hat))
        llr       = n * intercept + slope * chi2_sum

    ``llr`` is the block's own best-fit log-likelihood ratio: since ``r_hat``
    is the MLE for that block alone, ``llr`` is *always* the largest value
    achievable for that block's data, which would let any small subset of
    pairs with an above-average chi-square by pure chance look significant
    against itself. Charging a per-parameter BIC cost
    (``0.5 * log(n)``, for the one free parameter ``r_hat``) corrects that:
    a truly null block (``chi2`` averaging ~1, ``r_hat`` at its floor) then
    scores negative for any ``n > 1``, so extending a block into dead space
    still costs, exactly as the old global-``r`` score's negative intercept
    did -- but now the pressure comes from each block's own fit, not from a
    value borrowed from elsewhere in the scan.

    Returns ``gain = llr - log(n) / 2``, the *raw* (uncharged beyond BIC)
    per-block gain that ``_dp_segment_blocks`` sums directly, with no further
    per-block penalty (see ``_dp_segment_blocks`` for why: unlike a flat
    subtracted penalty, this makes choosing more blocks free provided each
    one clears its own admission bar, so two adjacent real domains are never
    cheaper to merge than to keep separate)."""
    import numpy as np

    with np.errstate(divide="ignore", invalid="ignore"):
        r_hat = np.maximum(chi2_sum / n - 1.0, _SCORE_MIN_SHAPE)
    intercept = -np.log1p(r_hat) / 2.0
    slope = r_hat / (2.0 * (1.0 + r_hat))
    llr = n * intercept + slope * chi2_sum
    bic = np.log(np.maximum(n, 1)) / 2.0
    return llr - bic


def _pair_band_row_cumsum(
    table: pd.DataFrame, total_end5: int, total_end3: int, value_col: str
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Build banded row-cumulative sums of the observable-pair count and
    of ``value_col`` (e.g. ``CHI_SQUARE_COL``, for domain calling) over
    positions [total_end5, total_end3], 0-indexed, so that the pair count
    and the summed value whose triangle lies within any interval [s, e]
    (both endpoints in [s, e]) can be read off via ``_triangle_sum_banded``
    -- see ``_calc_domains_by_dp_segmentation`` for how this drives the
    exact dynamic-program segmentation.

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
    ``_calc_domains_by_dp_segmentation``) and to carry ``value_col``: a pair
    filtered out there (including every pair touching a fully masked
    position) never reaches this function, so it is scattered into neither
    returned array and cannot contribute to any downstream count or block
    score.

    Returns ``(obs_row_cum, value_row_cum, n_positions, max_gap)``, where
    ``obs_row_cum[a, d]`` (for ``d`` in ``0..max_gap+1``) is the number of
    observable pairs ``(a, a+d')`` for ``d' = 0..d-1`` and
    ``value_row_cum[a, d]`` is the sum of their ``value_col`` values
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
    value = table[value_col].to_numpy(dtype=float)
    gaps = pos_b - pos_a
    max_gap = int(gaps.max()) if len(gaps) else 0

    obs_grid = np.zeros((n_positions, max_gap + 1), dtype=np.int64)
    value_grid = np.zeros((n_positions, max_gap + 1), dtype=float)
    np.add.at(obs_grid, (pos_a, gaps), 1)
    np.add.at(value_grid, (pos_a, gaps), value)

    obs_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
    value_row_cum = np.zeros((n_positions, max_gap + 2), dtype=float)
    obs_row_cum[:, 1:] = np.cumsum(obs_grid, axis=1)
    value_row_cum[:, 1:] = np.cumsum(value_grid, axis=1)
    return obs_row_cum, value_row_cum, n_positions, max_gap


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


def _calc_block_pvalue_cutoff(
    obs_row_cum: np.ndarray,
    chi2_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    detect_fdr: float,
) -> float:
    """Benjamini-Hochberg p-value cutoff over *every* candidate block the DP
    could ever consider, computed analytically with no null replicates.

    Each candidate block's own exact null p-value is ``chi2.sf(chi2_sum,
    df=n)``: valid with no simulation, since different pairs' chi-square
    values are independent under this null (confirmed empirically this
    session, both on average and conditional on other pairs looking
    elevated) and each pair's own marginal chi-square distribution is exact
    ``chi2(1)`` regardless of its correlation with other pairs -- so a
    block's chi-square sum is a sum of independent (though not identical)
    terms, exactly ``chi2(n)``-distributed to the same precision established
    for ``_block_score``.

    This complements, rather than replaces, ``_block_score``'s admission bar:
    a single "lucky" pair or small fragment can score positively enough under
    the BIC-corrected LLR to look like its own domain (verified this session
    against realistic per-pair noise), but is very unlikely to survive
    Benjamini-Hochberg correction against the full candidate pool an actual
    scan considers, since that pool is dominated by null-like candidates.
    Enumerates the same O(L^2) candidate space ``_dp_segment_blocks`` visits
    (one pass, before the DP itself runs), so this changes the constant
    factor of a scan, not its asymptotic cost.

    Returns the largest raw p-value whose BH-adjusted q-value is
    ``<= detect_fdr``, or ``-1.0`` if none survive (so that no p-value,
    always ``>= 0``, can ever be admitted -- i.e. "admit nothing")."""
    import numpy as np
    from scipy.stats import chi2 as chi2dist

    all_n = []
    all_chi2 = []
    for t in range(1, n_positions + 1):
        e = t - 1
        n_se = _triangle_sum_banded(obs_row_cum, e, max_gap)
        chi2_se = _triangle_sum_banded(chi2_row_cum, e, max_gap)
        mask = n_se > 0
        all_n.append(n_se[mask])
        all_chi2.append(chi2_se[mask])
    n_all = np.concatenate(all_n) if all_n else np.zeros(0)
    chi2_all = np.concatenate(all_chi2) if all_chi2 else np.zeros(0)
    if n_all.size == 0:
        return -1.0
    pvals = chi2dist.sf(chi2_all, df=n_all)
    qvals = calc_bh_adjusted_pvals(pvals)
    significant = pvals[qvals <= detect_fdr]
    if significant.size == 0:
        return -1.0
    return float(significant.max())


def _dp_segment_blocks(
    obs_row_cum: np.ndarray,
    chi2_row_cum: np.ndarray,
    n_positions: int,
    max_gap: int,
    p_cutoff: float,
) -> tuple[list[tuple[int, int]], float]:
    """Find the partition of [0, n_positions) into background and domain
    blocks that maximizes the total admitted block score, via exact dynamic
    programming: ``dp[t] = max(dp[t-1], max_s dp[s] + gain(s, t-1))`` over
    every candidate block ``(s, t-1)`` whose own ``gain`` (``_block_score``,
    fit locally from that block's own pair count and chi-square sum) is
    ``>= 0`` *and* whose own exact chi-square p-value clears the
    Benjamini-Hochberg cutoff ``p_cutoff`` (``_calc_block_pvalue_cutoff``);
    ``dp[t]`` is the best score using positions [0, t), and ``gain(s, t-1)``
    is the block's own locally-fit gain -- summed *raw*, with **no
    per-block subtraction**.

    ``gain >= 0`` needs no calibration: ``_block_score``'s BIC charge
    already makes a truly null block score negative for any ``n > 1``, so
    zero is the analytically-correct bar for "this block's own fit clears
    its own complexity cost".

    This is the key structural property that fixes over-merging: a flat
    subtracted penalty makes choosing one more block always cost a fixed
    fee, so two adjacent real domains can be cheaper to merge into one block
    than to keep separate whenever the dilution from a null bridge between
    them is milder than that fee. Gating on a fixed bar instead of
    subtracting a penalty means picking more blocks is never penalized --
    each one just has to independently clear the bar -- so two strong
    domains keep their own high scores intact when kept separate, while a
    diluted merged block (dragged down by an intervening null bridge) can
    only lose out.

    The p-value gate is a second, independent admission bar, not a
    replacement for the gain bar: a single "lucky" pair or small fragment
    can score positively enough under the BIC-corrected LLR to look like its
    own domain (verified against realistic per-pair noise), but is unlikely
    to survive Benjamini-Hochberg correction against the full candidate pool
    an actual scan considers. Eligibility gates never change which of two
    *already-eligible* alternatives the DP prefers -- that's still decided
    purely by comparing their ``_block_score`` values -- so this closes the
    small-fragment gap without touching the merge/split behavior above.

    Every interval length at every position is considered (a block may be
    far longer than ``max_gap`` even though the band restricts any single
    pair's span, so this loop is still O(n_positions^2), unchanged from a
    dense grid -- only the O(n_positions * max_gap) memory improves), so no
    window scale is ever chosen.

    Returns ``(blocks, objective)``: the chosen blocks as sorted,
    non-overlapping 0-indexed (start, end) pairs, inclusive, and the
    maximized objective ``dp[n_positions]`` (the total summed raw gain of
    the admitted blocks).
    """
    import numpy as np
    from scipy.stats import chi2 as chi2dist

    dp = np.zeros(n_positions + 1)
    backtrack = np.full(n_positions + 1, -1, dtype=np.int64)
    for t in range(1, n_positions + 1):
        e = t - 1
        s_vals = np.arange(t)
        n_se = _triangle_sum_banded(obs_row_cum, e, max_gap)
        chi2_se = _triangle_sum_banded(chi2_row_cum, e, max_gap)
        gain = _block_score(n_se, chi2_se)
        with np.errstate(divide="ignore", invalid="ignore"):
            pvals = chi2dist.sf(chi2_se, df=np.maximum(n_se, 1))
        eligible = (n_se > 0) & (gain >= 0.0) & (pvals <= p_cutoff)
        candidate_vals = np.where(eligible, dp[s_vals] + gain, -np.inf)
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


def _cut_crossing_scores(
    obs_table: pd.DataFrame, total_end5: int, total_end3: int
) -> tuple[np.ndarray, np.ndarray]:
    """For every cut -- the boundary between position ``m`` and ``m + 1``,
    for ``m`` in ``[total_end5, total_end3)`` -- return the raw
    ``(n_cross, chi2_cross)`` of every observable pair ``(a, b)`` with
    ``a <= m < b``: the pairs whose two endpoints straddle that specific
    cut.

    Computed as an ``O(n_pairs + L)`` difference-array sweep: each pair
    ``(a, b)`` straddles every cut ``m = a .. b - 1``, so it contributes
    ``+1``/``+chi_square`` at cut ``a`` and ``-1``/``-chi_square`` at cut
    ``b``, before a single cumulative sum -- not via inclusion-exclusion on
    ``_triangle_sum_banded`` triangle sums, for the same floating-point-
    cancellation reason this function's predecessor, ``_crossing_gain``,
    avoided that: a cut with truly zero crossing pairs must read back as
    exactly zero, not a rounding residue of the wrong sign.

    This is what ``_merge_connected_blocks`` uses to decide whether two
    DP-chosen blocks separated by dead space are still directly connected:
    unlike a block-to-block-only crossing test, every cut this function
    scores sees *every* pair that spans it, including one that reaches from
    deep inside one block, past an intervening block, into another --
    exactly the long-range evidence a block (which can only test its own
    contiguous span) cannot see on its own.

    Returns ``(n_cross, chi2_cross)``, each of length
    ``total_end3 - total_end5`` (one entry per cut, 0-indexed by
    ``m - total_end5`` for ``m`` in ``[total_end5, total_end3)``)."""
    import numpy as np

    n_positions = total_end3 - total_end5 + 1
    n_cuts = max(n_positions - 1, 0)
    # An empty table's index has no dtype to infer integer positions from
    # (it defaults to float), which np.add.at rejects as an indexer -- cast
    # explicitly so the empty case (as well as the normal one) works.
    pos_a = (
        obs_table.index.get_level_values(POSITION_A).to_numpy(dtype=np.int64)
        - total_end5
    )
    pos_b = (
        obs_table.index.get_level_values(POSITION_B).to_numpy(dtype=np.int64)
        - total_end5
    )
    value = obs_table[CHI_SQUARE_COL].to_numpy(dtype=float)
    diff_n = np.zeros(n_positions + 1)
    diff_chi2 = np.zeros(n_positions + 1)
    np.add.at(diff_n, pos_a, 1)
    np.add.at(diff_n, pos_b, -1)
    np.add.at(diff_chi2, pos_a, value)
    np.add.at(diff_chi2, pos_b, -value)
    n_cross = np.cumsum(diff_n)[:n_cuts]
    chi2_cross = np.cumsum(diff_chi2)[:n_cuts]
    return n_cross, chi2_cross


def _merge_connected_blocks(
    blocks: list[tuple[int, int]],
    n_cross: np.ndarray,
    chi2_cross: np.ndarray,
    merge_fdr: float,
) -> list[tuple[int, int]]:
    """Merge adjacent DP-chosen blocks whose intervening gap is, at *every*
    cut within it, directly connected by real long-range correlations --
    not a re-test of the DP's own whole-bridge objective (which the DP's
    exactness already makes unwinnable: whenever it leaves two blocks
    split, no single block spanning their whole union could have cleared
    the admission bar, or the DP would already have chosen it), but a
    targeted test for a real gap that nonetheless has real correlations
    crossing it (e.g. a helix with an unpaired linker or internal loop in
    the middle).

    A cut is connected iff its crossing pairs (``_cut_crossing_scores``)
    clear the same two-gate admission standard used to call a domain in
    the first place: ``_block_score(n_cross, chi2_cross) >= 0`` (effect
    size) *and* an exact chi-square p-value (``chi2.sf(chi2_cross,
    df=n_cross)``) that survives Benjamini-Hochberg correction across every
    cut in the scanned region, at ``merge_fdr``. Both gates are required:
    a cut's crossing pair count can be in the thousands (every pair that
    reaches across it, not just those touching the two specific blocks
    being tested), so BH-significance alone flags biologically-trivial
    elevations (mean chi-square barely above 1) as "significant" --
    verified this session to over-merge a real domain deep into adjacent
    structureless dead space. The effect-size floor is what keeps a
    genuinely null cut (chi-square averaging ~1) below zero and holds the
    boundary there.

    Because whether two blocks are joined depends only on the cuts *within
    their own gap* -- fixed values that merging never changes -- marking
    each of the ``len(blocks) - 1`` inter-block gaps connected/not and then
    grouping maximal runs of connected gaps is equivalent to a transitive,
    order-independent merge: sweeping left-to-right, right-to-left, or via
    union-find all yield the identical partition."""
    import numpy as np
    from scipy.stats import chi2 as chi2dist

    if len(blocks) < 2:
        return blocks
    valid = n_cross > 0
    gains = np.where(valid, _block_score(np.maximum(n_cross, 1), chi2_cross), -1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pvals = np.where(valid, chi2dist.sf(chi2_cross, df=np.maximum(n_cross, 1)), 1.0)
    qvals = np.ones_like(pvals)
    if valid.any():
        qvals[valid] = calc_bh_adjusted_pvals(pvals[valid])
    connected_cut = valid & (gains >= 0.0) & (qvals <= merge_fdr)
    gap_connected = [
        bool(connected_cut[blocks[i][1] : blocks[i + 1][0]].all())
        for i in range(len(blocks) - 1)
    ]
    logger.debug(
        "Cut-crossing merge: {}/{} adjacent block pair(s) fully connected "
        "(target FDR {})",
        sum(gap_connected),
        len(gap_connected),
        merge_fdr,
    )
    merged: list[tuple[int, int]] = [blocks[0]]
    for i, connected in enumerate(gap_connected):
        if connected:
            s0, _ = merged[-1]
            _, e1 = blocks[i + 1]
            merged[-1] = (s0, e1)
        else:
            merged.append(blocks[i + 1])
    return merged


def _calc_domains_by_dp_segmentation(
    table: pd.DataFrame,
    total_end5: int,
    total_end3: int,
    detect_fdr: float,
    merge_fdr: float,
    min_pair_coverage: int,
    min_expect_both: float,
) -> list[tuple[int, int]]:
    """Call domains by exact dynamic-program segmentation of the per-pair
    chi-square contact map into background/domain blocks (see
    ``_dp_segment_blocks``): partition [total_end5, total_end3] into
    consecutive intervals and keep as a domain each interval whose own
    locally-fit log-likelihood ratio (``_block_score``: its own variance-
    inflation shape ``r``, fit from just that interval's pairs, BIC-charged
    for that one free parameter) is ``>= 0`` *and* whose own exact
    chi-square p-value clears a Benjamini-Hochberg cutoff calibrated over
    every candidate this scan considers (``_calc_block_pvalue_cutoff``, no
    null replicates needed for either gate), maximizing the total of the
    admitted blocks' raw gains (no per-block subtraction -- see
    ``_dp_segment_blocks`` for why admitting more blocks is never itself
    penalized).

    Each candidate block fits its *own* variance-inflation shape rather than
    sharing one estimated across the whole scanned region: a strong domain
    elsewhere in the region can no longer inflate the apparent significance
    of an unrelated, genuinely null bridge, and merging two real domains
    across such a bridge can only dilute their combined score (never make it
    cheaper than keeping them separate, since there is no separate-block
    fee to avoid by merging). No window scale is ever chosen: the DP
    considers every possible interval length at every position.

    Only pairs that pass ``_analyzed_pairs_mask`` (enough joint coverage,
    ``min_pair_coverage``, and a large enough independence-expected
    both-mutated count, ``min_expect_both``, for the chi-square
    approximation to be reliable) are used. This is the *only* place
    masking is applied, and it is applied before anything else: a pair that
    fails either check -- including every pair touching a fully masked
    position, which has none above it -- is dropped from ``obs_table`` here
    and so never reaches ``_pair_band_row_cumsum``. It cannot contribute to
    any triangle count or block score downstream; it is not merely
    down-weighted, it simply does not exist for the rest of this function.

    Represents the grid banded (see ``_pair_band_row_cumsum``), using
    O(L * max_gap) memory where L is the scanned region's length and
    max_gap is the largest separation of any observable pair -- tens of
    megabytes instead of the tens of gigabytes a dense O(L^2) grid would
    need. The DP itself is still O(L^2) time (a domain may be far longer
    than the band), unchanged from a dense grid.

    After the DP, adjacent blocks are merged when every cut in the gap
    between them is directly connected by crossing pairs alone
    (``_merge_connected_blocks``, at ``merge_fdr``) -- not a
    post-hoc re-test of the DP's own whole-bridge objective, which the DP's
    exactness already makes unwinnable, but a targeted test for a real gap
    that nonetheless has real long-range correlations crossing it (a single
    per-block score, fit from only that block's own contained pairs, cannot
    see this: it is always cheaper to carve out the densest sub-interval of
    a heterogeneous domain than to fit the whole thing, so the fix for
    fragmentation has to live at the cuts, not in the block score).
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
    obs_table = table[_analyzed_pairs_mask(table, min_pair_coverage, min_expect_both)]
    with logger.debug.single_context(
        "Identifying domains from {} observable pair(s) (of {} in-band) "
        "by DP block-diagonal segmentation",
        len(obs_table.index),
        len(table.index),
    ):
        obs_row_cum, chi2_row_cum, n_positions, max_gap = _pair_band_row_cumsum(
            obs_table, total_end5, total_end3, value_col=CHI_SQUARE_COL
        )
        p_cutoff = _calc_block_pvalue_cutoff(
            obs_row_cum, chi2_row_cum, n_positions, max_gap, detect_fdr
        )
        blocks, _ = _dp_segment_blocks(
            obs_row_cum, chi2_row_cum, n_positions, max_gap, p_cutoff
        )
        n_cross, chi2_cross = _cut_crossing_scores(obs_table, total_end5, total_end3)
        blocks = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr)
        domains = sorted((total_end5 + s, total_end5 + e) for s, e in blocks)
        logger.debug("Calculated domains: {}", domains)
        return domains


def _graph_pairs_and_domains(
    pairs: list[tuple[int, int]],
    chi_squares: list[float],
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
    # (#D55E00), with opacity set by their chi-square on a square-root scale
    # anchored at zero (0% opacity at chi-square 0, regardless of the minimum
    # chi-square actually present, 100% opacity at the maximum chi-square in
    # the data). A square-root scale sits between linear (too many low
    # chi-squares crowded near 0% opacity) and log (too many points pulled up
    # toward 100%). Plotly ties a colorbar to marker.color + colorscale,
    # not marker.opacity, so the opacity ramp is expressed as a two-stop
    # colorscale in that same color (0% to 100% alpha) applied to
    # sqrt(chi-square), which renders identically to varying opacity while
    # keeping the colorbar legend; the colorbar ticks are relabeled back
    # to chi-square units since the underlying values are sqrt-scaled.
    pos5s, pos3s = _tuples_to_ends_arrays(pairs)
    pairs_midpoints, pairs_distances = _calc_midpoints_distances(pos5s, pos3s)
    chi2_arr = np.asarray(chi_squares, dtype=float)
    finite_chi2 = chi2_arr[np.isfinite(chi2_arr)]
    pairs_color = "rgb(213,94,0)"  # #D55E00
    if finite_chi2.size and finite_chi2.max() > 0:
        max_chi2 = float(finite_chi2.max())
        sqrt_chi2 = np.sqrt(np.clip(chi2_arr, 0.0, None))
        sqrt_max = float(np.sqrt(max_chi2))
        tick_vals = np.linspace(0.0, sqrt_max, num=6)
        marker = dict(
            color=sqrt_chi2.tolist(),
            colorscale=[[0.0, "rgba(213,94,0,0)"], [1.0, "rgba(213,94,0,1)"]],
            cmin=0.0,
            cmax=sqrt_max,
            showscale=True,
            colorbar=dict(
                title="Chi-Square",
                tickvals=tick_vals.tolist(),
                ticktext=[f"{v:.2g}" for v in (tick_vals**2)],
            ),
        )
    else:
        marker = dict(color=pairs_color)
    pairs_trace_index = len(domains)
    pairs_text = [
        f"Pair {pair}, chi-square={chi2:.2f}"
        for pair, chi2 in zip(pairs, chi_squares, strict=True)
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
    # overlap the chi-square colorbar, which sits to the right of the plot.
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
    # Add a slider that hides pairs with chi-square below a chosen minimum,
    # by restyling the pairs trace's x/y/text/marker.color to the subset at
    # or above each threshold (cmin/cmax stay fixed at the full-data
    # range so the color/opacity scale doesn't shift as points drop out).
    if finite_chi2.size:
        thresholds = np.square(np.linspace(0.0, np.sqrt(finite_chi2.max()), 31))
        steps = []
        for t in thresholds:
            keep = chi2_arr >= t
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
                    currentvalue=dict(prefix="Min chi-square: "),
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
    """Write the given pairs' 2x2 confusion-matrix counts to a CSV file.
    Chi-square (and every other derived quantity) is recomputable from
    the four raw counts, so only these are written."""
    import pandas as pd

    df = pd.DataFrame(
        {
            POSITION_A: pos_table.index.get_level_values(POSITION_A).to_numpy(),
            POSITION_B: pos_table.index.get_level_values(POSITION_B).to_numpy(),
        }
        | {col: pos_table[col].to_numpy() for col in CONFUSION_COLS}
    )
    df.to_csv(csv_file, index=False)


def _calc_cluster_domains(
    filter_dirs: list[Path],
    report_dir: Path,
    num_cpus: int,
    band_width: int,
    min_length: int,
    gap_mode: str,
    min_pair_coverage: int,
    min_expect_both: float,
    detect_fdr: float,
    merge_fdr: float,
):
    """Calculate the cluster regions for all tiles of one reference.

    Domains are called by an exact chi-square p-value per candidate block
    (``_calc_block_pvalue_cutoff``), Benjamini-Hochberg-adjusted at false-
    discovery rate ``detect_fdr``, then merged across gaps at false-
    discovery rate ``merge_fdr`` -- no null-replicate simulation needed."""
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
    # Build the banded per-pair chi-square table across all tiles, from a
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
    # Restrict to the same observable pairs the domain caller uses, so the
    # displayed/saved pairs and n_positive_pairs match what was actually
    # analyzed (a pair failing _analyzed_pairs_mask is invisible to the DP
    # and can never lie inside a called domain). A pair is "positive" when
    # its chi-square exceeds the chi^2(1) null's expectation of 1 -- the one
    # criterion for what counts as a correlated pair anywhere in this
    # pipeline's reporting.
    obs_table = table[_analyzed_pairs_mask(table, min_pair_coverage, min_expect_both)]
    pos_table = obs_table[obs_table[CHI_SQUARE_COL] > 1]
    pairs = pos_table.index.to_list()
    pair_chi2 = pos_table[CHI_SQUARE_COL].to_list()
    # Find domains via DP block-diagonal segmentation of the per-pair
    # chi-square contact map.
    domains = _calc_domains_by_dp_segmentation(
        table,
        total_end5=region.end5,
        total_end3=region.end3,
        min_pair_coverage=min_pair_coverage,
        min_expect_both=min_expect_both,
        detect_fdr=detect_fdr,
        merge_fdr=merge_fdr,
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
            pairs, pair_chi2, domains, region.end5, region.end3, html_file
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
    detect_fdr: float,
    merge_fdr: float,
    min_pair_coverage: int,
    min_expect_both: float,
    min_domain_length: int,
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
            detect_fdr=detect_fdr,
            merge_fdr=merge_fdr,
            min_pair_coverage=min_pair_coverage,
            min_expect_both=min_expect_both,
            min_length=min_domain_length,
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
            min_expect_both=min_expect_both,
            detect_fdr=detect_fdr,
            merge_fdr=merge_fdr,
            min_domain_length=min_domain_length,
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
