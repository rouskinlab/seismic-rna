import tempfile
import unittest as ut
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.batch.confusion import (
    POSITION_A,
    POSITION_B,
    calc_bh_adjusted_pvals,
)
from seismicrna.cluster.data import ClusterMutsDataset, get_clust_params
from seismicrna.cluster.main import run as run_cluster
from seismicrna.core.error import IncompatibleValuesError, OutOfBoundsError
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.report import DomainCoordsF
from seismicrna.filter.main import run as run_filter
from seismicrna.main import cli as seismic_cli
from seismicrna.core.rna.convert import db_to_ct
from seismicrna.core.seq.fasta import write_fasta
from seismicrna.core.seq.region import FULL_NAME
from seismicrna.core.seq.xna import DNA
from seismicrna.filterscan.main import run as run_filterscan
from seismicrna.filterscan.report import FilterScanReport
from seismicrna.filterscan.write import (
    N_COL,
    CONFUSION_COLS,
    NEITHER_COL,
    ONLY_A_COL,
    ONLY_B_COL,
    BOTH_COL,
    BRIDGE_COL,
    P_VALUE_COL,
    ADJ_P_VALUE_COL,
    FOLD_CHANGE_COL,
    _calc_tiles,
    _build_banded_table,
    _write_pairs_with_confusion,
    _bridge_mask,
    _block_score,
    _pair_band_row_cumsum,
    _triangle_sum_banded,
    _calc_block_pvalue_cutoff,
    _dp_segment_blocks,
    _cut_crossing_scores,
    _merge_connected_blocks,
    _extend_domains_by_bridges,
    _estimate_null_bridge_rate,
    _calc_domains_by_dp_segmentation,
    _filter_domains_length,
    _split_gap_evenly,
    _fill_domains_into_gaps,
    _widen_domains_into_gaps,
    _label_widened,
)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.idmut import run as sim_idmut


def _make_table(rows: list[tuple[int, int, float]]):
    """Build a per-pair table like ``_confusion_to_table`` returns, from an
    iterable of (pos_a, pos_b, n).

    The four confusion cells are synthesized as a simple, fixed split of
    ``n`` purely so ``_analyzed_pairs_mask`` (which needs them to compute the
    independence-expected both-mutated count) has columns to read; tests
    that care about the exact confusion counts (e.g. ``min_expect_both``
    filtering) build their own table."""
    import pandas as pd

    if not rows:
        pos_a, pos_b, ns = (), (), ()
    else:
        pos_a, pos_b, ns = zip(*rows)
    index = pd.MultiIndex.from_arrays([pos_a, pos_b], names=[POSITION_A, POSITION_B])
    only_a = [n / 2 for n in ns]
    only_b = [n / 2 for n in ns]
    both = [n / 4 for n in ns]
    neither = [n - (a + b - ab) for n, a, b, ab in zip(ns, only_a, only_b, both)]
    return pd.DataFrame(
        {
            N_COL: ns,
            NEITHER_COL: neither,
            ONLY_A_COL: only_a,
            ONLY_B_COL: only_b,
            BOTH_COL: both,
        },
        index=index,
    )


class TestCalcTiles(ut.TestCase):
    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR)

    def tearDown(self):
        set_config(*self._config)

    def test_tile_min_overlap_25(self):
        result = _calc_tiles(41, 145, 60, 0.25)
        expect = [(41, 100), (86, 145)]
        self.assertEqual(result, expect)

    def test_tile_min_overlap_75(self):
        result = _calc_tiles(41, 145, 60, 0.75)
        expect = [(41, 100), (56, 115), (71, 130), (86, 145)]
        self.assertEqual(result, expect)

    def test_region_not_divisible(self):
        result = _calc_tiles(41, 170, 60, 0.25)
        expect = [(41, 100), (76, 135), (111, 170)]
        self.assertEqual(result, expect)

    def test_tile_length_larger(self):
        result = _calc_tiles(41, 49, 52, 0.9)
        expect = [(41, 49)]
        self.assertEqual(result, expect)

    def test_total_tile_length_1(self):
        result = _calc_tiles(41, 41, 52, 0.9)
        expect = [(41, 41)]
        self.assertEqual(result, expect)


class TestBuildBandedTable(ut.TestCase):
    def test_empty(self):
        result = _build_banded_table([], band_width=0)
        self.assertEqual(len(result.index), 0)
        self.assertListEqual(list(result.columns), [N_COL, *CONFUSION_COLS])

    def test_dedup_keeps_max_n(self):
        # The same pair observed in two overlapping tiles: keep the
        # observation with the greater coverage (N), and its whole row
        # (not merely N) -- only_a = n/2 distinguishes the two tiles' rows.
        table1 = _make_table([(1, 5, 10.0)])
        table2 = _make_table([(1, 5, 50.0)])
        result = _build_banded_table([table1, table2], band_width=0)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result[N_COL].iloc[0], 50.0)
        self.assertEqual(result[ONLY_A_COL].iloc[0], 25.0)

    def test_band_filter(self):
        rows = [(1, 5, 10.0), (1, 50, 10.0)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=10)
        self.assertListEqual(result.index.to_list(), [(1, 5)])

    def test_band_width_zero_applies_no_extra_cap(self):
        rows = [(1, 5, 10.0), (1, 50, 10.0)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=0)
        self.assertListEqual(result.index.to_list(), [(1, 5), (1, 50)])


class TestWritePairsWithConfusion(ut.TestCase):
    """pairs.csv must carry each written pair's 2x2 confusion-matrix counts,
    plus the p-value and fold change derived from them."""

    @staticmethod
    def _table(rows):
        # rows: (pos_a, pos_b, neither, only_a, only_b, both)
        import pandas as pd

        pa, pb, ne, oa, ob, bo = zip(*rows)
        index = pd.MultiIndex.from_arrays([pa, pb], names=[POSITION_A, POSITION_B])
        n = [a + b + c + d for a, b, c, d in zip(ne, oa, ob, bo)]
        return pd.DataFrame(
            {
                N_COL: n,
                NEITHER_COL: list(ne),
                ONLY_A_COL: list(oa),
                ONLY_B_COL: list(ob),
                BOTH_COL: list(bo),
            },
            index=index,
        )

    def test_writes_confusion_counts(self):
        import pandas as pd
        from scipy.stats import hypergeom

        pos_table = self._table([(1, 5, 900, 30, 20, 50), (2, 9, 800, 60, 40, 100)])
        # pvalue and adj_pvalue are passed in, as they would be from
        # _bridge_mask, rather than recomputed here.
        expected_pvals = [
            hypergeom.cdf(50, 1000, 80, 70),
            hypergeom.cdf(100, 1000, 160, 140),
        ]
        expected_adj_pvals = calc_bh_adjusted_pvals(np.array(expected_pvals))
        pvalue = pd.Series(expected_pvals, index=pos_table.index)
        adj_pvalue = pd.Series(expected_adj_pvals, index=pos_table.index)
        with tempfile.TemporaryDirectory() as tmp:
            csv_file = Path(tmp) / "pairs.csv"
            _write_pairs_with_confusion(pos_table, pvalue, adj_pvalue, csv_file)
            out = pd.read_csv(csv_file)
        # The position columns, the four confusion cells, and the p-value,
        # BH-adjusted p-value, and fold change derived from them.
        self.assertListEqual(
            list(out.columns),
            [
                POSITION_A,
                POSITION_B,
                *CONFUSION_COLS,
                P_VALUE_COL,
                ADJ_P_VALUE_COL,
                FOLD_CHANGE_COL,
            ],
        )
        self.assertListEqual(
            list(zip(out[POSITION_A], out[POSITION_B])), [(1, 5), (2, 9)]
        )
        self.assertListEqual(list(out[BOTH_COL]), [50, 100])
        self.assertListEqual(list(out[ONLY_A_COL]), [30, 60])
        # N is recoverable from the written cells.
        recovered_n = (
            out[NEITHER_COL] + out[ONLY_A_COL] + out[ONLY_B_COL] + out[BOTH_COL]
        )
        self.assertListEqual(list(recovered_n), [1000, 1000])
        # P-value: exact hypergeometric left-tail P(X <= both).
        np.testing.assert_allclose(out[P_VALUE_COL], expected_pvals)
        # BH-adjusted p-value: the value passed in.
        np.testing.assert_allclose(out[ADJ_P_VALUE_COL], expected_adj_pvals)
        # Fold change: independence expectation over the observed count.
        expected_fcs = [(80 * 70 / 1000) / 50, (160 * 140 / 1000) / 100]
        np.testing.assert_allclose(out[FOLD_CHANGE_COL], expected_fcs)


def _make_confusion_table(rows: list[tuple[int, int, float, float, float, float]]):
    """Build a per-pair table directly from (pos_a, pos_b, n, only_a, only_b,
    both), for tests that need precise control over the confusion cells (and
    so the exact fold change/hypergeometric significance) rather than the
    fixed ``n``-proportional split ``_make_table`` uses."""
    import pandas as pd

    pos_a: tuple[int, ...]
    pos_b: tuple[int, ...]
    ns: tuple[float, ...]
    only_a: tuple[float, ...]
    only_b: tuple[float, ...]
    both: tuple[float, ...]
    if not rows:
        pos_a, pos_b, ns, only_a, only_b, both = (), (), (), (), (), ()
    else:
        pos_a, pos_b, ns, only_a, only_b, both = zip(*rows)
    index = pd.MultiIndex.from_arrays([pos_a, pos_b], names=[POSITION_A, POSITION_B])
    neither = [n - (a + b - ab) for n, a, b, ab in zip(ns, only_a, only_b, both)]
    return pd.DataFrame(
        {
            N_COL: ns,
            NEITHER_COL: neither,
            ONLY_A_COL: only_a,
            ONLY_B_COL: only_b,
            BOTH_COL: both,
        },
        index=index,
    )


class TestBridgeMask(ut.TestCase):
    """``_bridge_mask`` marks a pair ELIGIBLE (enough coverage and
    expected-both, regardless of correlation direction) and, among those, a
    BRIDGE if it is also exact-hypergeometric-significant (``pair_fdr``) and
    depleted by at least ``min_fold_change`` relative to independence."""

    def test_strong_depletion_is_a_bridge(self):
        # n=100000, a=b=5000 (independence expects ab=250); observed ab=100
        # is a 2.5x depletion, and hugely significant at this scale.
        table = _make_confusion_table([(1, 2, 100000.0, 4900.0, 4900.0, 100.0)])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertListEqual(list(eligible), [True])
        self.assertListEqual(list(bridge), [True])

    def test_below_fold_change_floor_is_eligible_not_bridge(self):
        # n=10000, a=b=500 (independence expects ab=25); observed ab=20 is
        # anti-correlated and mildly depleted (fold change 1.25), below the
        # default 2.0 floor.
        table = _make_confusion_table([(1, 2, 10000.0, 480.0, 480.0, 20.0)])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertListEqual(list(eligible), [True])
        self.assertListEqual(list(bridge), [False])

    def test_positive_correlation_is_eligible_but_not_bridge(self):
        # observed ab=30 exceeds the independence expectation of 25: enough
        # coverage and expected-both to be eligible (and so to count in the
        # one-sided p-value family), but co-occurring rather than depleted, so
        # it fails both the left-tail significance and the effect-size floor
        # and is never a bridge.
        table = _make_confusion_table([(1, 2, 10000.0, 470.0, 470.0, 30.0)])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertListEqual(list(eligible), [True])
        self.assertListEqual(list(bridge), [False])

    def test_low_coverage_excluded_regardless_of_effect(self):
        table = _make_confusion_table([(1, 2, 500.0, 245.0, 245.0, 5.0)])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertListEqual(list(eligible), [False])
        self.assertListEqual(list(bridge), [False])

    def test_zero_both_mutated_is_an_infinite_fold_change_bridge(self):
        # Independence expects ab=9 (>= min_expect_both); observing zero
        # co-occurrences is an infinite fold change and, at this scale,
        # extremely significant.
        table = _make_confusion_table([(1, 2, 10000.0, 300.0, 300.0, 0.0)])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertListEqual(list(eligible), [True])
        self.assertListEqual(list(bridge), [True])

    def test_empty_table(self):
        table = _make_confusion_table([])
        eligible, bridge, _, _ = _bridge_mask(
            table,
            min_pair_coverage=1000,
            min_expect_both=5.0,
            pair_fdr=0.05,
            min_fold_change=2.0,
        )
        self.assertEqual(len(eligible), 0)
        self.assertEqual(len(bridge), 0)


def _eb(pairs):
    """Build ``(table, eligible, bridge)`` from ``pairs`` -- a list of
    ``(a, b, is_bridge)`` -- treating every listed pair as eligible."""
    index = pd.MultiIndex.from_tuples(
        [(a, b) for a, b, _ in pairs], names=[POSITION_A, POSITION_B]
    )
    table = pd.DataFrame(index=index)
    eligible = pd.Series(True, index=index)
    bridge = pd.Series([bool(x) for *_, x in pairs], index=index)
    return table, eligible, bridge


def _dense_block(lo, hi, band, rng, p_in, p_bg, total, cross=None, p_cross=0.0):
    """Yield ``(a, b, is_bridge)`` for every in-band pair of [1, total], with
    bridge probability ``p_in`` inside [lo, hi], ``p_bg`` in the background, and
    ``p_cross`` for pairs straddling the ``cross`` region ``(clo, chi)``."""
    for a in range(1, total + 1):
        for b in range(a + 1, min(a + band, total) + 1):
            if lo <= a <= b <= hi:
                p = p_in
            else:
                p = p_bg
            if cross is not None and a <= cross[0] and b >= cross[1]:
                p = p_cross
            yield (a, b, rng.random() < p)


class TestBlockScore(ut.TestCase):
    """``_block_score`` is the BIC-corrected binomial log-likelihood ratio of a
    block's own bridge rate against the null rate ``pi0``."""

    def test_matches_hand_computed_llr_minus_bic(self):
        pi0 = 0.01
        ne, nb = 1000, 50
        pi_hat = nb / ne
        llr = nb * np.log(pi_hat / pi0) + (ne - nb) * np.log((1 - pi_hat) / (1 - pi0))
        expected = llr - np.log(ne) / 2
        got = _block_score(np.array([ne]), np.array([nb]), pi0)[0]
        self.assertAlmostEqual(got, expected, places=6)

    def test_block_at_or_below_null_scores_negative(self):
        pi0 = 0.02
        # pi_hat = pi0 exactly -> llr 0 -> gain = -bic < 0.
        at_null = _block_score(np.array([500]), np.array([10]), pi0)[0]
        self.assertLess(at_null, 0.0)
        # pi_hat below pi0 -> clamped to pi0 -> still negative.
        below = _block_score(np.array([500]), np.array([2]), pi0)[0]
        self.assertLess(below, 0.0)

    def test_enrichment_increases_gain(self):
        pi0 = 0.01
        low = _block_score(np.array([1000]), np.array([20]), pi0)[0]
        high = _block_score(np.array([1000]), np.array([80]), pi0)[0]
        self.assertGreater(high, low)

    def test_empty_block_is_minus_inf(self):
        self.assertEqual(_block_score(np.array([0]), np.array([0]), 0.01)[0], -np.inf)

    def test_locality_independent_of_other_blocks(self):
        # The score of a block depends only on its own (n_elig, n_bridge) and
        # pi0 -- a fixed constant -- never on any other block. Directly
        # regression-tests failure #1 (global-rate contamination).
        pi0 = 0.01
        s1 = _block_score(np.array([800]), np.array([12]), pi0)[0]
        s2 = _block_score(np.array([800]), np.array([12]), pi0)[0]
        self.assertEqual(s1, s2)


class TestTriangleSumBanded(ut.TestCase):
    """``_triangle_sum_banded`` reads triangle (contained-pair) sums off the
    banded row-cumulative arrays; the ``s_min`` bound matches the unbounded
    result on the overlapping range."""

    def _grids(self, total=40, band=15, seed=0):
        rng = np.random.default_rng(seed)
        pairs = [
            (a, b, rng.random() < 0.25)
            for a in range(1, total + 1)
            for b in range(a + 1, min(a + band, total) + 1)
        ]
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        return _pair_band_row_cumsum(table, 1, total, value_col=BRIDGE_COL), pairs

    def test_matches_brute_force_contained_counts(self):
        (count_rc, value_rc, npos, mg), pairs = self._grids()
        pa = np.array([a for a, _, _ in pairs])
        pb = np.array([b for _, b, _ in pairs])
        bv = np.array([float(x) for *_, x in pairs])
        for s, e in [(1, 10), (5, 25), (12, 40)]:
            m = (pa >= s) & (pb <= e)
            ne_bf, nb_bf = int(m.sum()), float(bv[m].sum())
            ne_ts = _triangle_sum_banded(count_rc, e - 1, mg, s_min=s - 1)[s - 1]
            nb_ts = _triangle_sum_banded(value_rc, e - 1, mg, s_min=s - 1)[s - 1]
            self.assertEqual(ne_bf, int(ne_ts))
            self.assertTrue(np.isclose(nb_bf, nb_ts))

    def test_s_min_equals_unbounded_on_overlap(self):
        (count_rc, _, npos, mg), _ = self._grids()
        for e in [4, 18, 39]:
            full = _triangle_sum_banded(count_rc, e, mg, s_min=0)
            for s_min in [0, max(0, e - 8), max(0, e - 2)]:
                bounded = _triangle_sum_banded(count_rc, e, mg, s_min=s_min)
                self.assertTrue(
                    np.array_equal(full[s_min : e + 1], bounded[s_min : e + 1])
                )


class TestCutCrossingScores(ut.TestCase):
    """``_cut_crossing_scores`` counts, per cut, the eligible pairs straddling
    it and how many are bridges (difference-array sweep)."""

    def test_matches_brute_force(self):
        pairs = [(1, 4, True), (2, 6, False), (3, 5, True), (5, 8, True)]
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        ne_cross, nb_cross = _cut_crossing_scores(table, 1, 8)
        # cut m (index m-1) is between positions m and m+1.
        for m in range(1, 8):
            elig_bf = sum(1 for a, b, _ in pairs if a <= m < b)
            brid_bf = sum(1 for a, b, x in pairs if a <= m < b and x)
            self.assertEqual(int(ne_cross[m - 1]), elig_bf)
            self.assertEqual(int(round(nb_cross[m - 1])), brid_bf)

    def test_zero_crossing_reads_exactly_zero(self):
        pairs = [(1, 3, True), (5, 8, True)]  # nothing crosses cut 4 (3|4)
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        ne_cross, nb_cross = _cut_crossing_scores(table, 1, 8)
        self.assertEqual(int(ne_cross[3]), 0)  # cut between 4 and 5
        self.assertEqual(nb_cross[3], 0.0)


class TestCalcBlockPvalueCutoff(ut.TestCase):
    """``_calc_block_pvalue_cutoff`` returns the BH cutoff over every candidate
    block; -1 (admit nothing) when no block is significant."""

    def _grids(self, pairs, total):
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        return _pair_band_row_cumsum(table, 1, total, value_col=BRIDGE_COL)

    def test_admit_nothing_on_pure_background(self):
        rng = np.random.default_rng(0)
        total, band = 60, 20
        pairs = list(_dense_block(0, 0, band, rng, 0.0, 0.0, total))  # no bridges
        count_rc, value_rc, npos, mg = self._grids(pairs, total)
        cutoff = _calc_block_pvalue_cutoff(
            count_rc, value_rc, npos, mg, 0.01, 0.05, 200
        )
        self.assertEqual(cutoff, -1.0)

    def test_dense_block_admitted(self):
        rng = np.random.default_rng(1)
        total, band = 60, 20
        pairs = list(_dense_block(20, 40, band, rng, 0.4, 0.0, total))
        count_rc, value_rc, npos, mg = self._grids(pairs, total)
        cutoff = _calc_block_pvalue_cutoff(
            count_rc, value_rc, npos, mg, 0.01, 0.05, 200
        )
        self.assertGreater(cutoff, 0.0)


class TestDpSegmentBlocks(ut.TestCase):
    """``_dp_segment_blocks`` finds the globally-optimal admitted partition and
    honours ``max_domain_length``."""

    def _grids(self, pairs, total):
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        return _pair_band_row_cumsum(table, 1, total, value_col=BRIDGE_COL)

    def test_finds_dense_block(self):
        rng = np.random.default_rng(2)
        total, band = 80, 20
        pairs = list(_dense_block(30, 55, band, rng, 0.4, 0.002, total))
        count_rc, value_rc, npos, mg = self._grids(pairs, total)
        cutoff = _calc_block_pvalue_cutoff(
            count_rc, value_rc, npos, mg, 0.005, 0.05, 200
        )
        blocks, _ = _dp_segment_blocks(count_rc, value_rc, npos, mg, 0.005, cutoff, 200)
        # exactly one block, overlapping the dense region [30,55] (0-indexed
        # 29..54).
        self.assertEqual(len(blocks), 1)
        s, e = blocks[0]
        self.assertLessEqual(s, 29)
        self.assertGreaterEqual(e, 54)

    def test_respects_max_domain_length(self):
        rng = np.random.default_rng(3)
        total, band = 80, 30
        pairs = list(_dense_block(10, 70, band, rng, 0.4, 0.002, total))
        count_rc, value_rc, npos, mg = self._grids(pairs, total)
        cutoff = _calc_block_pvalue_cutoff(
            count_rc, value_rc, npos, mg, 0.005, 0.05, 25
        )
        blocks, _ = _dp_segment_blocks(count_rc, value_rc, npos, mg, 0.005, cutoff, 25)
        self.assertTrue(blocks)
        for s, e in blocks:
            self.assertLessEqual(e - s + 1, 25)


class TestMergeConnectedBlocks(ut.TestCase):
    """``_merge_connected_blocks`` joins adjacent blocks whose gap is crossed by
    enriched bridges at every cut, subject to ``max_domain_length``."""

    def _cross(self, total, per_cut_bridge, per_cut_elig, gap):
        # Build crossing arrays where every cut in ``gap`` (inclusive range of
        # cut indices) has (per_cut_elig, per_cut_bridge), others are null.
        n_cuts = total - 1
        ne = np.zeros(n_cuts, dtype=np.int64)
        nb = np.zeros(n_cuts, dtype=float)
        for c in range(gap[0], gap[1] + 1):
            ne[c] = per_cut_elig
            nb[c] = per_cut_bridge
        return ne, nb

    def test_merges_when_every_gap_cut_connected(self):
        blocks = [(0, 20), (30, 50)]  # gap cuts 20..29
        ne, nb = self._cross(60, per_cut_bridge=40, per_cut_elig=100, gap=(20, 29))
        merged = _merge_connected_blocks(blocks, ne, nb, 0.01, 0.05, 200)
        self.assertEqual(merged, [(0, 50)])

    def test_keeps_split_when_a_cut_is_disconnected(self):
        blocks = [(0, 20), (30, 50)]
        ne, nb = self._cross(60, per_cut_bridge=40, per_cut_elig=100, gap=(20, 29))
        # blank out one interior cut -> not connected everywhere.
        nb[25] = 0.0
        merged = _merge_connected_blocks(blocks, ne, nb, 0.01, 0.05, 200)
        self.assertEqual(merged, [(0, 20), (30, 50)])

    def test_max_domain_length_guard_blocks_merge(self):
        blocks = [(0, 20), (30, 50)]
        ne, nb = self._cross(60, per_cut_bridge=40, per_cut_elig=100, gap=(20, 29))
        merged = _merge_connected_blocks(blocks, ne, nb, 0.01, 0.05, 40)
        self.assertEqual(merged, [(0, 20), (30, 50)])

    def test_faint_uniform_gap_does_not_merge(self):
        # Crossing rate barely above pi0 but not significantly enriched -> the
        # effect-size/BH gate holds the boundary (no over-merge).
        blocks = [(0, 20), (30, 50)]
        ne, nb = self._cross(60, per_cut_bridge=1, per_cut_elig=100, gap=(20, 29))
        merged = _merge_connected_blocks(blocks, ne, nb, 0.01, 0.05, 200)
        self.assertEqual(merged, [(0, 20), (30, 50)])


class TestExtendDomainsByBridges(ut.TestCase):
    """``_extend_domains_by_bridges`` grows a domain outward to absorb a thin
    line of edge bridges, but ignores a lone distant stray bridge."""

    def _table(self, pairs):
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        return table

    def _core(self):
        # Dense core [50, 80], all pairs bridges (band 20).
        return [
            (a, b, True) for a in range(50, 81) for b in range(a + 1, min(a + 20, 81))
        ]

    def _sparse_bg(self, lo, hi, exclude):
        # A sparse eligible (non-bridge) background between lo and hi: every 7th
        # in-band pair, mimicking real data where few pairs clear eligibility.
        out = []
        k = 0
        for a in range(lo, hi + 1):
            for b in range(a + 1, min(a + 20, 90) + 1):
                k += 1
                if k % 7 == 0 and (a, b) not in exclude:
                    out.append((a, b, False))
        return out

    def test_absorbs_thin_line_of_edge_bridges(self):
        # A thin line of 8 bridges reaches into the core from positions 10-13,
        # across sparse empty space 14-49.
        thin = [
            (10, 55),
            (11, 56),
            (12, 57),
            (13, 58),
            (10, 60),
            (11, 61),
            (12, 62),
            (13, 63),
        ]
        pairs = (
            self._core()
            + self._sparse_bg(1, 49, {(a, b) for a, b in thin})
            + [(a, b, True) for a, b in thin]
        )
        got = _extend_domains_by_bridges(
            [(50, 80)],
            self._table(pairs),
            pi0=0.01,
            extend_fdr=0.05,
            max_domain_length=200,
            total_end5=1,
            total_end3=90,
        )
        ((s, e),) = got
        self.assertLessEqual(s, 13, f"did not absorb the thin line: {got}")
        self.assertGreaterEqual(s, 10)

    def test_ignores_lone_distant_bridge(self):
        pairs = (
            self._core()
            + self._sparse_bg(1, 49, {(12, 60)})
            + [(12, 60, True)]  # a single stray bridge into the core
        )
        got = _extend_domains_by_bridges(
            [(50, 80)],
            self._table(pairs),
            pi0=0.01,
            extend_fdr=0.05,
            max_domain_length=200,
            total_end5=1,
            total_end3=90,
        )
        ((s, e),) = got
        self.assertGreater(s, 13, f"a lone stray bridge should not extend: {got}")

    def test_extension_respects_max_domain_length(self):
        thin = [
            (10, 55),
            (11, 56),
            (12, 57),
            (13, 58),
            (10, 60),
            (11, 61),
            (12, 62),
            (13, 63),
        ]
        pairs = (
            [(a, b, True) for a in range(50, 71) for b in range(a + 1, min(a + 20, 71))]
            + self._sparse_bg(1, 49, {(a, b) for a, b in thin})
            + [(a, b, True) for a, b in thin]
        )
        got = _extend_domains_by_bridges(
            [(50, 70)],
            self._table(pairs),
            pi0=0.01,
            extend_fdr=0.05,
            max_domain_length=25,
            total_end5=1,
            total_end3=90,
        )
        ((s, e),) = got
        self.assertLessEqual(e - s + 1, 25)  # capped, never longer


class TestEstimateNullBridgeRate(ut.TestCase):
    """``_estimate_null_bridge_rate`` estimates pi0 as the out-of-domain bridge
    rate, floored."""

    def test_out_of_domain_rate(self):
        # 100 background eligible pairs, 4 bridges; all in-domain bridges
        # excluded.
        pairs = []
        # domain [10,20]: many bridges (should be excluded)
        for k in range(20):
            pairs.append((10 + (k % 5), 16 + (k % 5), True))
        # background: 100 pairs, 4 bridges
        for k in range(100):
            pairs.append((30 + k % 20, 55 + k % 20, k < 4))
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        rate = _estimate_null_bridge_rate(table, [(10, 20)], floor=1e-6)
        self.assertAlmostEqual(rate, 4 / 100, places=6)

    def test_floor_applied_when_background_clean(self):
        pairs = [(30, 40, False), (31, 41, False)]
        table, _, bridge = _eb(pairs)
        table[BRIDGE_COL] = bridge.to_numpy(dtype=float)
        rate = _estimate_null_bridge_rate(table, [], floor=0.003)
        self.assertEqual(rate, 0.003)


class TestGapUtilities(ut.TestCase):
    """The gap-handling helpers (``_split``/``_fill``/``_widen``/``_label``/
    ``_filter``)."""

    def test_filter_drops_short_domains(self):
        got = _filter_domains_length([(1, 5), (10, 40), (45, 46)], min_length=10)
        self.assertEqual(got, [(10, 40)])

    def test_split_gap_evenly_divides_exactly(self):
        got = _split_gap_evenly(1, 20, 10)
        self.assertEqual(got, [(1, 10), (11, 20)])

    def test_split_gap_evenly_minimum_count_as_equal_as_possible(self):
        got = _split_gap_evenly(1, 25, 10)
        # ceil(25 / 10) == 3 pieces; sizes differ by at most 1 (9, 8, 8).
        self.assertEqual(got, [(1, 9), (10, 17), (18, 25)])
        # No piece exceeds the cap, and the pieces exactly reconstitute the
        # input range with no gaps or overlaps.
        self.assertTrue(all(e3 - e5 + 1 <= 10 for e5, e3 in got))
        self.assertEqual(got[0][0], 1)
        self.assertEqual(got[-1][1], 25)
        for (_, e3), (e5, _) in zip(got, got[1:]):
            self.assertEqual(e3 + 1, e5)

    def test_fill_matches_old_insert_behavior_when_cap_is_non_binding(self):
        got = _fill_domains_into_gaps([(10, 20), (31, 40)], 1, 50, 100)
        self.assertEqual(got, [(1, 9), (10, 20), (21, 30), (31, 40), (41, 50)])

    def test_fill_splits_oversized_leading_gap(self):
        got = _fill_domains_into_gaps([(21, 30)], 1, 30, 10)
        self.assertEqual(got, [(1, 10), (11, 20), (21, 30)])

    def test_fill_splits_oversized_interior_gap(self):
        got = _fill_domains_into_gaps([(1, 10), (41, 50)], 1, 50, 15)
        self.assertEqual(got, [(1, 10), (11, 25), (26, 40), (41, 50)])

    def test_fill_leaves_adjacent_domains_untouched(self):
        got = _fill_domains_into_gaps([(1, 10), (11, 20)], 1, 20, 5)
        self.assertEqual(got, [(1, 10), (11, 20)])

    def test_fill_covers_whole_region_when_no_domains(self):
        got = _fill_domains_into_gaps([], 1, 25, 10)
        self.assertEqual(got, [(1, 9), (10, 17), (18, 25)])

    def test_widen_matches_old_expand_behavior_when_cap_is_non_binding(self):
        got = _widen_domains_into_gaps([(10, 20), (31, 40)], 1, 50, 100)
        self.assertEqual(got, [(1, 25), (26, 50)])

    def test_widen_binding_cap_leaves_leading_leftover(self):
        got = _widen_domains_into_gaps([(10, 20)], 1, 100, 15)
        # Budget is only 15 - 11 = 4, so the domain can absorb only 4 of the
        # 9 positions in the leading gap, and has nothing left for the
        # (much larger) trailing gap.
        self.assertEqual(got, [(6, 20)])

    def test_widen_binding_cap_leaves_interior_leftover(self):
        got = _widen_domains_into_gaps([(1, 10), (41, 50)], 1, 50, 15)
        # Each domain's budget (5) covers only half of its half of the gap,
        # leaving a residual gap of (16, 35) between the widened domains.
        self.assertEqual(got, [(1, 15), (36, 50)])

    def test_widen_then_fill_mops_up_residual_gap(self):
        widened = _widen_domains_into_gaps([(1, 10), (41, 50)], 1, 50, 15)
        got = _fill_domains_into_gaps(widened, 1, 50, 15)
        self.assertEqual(got, [(1, 15), (16, 25), (26, 35), (36, 50)])

    def test_label_widened_marks_unchanged_and_grown(self):
        raw = [(1, 10), (20, 30)]
        widened = [(1, 10), (15, 30)]
        got = _label_widened(raw, widened)
        self.assertEqual(got, {(1, 10): "original", (15, 30): "widened"})


class TestCalcDomainsByDpSegmentation(ut.TestCase):
    """End-to-end DP domain calling on constructed bridge tables."""

    def test_two_domains_recovered_and_not_merged(self):
        rng = np.random.default_rng(10)
        total, band = 90, 20
        pairs = []
        for a in range(1, total + 1):
            for b in range(a + 1, min(a + band, total) + 1):
                if 10 <= a <= b <= 30 or 51 <= a <= b <= 70:
                    p = 0.3
                else:
                    p = 0.003
                pairs.append((a, b, rng.random() < p))
        table, eligible, bridge = _eb(pairs)
        domains, pi0 = _calc_domains_by_dp_segmentation(
            table,
            eligible,
            bridge,
            1,
            total,
            pair_fdr=0.05,
            detect_fdr=0.05,
            merge_fdr=0.05,
            max_domain_length=200,
        )
        big = [(s, e) for s, e in domains if e - s + 1 >= 10]
        self.assertTrue(any(s <= 10 and e >= 30 for s, e in big), big)
        self.assertTrue(any(s <= 51 and e >= 70 for s, e in big), big)
        self.assertFalse(any(s <= 10 and e >= 70 for s, e in big), big)

    def test_over_fragmentation_merge_across_spanned_void(self):
        # Two dense blocks separated by a void that long-range bridges span at
        # every cut -> the cut-crossing merge rejoins them into one domain.
        rng = np.random.default_rng(11)
        total, band = 90, 40
        pairs = []
        A, B = (10, 30), (46, 66)
        for a in range(1, total + 1):
            for b in range(a + 1, min(a + band, total) + 1):
                if A[0] <= a <= b <= A[1] or B[0] <= a <= b <= B[1]:
                    p = 0.3
                elif a <= A[1] and b >= B[0]:  # long-range crossing bridges
                    p = 0.5
                else:
                    p = 0.003
                pairs.append((a, b, rng.random() < p))
        table, eligible, bridge = _eb(pairs)
        domains, _ = _calc_domains_by_dp_segmentation(
            table,
            eligible,
            bridge,
            1,
            total,
            pair_fdr=0.05,
            detect_fdr=0.05,
            merge_fdr=0.05,
            max_domain_length=200,
        )
        big = [(s, e) for s, e in domains if e - s + 1 >= 10]
        self.assertEqual(len(big), 1, f"expected one merged domain, got {big}")
        s, e = big[0]
        self.assertLessEqual(s, A[0])
        self.assertGreaterEqual(e, B[1])

    def test_separate_domains_stay_split_without_crossing(self):
        rng = np.random.default_rng(11)
        total, band = 90, 40
        pairs = []
        A, B = (10, 30), (46, 66)
        for a in range(1, total + 1):
            for b in range(a + 1, min(a + band, total) + 1):
                if A[0] <= a <= b <= A[1] or B[0] <= a <= b <= B[1]:
                    p = 0.3
                else:
                    p = 0.003  # no crossing bridges
                pairs.append((a, b, rng.random() < p))
        table, eligible, bridge = _eb(pairs)
        domains, _ = _calc_domains_by_dp_segmentation(
            table,
            eligible,
            bridge,
            1,
            total,
            pair_fdr=0.05,
            detect_fdr=0.05,
            merge_fdr=0.05,
            max_domain_length=200,
        )
        big = [(s, e) for s, e in domains if e - s + 1 >= 10]
        self.assertEqual(len(big), 2, f"expected two domains, got {big}")

    def test_block_locality_distant_dense_block_does_not_change_call(self):
        # A near real domain must still be called (with the same dense core)
        # whether or not a distant dense block is present -- regression for
        # failure #1 (a distant domain cannot corrupt an unrelated region's
        # call). The exact extended boundary may shift slightly because pi0 is a
        # global background estimate, so assert the core is robustly recovered
        # rather than byte-identical coordinates.
        def run(with_distant):
            rng = np.random.default_rng(12)
            total, band = 120, 20
            pairs = []
            for a in range(1, total + 1):
                for b in range(a + 1, min(a + band, total) + 1):
                    if 10 <= a <= b <= 30:  # the near real domain
                        p = 0.3
                    elif with_distant and 90 <= a <= b <= 115:  # distant dense
                        p = 0.5
                    else:
                        p = 0.003
                    pairs.append((a, b, rng.random() < p))
            table, eligible, bridge = _eb(pairs)
            domains, _ = _calc_domains_by_dp_segmentation(
                table,
                eligible,
                bridge,
                1,
                total,
                pair_fdr=0.05,
                detect_fdr=0.05,
                merge_fdr=0.05,
                max_domain_length=200,
            )
            # the near domain (overlapping [10,30]), restricted below position 80
            return [(s, e) for s, e in domains if e < 80 and e - s + 1 >= 10]

        for res in (run(False), run(True)):
            self.assertEqual(len(res), 1, res)
            s, e = res[0]
            self.assertLessEqual(s, 12, res)  # core start recovered
            self.assertGreaterEqual(e, 28, res)  # core end recovered

    def test_max_domain_length_bounds_final_domains(self):
        rng = np.random.default_rng(13)
        total, band = 100, 30
        pairs = []
        for a in range(1, total + 1):
            for b in range(a + 1, min(a + band, total) + 1):
                p = 0.35 if 10 <= a <= b <= 90 else 0.003
                pairs.append((a, b, rng.random() < p))
        table, eligible, bridge = _eb(pairs)
        domains, _ = _calc_domains_by_dp_segmentation(
            table,
            eligible,
            bridge,
            1,
            total,
            pair_fdr=0.05,
            detect_fdr=0.05,
            merge_fdr=0.05,
            max_domain_length=30,
        )
        for s, e in domains:
            self.assertLessEqual(e - s + 1, 30)

    def test_out_of_order_ends_raise(self):
        table, eligible, bridge = _eb([(1, 3, True)])
        with self.assertRaises(IncompatibleValuesError):
            _calc_domains_by_dp_segmentation(
                table,
                eligible,
                bridge,
                10,
                1,
                pair_fdr=0.05,
                detect_fdr=0.05,
                merge_fdr=0.05,
                max_domain_length=100,
            )

    def test_bad_max_domain_length_raises(self):
        table, eligible, bridge = _eb([(1, 3, True)])
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                eligible,
                bridge,
                1,
                10,
                pair_fdr=0.05,
                detect_fdr=0.05,
                merge_fdr=0.05,
                max_domain_length=0,
            )


class ScanTestBase(ut.TestCase):
    """Shared simulation infrastructure for the scan tests."""

    SAMPLE = "test_sample"
    REFS = "test_refs"
    REF = "test_ref"
    PROFILE = "scan"

    # Folding domains of the reference sequence (each 60 nt).
    DOMAINS = [
        (
            "TGACGAACAACGTGTTTGTGAACCATATAGGTAAACGCTGAATGCGTTCGCGCGGAGGGT",
            [
                "..(((((((...)))))))..(((......((.(((((.....))))).))......)))",
                "...(((...(((((((((((.(((.....)))...))).))))))))))).((.....))",
            ],
        ),
        (
            "TTTGCAGGAAGATGGTCAACTCTACACCTAGTTTTTACCAGTCCACAAGAGTTTGAACTG",
            [".(((..(((...((((.((..(((....)))..)).)))).))).))).(((....)))."],
        ),
        (
            "GTGCCTTAACCTGAGTACGCCCATATCATGGGAGACATTACAACTCAAATTCTAGGTGTG",
            [
                "..((((.....(((((...(((((...)))))..........)))))......))))...",
                "((((.(((...)))))))...((((((.((((..................))))))))))",
            ],
        ),
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None
        self._tmpdir = None

    @property
    def sim_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name)

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)
        self._tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmpdir.cleanup()
        self._tmpdir = None
        set_config(*self._config)

    def sim_data(self, domain_nums: list[int], read_length: int, seed: int):
        # Assemble and write the reference sequence.
        domains = dict(self.DOMAINS[m] for m in domain_nums)
        refseq = DNA("".join(domains.keys()))
        refs_dir = self.sim_dir.joinpath("refs")
        refs_dir.mkdir()
        fasta = refs_dir.joinpath(f"{self.REFS}.fa")
        write_fasta(fasta, [(self.REF, refseq)])
        # Assemble and write the secondary structures.
        structures = list(map("".join, product(*domains.values())))
        param_dir = self.sim_dir.joinpath("params", self.REF, FULL_NAME)
        param_dir.mkdir(parents=True)
        db_file = param_dir.joinpath(f"{self.PROFILE}.db")
        with open(db_file, "x") as f:
            for i, struct in enumerate(structures):
                if i == 0:
                    f.write(f">structure0\n{refseq.tr()}\n{struct}\n")
                else:
                    f.write(f">structure{i}\n{struct}\n")
        ct_file = db_to_ct(db_file)
        # Simulate data.
        sim_params(
            ct_file=[ct_file],
            # Make pmut_unpaired for A and C large (17%) so that
            # most reads get at least two mutations despite being
            # short and are thus useful for clustering.
            pmut_unpaired=[("am", 1 / 6), ("cm", 1 / 6)],
            # Make all reads the same length.
            length_fmean=(read_length / len(refseq)),
            length_fvar=0.0,
            # Make clust_conc very large so that the proportion
            # of each cluster is approximately equal, which makes
            # clustering easier.
            clust_conc=1000.0,
            seed=seed,
        )
        idmut_dirs = sim_idmut(
            param_dir=[param_dir],
            sample=self.SAMPLE,
            profile_name=self.PROFILE,
            num_reads=200000,
            paired_end=False,
            brotli_level=0,
            seed=seed,
        )
        return idmut_dirs


class TestFilterScan(ScanTestBase):
    """Test that filterscan identifies domains without clustering."""

    def run_filterscan_check(
        self, idmut_dirs: list[Path], expect_regions: list[tuple[int, int]], **kwargs
    ):
        filterscan_dirs = run_filterscan(
            idmut_dirs,
            brotli_level=0,
            filter_pos_table=False,
            filter_read_table=False,
            **kwargs,
        )
        domains = list()
        for filterscan_dir in filterscan_dirs:
            for report_file in path.find_files_chain(
                [filterscan_dir], FilterScanReport.get_path_seg_types()
            ):
                report = FilterScanReport.load(report_file)
                for reg5, reg3 in report.get_field(DomainCoordsF):
                    domains.append((int(reg5), int(reg3)))
        for exp5, exp3 in expect_regions:
            expect_length = exp3 - exp5 + 1
            for reg5, reg3 in domains:
                overlap = max(min(reg3, exp3) - max(reg5, exp5), 0)
                if overlap >= expect_length / 2:
                    break
            else:
                raise ValueError(
                    f"Expected region {exp5, exp3} does not overlap at "
                    f"least 50% of any detected domain among {sorted(domains)}"
                )

    def test_domains012_read180(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        self.run_filterscan_check(
            idmut_dirs, [(1, 60), (121, 180)], min_domain_length=20
        )

    def test_domains012_read120(self):
        idmut_dirs = self.sim_data([0, 1, 2], 120, seed=0)
        self.run_filterscan_check(
            idmut_dirs, [(1, 60), (121, 180)], min_domain_length=20
        )

    def test_domains012_read60(self):
        idmut_dirs = self.sim_data([0, 1, 2], 60, seed=0)
        self.run_filterscan_check(
            idmut_dirs, [(1, 60), (121, 180)], min_domain_length=20
        )

    def test_domains02_read60(self):
        idmut_dirs = self.sim_data([0, 2], 60, seed=0)
        self.run_filterscan_check(
            idmut_dirs, [(1, 60), (61, 120)], min_domain_length=20
        )

    def test_domains012_read180_cli(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        runner = CliRunner()
        args = (
            ["-qq", "--exit-on-error", "filterscan"]
            + [str(d) for d in idmut_dirs]
            + ["--brotli-level", "0"]
        )
        result = runner.invoke(seismic_cli, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, msg=result.output)
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)


class TestSplitAuthenticity(ScanTestBase):
    """Validate that a genuine domain split is authentic: if two domains are
    truly independent structural units, clustering them separately should
    recover -- up to EM/sampling noise -- the same joint mixture as
    clustering the merged region. Concretely, the merged region's cluster
    count and proportions should match the CARTESIAN PRODUCT of the two
    domains' independent clusterings: ``K_merged == K_A * K_B`` and the
    merged proportions should equal the outer product of the two domains'
    marginal proportions. This is a general check of the splitting concept
    (not specific to the insulation-based caller): an authentic boundary is
    one where the two sides are statistically independent, which is exactly
    what the cartesian-product prediction tests."""

    def _cluster_region(self, idmut_dirs, end5: int, end3: int, max_clusters: int):
        filter_dirs = run_filter(
            idmut_dirs,
            region_coords=[(self.REF, end5, end3)],
            filter_pos_table=False,
            filter_read_table=False,
            brotli_level=0,
        )
        cluster_dirs = run_cluster(
            filter_dirs,
            max_clusters=max_clusters,
            min_em_runs=1,
            max_em_runs=1,
            jackpot=False,
            cluster_pos_table=False,
            cluster_abundance_table=False,
            brotli_level=0,
            seed=0,
        )
        report_file = cluster_dirs[0].joinpath("cluster-report.json")
        dataset = ClusterMutsDataset(report_file)
        params = get_clust_params(dataset)
        proportions = params.loc[(0, "p")].to_numpy(dtype=float)
        return dataset.best_k, proportions

    def test_merged_region_matches_cartesian_product_of_separate_domains(self):
        # DOMAINS[0] and DOMAINS[2] each independently fold into 2
        # alternative structures, and are structurally unrelated to each
        # other (no shared base pairing), so they are a genuine authentic-
        # split pair. Use full-length reads (read_length == the merged
        # reference's own length) so every read spans BOTH domains and can
        # carry linking evidence about their JOINT cluster assignment --
        # without that, no single read could ever reveal which domain-A
        # structure co-occurs with which domain-B structure, and the merged
        # clustering could not resolve a 4-way joint mixture at all.
        idmut_dirs = self.sim_data([0, 2], 120, seed=0)
        k_a, pi_a = self._cluster_region(idmut_dirs, 1, 60, max_clusters=4)
        k_b, pi_b = self._cluster_region(idmut_dirs, 61, 120, max_clusters=4)
        k_merged, pi_merged = self._cluster_region(idmut_dirs, 1, 120, max_clusters=8)
        self.assertEqual(
            k_merged,
            k_a * k_b,
            f"Expected the merged region's cluster count ({k_merged}) to "
            f"equal the product of the separate domains' ({k_a} * {k_b})",
        )
        # Cluster labels are not aligned across the separate and merged
        # runs, so compare the proportions as sorted, unlabeled sets: an
        # authentic split predicts the merged proportions equal the outer
        # product of the two domains' marginal proportions.
        outer = np.outer(pi_a, pi_b).flatten()
        np.testing.assert_allclose(np.sort(pi_merged), np.sort(outer), atol=0.05)


if __name__ == "__main__":
    ut.main(verbosity=2)
