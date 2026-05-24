import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.confusion import (
    init_confusion_matrix,
    calc_confusion_matrix,
    calc_confusion_phi,
    calc_bh_adjusted_pvals,
    _count_intersection,
)
from seismicrna.core.batch.count import calc_covered_reads_per_pos, calc_reads_per_pos
from seismicrna.core.header import ClustHeader
from seismicrna.core.rel.code import DELET, SUB_A, SUB_C, SUB_G, SUB_T
from seismicrna.core.rel.pattern import RelPattern
from seismicrna.core.seq.region import POS_NAME, BASE_NAME


class TestCountIntersection(ut.TestCase):
    def test_empty(self):
        x = np.array([])
        y = np.array([])
        self.assertEqual(_count_intersection(x, y), 0)

    def test_0_overlap(self):
        x = np.array([0, 2, 4, 6])
        y = np.array([1, 3, 5])
        self.assertEqual(_count_intersection(x, y), 0)
        self.assertEqual(_count_intersection(y, x), 0)

    def test_first_overlap(self):
        x = np.array([1, 2, 4, 6])
        y = np.array([1, 3, 5])
        self.assertEqual(_count_intersection(x, y), 1)
        self.assertEqual(_count_intersection(y, x), 1)

    def test_last_overlap(self):
        x = np.array([0, 2, 4, 6])
        y = np.array([1, 3, 6])
        self.assertEqual(_count_intersection(x, y), 1)
        self.assertEqual(_count_intersection(y, x), 1)

    def test_middle_overlap(self):
        x = np.array([0, 2, 4, 6, 7, 8, 9, 10, 12, 14])
        y = np.array([1, 3, 6, 7, 8, 9, 11, 13])
        self.assertEqual(_count_intersection(x, y), 4)
        self.assertEqual(_count_intersection(y, x), 4)

    def test_count_intersect(self):
        rng = np.random.default_rng(seed=0)
        x = np.unique(rng.integers(100, size=100))
        y = np.unique(rng.integers(100, size=100))
        count = np.intersect1d(x, y).size
        self.assertEqual(_count_intersection(x, y), count)
        self.assertEqual(_count_intersection(y, x), count)


class TestInitConfusionMatrix(ut.TestCase):
    def test_no_clusters(self):
        pos_index = pd.MultiIndex.from_arrays(
            [[2, 4, 7, 9], ["C", "C", "A", "A"]], names=[POS_NAME, BASE_NAME]
        )
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 2, 4, 4, 7], [4, 7, 9, 7, 9, 9]], names=["Position A", "Position B"]
        )
        n, a, b, ab = init_confusion_matrix(pos_index)
        for entry in [n, a, b, ab]:
            self.assertIsInstance(entry, pd.Series)
            self.assertTrue(entry.index.equals(expect_index))
            self.assertTrue(entry.index.equal_levels(expect_index))
            self.assertTrue(np.all(entry == 0))

    def test_no_clusters_min_gap(self):
        pos_index = pd.MultiIndex.from_arrays(
            [[2, 4, 7, 9], ["C", "C", "A", "A"]], names=[POS_NAME, BASE_NAME]
        )
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 4, 4], [7, 9, 7, 9]], names=["Position A", "Position B"]
        )
        n, a, b, ab = init_confusion_matrix(pos_index, min_gap=2)
        for entry in [n, a, b, ab]:
            self.assertIsInstance(entry, pd.Series)
            self.assertTrue(entry.index.equals(expect_index))
            self.assertTrue(entry.index.equal_levels(expect_index))
            self.assertTrue(np.all(entry == 0))

    def test_clusters(self):
        pos_index = pd.MultiIndex.from_arrays(
            [[2, 4, 7, 9], ["C", "C", "A", "A"]], names=[POS_NAME, BASE_NAME]
        )
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 2, 4, 4, 7], [4, 7, 9, 7, 9, 9]], names=["Position A", "Position B"]
        )
        clusters = ClustHeader(ks=[2]).index
        n, a, b, ab = init_confusion_matrix(pos_index, clusters)
        for entry in [n, a, b, ab]:
            self.assertIsInstance(entry, pd.DataFrame)
            self.assertTrue(entry.index.equals(expect_index))
            self.assertTrue(entry.index.equal_levels(expect_index))
            self.assertTrue(entry.columns.equals(clusters))
            self.assertTrue(np.all(entry == 0))


class TestCalcConfusionMatrix(ut.TestCase):
    def test_no_clusters(self):
        """
        .. code-block:: none

            Pos.
            Rd. 2579
            0   ====
            2   A=C=
            3   =G=T
            5   .G=.
            8   ..==
            9   T-=.
        """
        pattern = RelPattern.from_counts(count_del=False, count_ins=False)
        read_nums = np.array([0, 2, 3, 5, 8, 9])
        pos_index = pd.MultiIndex.from_arrays(
            [[2, 5, 7, 9], ["C", "C", "A", "A"]], names=[POS_NAME, BASE_NAME]
        )
        mutations = {
            2: {SUB_A: np.array([2]), SUB_T: np.array([9])},
            5: {SUB_G: np.array([3, 5]), DELET: np.array([9])},
            7: {SUB_C: np.array([2])},
            9: {SUB_T: np.array([3])},
        }
        seg_end5s = np.array([[2, 7], [2, 5], [2, 5], [5, 5], [7, 1], [2, 5]])
        seg_end3s = np.array([[5, 9], [7, 9], [5, 9], [7, 7], [9, 0], [5, 7]])
        seg_ends_mask = seg_end5s > seg_end3s
        covering_reads = calc_covered_reads_per_pos(
            pos_index, read_nums, seg_end5s, seg_end3s, seg_ends_mask
        )
        mutated_reads = calc_reads_per_pos(pattern, mutations, pos_index)
        # min_gap = 0
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 2, 5, 5, 7], [5, 7, 9, 7, 9, 9]], names=["Position A", "Position B"]
        )
        en = pd.Series([4, 4, 3, 5, 3, 4], expect_index)
        ea = pd.Series([2, 2, 1, 2, 1, 1], expect_index)
        eb = pd.Series([1, 1, 1, 1, 1, 1], expect_index)
        eab = pd.Series([0, 1, 0, 0, 1, 0], expect_index)
        for result, expect in zip(
            calc_confusion_matrix(pos_index, covering_reads, mutated_reads),
            [en, ea, eb, eab],
            strict=True,
        ):
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.index.equals(expect.index))
            self.assertTrue(result.index.equal_levels(expect.index))
            self.assertTrue(result.equals(expect))
        # min_gap = 3
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 5], [7, 9, 9]], names=["Position A", "Position B"]
        )
        en = pd.Series([4, 3, 3], expect_index)
        ea = pd.Series([2, 1, 1], expect_index)
        eb = pd.Series([1, 1, 1], expect_index)
        eab = pd.Series([1, 0, 1], expect_index)
        for result, expect in zip(
            calc_confusion_matrix(pos_index, covering_reads, mutated_reads, min_gap=3),
            [en, ea, eb, eab],
            strict=True,
        ):
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.index.equals(expect.index))
            self.assertTrue(result.index.equal_levels(expect.index))
            self.assertTrue(result.equals(expect))

    def test_clusters(self):
        """
        .. code-block:: none

            Pos.
            Rd. 2579
            0   ====
            2   A=C=
            3   =G=T
            5   .G=.
            8   ..==
            9   T-=.
        """
        pattern = RelPattern.from_counts(count_del=False, count_ins=False)
        read_nums = np.array([0, 2, 3, 5, 8, 9])
        pos_index = pd.MultiIndex.from_arrays(
            [[2, 5, 7, 9], ["C", "C", "A", "A"]], names=[POS_NAME, BASE_NAME]
        )
        clusters = ClustHeader(ks=[2]).index
        read_weights = pd.DataFrame(
            np.stack(
                [
                    np.linspace(0.9, 0.4, read_nums.size),
                    np.linspace(0.1, 0.6, read_nums.size),
                ],
                axis=1,
            ),
            read_nums,
            clusters,
        )
        mutations = {
            2: {SUB_A: np.array([2]), SUB_T: np.array([9])},
            5: {SUB_G: np.array([3, 5]), DELET: np.array([9])},
            7: {SUB_C: np.array([2])},
            9: {SUB_T: np.array([3])},
        }
        seg_end5s = np.array([[2, 7], [2, 5], [2, 5], [5, 5], [7, 1], [2, 5]])
        seg_end3s = np.array([[5, 9], [7, 9], [5, 9], [7, 7], [9, 0], [5, 7]])
        seg_ends_mask = seg_end5s > seg_end3s
        covering_reads = calc_covered_reads_per_pos(
            pos_index, read_nums, seg_end5s, seg_end3s, seg_ends_mask
        )
        mutated_reads = calc_reads_per_pos(pattern, mutations, pos_index)
        expect_index = pd.MultiIndex.from_arrays(
            [[2, 2, 2, 5, 5, 7], [5, 7, 9, 7, 9, 9]], names=["Position A", "Position B"]
        )
        en = pd.DataFrame(
            np.stack(
                [[2.8, 2.8, 2.4, 3.4, 2.4, 2.9], [1.2, 1.2, 0.6, 1.6, 0.6, 1.1]], axis=1
            ),
            expect_index,
            clusters,
        )
        ea = pd.DataFrame(
            np.stack(
                [[1.2, 1.2, 0.8, 1.3, 0.7, 0.8], [0.8, 0.8, 0.2, 0.7, 0.3, 0.2]], axis=1
            ),
            expect_index,
            clusters,
        )
        eb = pd.DataFrame(
            np.stack(
                [[0.7, 0.8, 0.7, 0.8, 0.7, 0.7], [0.3, 0.2, 0.3, 0.2, 0.3, 0.3]], axis=1
            ),
            expect_index,
            clusters,
        )
        eab = pd.DataFrame(
            np.stack(
                [[0.0, 0.8, 0.0, 0.0, 0.7, 0.0], [0.0, 0.2, 0.0, 0.0, 0.3, 0.0]], axis=1
            ),
            expect_index,
            clusters,
        )
        for result, expect in zip(
            calc_confusion_matrix(
                pos_index, covering_reads, mutated_reads, read_weights
            ),
            [en, ea, eb, eab],
            strict=True,
        ):
            self.assertIsInstance(result, pd.DataFrame)
            self.assertTrue(result.index.equals(expect.index))
            self.assertTrue(result.index.equal_levels(expect.index))
            self.assertTrue(result.columns.equals(expect.columns))
            self.assertTrue(np.allclose(result, expect))


class TestCalcConfusionPhi(ut.TestCase):
    def test_series(self):
        n = pd.Series([449491, 1000000, 8, 100, 100, 100])
        a = pd.Series([17259, 200000, 4, 100, 20, 50])
        b = pd.Series([5571, 100000, 2, 100, 20, 50])
        ab = pd.Series([0, 10000, 1, 100, 20, 0])
        expect = pd.Series([-0.022385334, -1 / 12, 0.0, np.nan, 1.0, -1.0])
        result = calc_confusion_phi(n, a, b, ab)
        self.assertTrue(np.allclose(result, expect, equal_nan=True))
        # min_cover = 100
        expect = pd.Series([-0.022385334, -1 / 12, np.nan, np.nan, 1.0, -1.0])
        result = calc_confusion_phi(n, a, b, ab, min_cover=100)
        self.assertTrue(np.allclose(result, expect, equal_nan=True))

    def test_dataframe(self):
        n = pd.DataFrame([[449491, 1000000, 8], [100, 100, 100]])
        a = pd.DataFrame([[17259, 200000, 4], [100, 20, 50]])
        b = pd.DataFrame([[5571, 100000, 2], [100, 20, 50]])
        ab = pd.DataFrame([[0, 10000, 1], [100, 20, 0]])
        expect = pd.DataFrame([[-0.022385334, -1 / 12, 0.0], [np.nan, 1.0, -1.0]])
        result = calc_confusion_phi(n, a, b, ab)
        self.assertTrue(np.allclose(result, expect, equal_nan=True))
        # min_cover = 100
        expect = pd.DataFrame([[-0.022385334, -1 / 12, np.nan], [np.nan, 1.0, -1.0]])
        result = calc_confusion_phi(n, a, b, ab, min_cover=100)
        self.assertTrue(np.allclose(result, expect, equal_nan=True))


def bh_ref_adjusted(p):
    """Reference BH adjusted p-values using scipy, NaN-safe."""
    from scipy.stats import false_discovery_control

    p = np.asarray(p, float)
    is_valid = ~np.isnan(p)
    adj_out = np.full(p.size, np.nan)
    if is_valid.any():
        adj_out[is_valid] = false_discovery_control(p[is_valid], method="bh")
    return adj_out


class TestCalcBhAdjustedPvals(ut.TestCase):
    def test_sorted_input_hand_calc(self):
        # No monotonicity enforcement needed: raw_adj already non-decreasing.
        p = np.array([0.001, 0.009, 0.020, 0.051, 0.200])
        m = 5
        expected = np.array(
            [0.001 * m / 1, 0.009 * m / 2, 0.020 * m / 3, 0.051 * m / 4, 0.200 * m / 5]
        )
        got = calc_bh_adjusted_pvals(p)
        self.assertTrue(np.allclose(got, expected))

    def test_unsorted_input_hand_calc(self):
        # Permutation of the sorted case: each element maps to the same adjusted value.
        p_sorted = np.array([0.001, 0.009, 0.020, 0.051, 0.200])
        m = 5
        expected_sorted = np.array(
            [0.001 * m / 1, 0.009 * m / 2, 0.020 * m / 3, 0.051 * m / 4, 0.200 * m / 5]
        )
        perm = np.array([2, 0, 4, 1, 3])
        p = p_sorted[perm]
        expected = expected_sorted[perm]
        got = calc_bh_adjusted_pvals(p)
        self.assertTrue(np.allclose(got, expected))

    def test_monotonicity_enforcement(self):
        # rank 3 raw_adj (0.05 * 4/3 ≈ 0.0667) > rank 4 raw_adj (0.06 * 4/4 = 0.06),
        # so cummin must pull rank 3 down to 0.06.
        p = np.array([0.001, 0.010, 0.050, 0.060])
        expected = np.array([0.001 * 4 / 1, 0.010 * 4 / 2, 0.06, 0.06])
        got = calc_bh_adjusted_pvals(p)
        self.assertTrue(np.allclose(got, expected))

    def test_nan_propagates_and_valid_entries_computed(self):
        # NaN inputs → NaN outputs; non-NaN entries are computed over the valid subset.
        p = np.array([0.001, np.nan, 0.02, np.nan, 0.9])
        got = calc_bh_adjusted_pvals(p)
        self.assertTrue(np.isnan(got[1]))
        self.assertTrue(np.isnan(got[3]))
        # Valid subset [0.001, 0.02, 0.9] has m=3; adj = [0.003, 0.03, 0.9].
        self.assertAlmostEqual(got[0], 0.003)
        self.assertAlmostEqual(got[2], 0.03)
        self.assertAlmostEqual(got[4], 0.9)

    def test_all_nan(self):
        p = np.array([np.nan, np.nan, np.nan])
        got = calc_bh_adjusted_pvals(p)
        self.assertTrue(np.all(np.isnan(got)))

    def test_single_value(self):
        # m=1: adjusted p-value equals the raw p-value.
        p = np.array([0.03])
        got = calc_bh_adjusted_pvals(p)
        self.assertAlmostEqual(got[0], 0.03)

    def test_order_invariance(self):
        rng = np.random.default_rng(seed=42)
        p = rng.uniform(size=100)
        p[rng.choice(100, size=10, replace=False)] = np.nan
        perm = rng.permutation(100)
        adj1 = calc_bh_adjusted_pvals(p)
        adj2 = calc_bh_adjusted_pvals(p[perm])
        self.assertTrue(np.allclose(adj1[perm], adj2, equal_nan=True))

    def test_pandas_series_preserves_index(self):
        s = pd.Series([0.001, 0.02, 0.3], index=list("abc"))
        got = calc_bh_adjusted_pvals(s)
        self.assertIsInstance(got, pd.Series)
        self.assertListEqual(got.index.tolist(), list("abc"))
        # m=3: adj = [0.001*3, 0.02*3/2, 0.3] = [0.003, 0.03, 0.3]
        self.assertTrue(np.allclose(got.values, [0.003, 0.03, 0.3]))

    def test_2d_ndarray_per_column(self):
        p = np.array(
            [[0.001, 0.009, 0.020, 0.051, 0.200], [0.005, 0.010, 0.010, 0.040, 0.900]]
        ).T
        got = calc_bh_adjusted_pvals(p)
        self.assertEqual(got.shape, p.shape)
        for i in range(p.shape[1]):
            self.assertTrue(np.allclose(got[:, i], calc_bh_adjusted_pvals(p[:, i])))

    def test_pandas_dataframe(self):
        p = pd.DataFrame({"a": [0.001, 0.009, 0.020], "b": [0.010, 0.050, 0.900]})
        got = calc_bh_adjusted_pvals(p)
        self.assertIsInstance(got, pd.DataFrame)
        self.assertEqual(got.shape, p.shape)
        self.assertListEqual(got.columns.tolist(), ["a", "b"])
        self.assertListEqual(got.index.tolist(), p.index.tolist())

    def test_against_reference_implementation(self):
        rng = np.random.default_rng(seed=0)
        p = rng.uniform(size=500)
        p[rng.choice(500, size=20, replace=False)] = np.nan
        got = calc_bh_adjusted_pvals(p)
        ref = bh_ref_adjusted(p)
        self.assertTrue(np.allclose(got, ref, equal_nan=True))

    def test_invalid_pvalues_raise(self):
        with self.assertRaises(ValueError):
            calc_bh_adjusted_pvals(np.array([0.1, 2.0, 0.3]))
        with self.assertRaises(ValueError):
            calc_bh_adjusted_pvals(np.array([-0.1, 0.5]))


if __name__ == "__main__":
    ut.main(verbosity=2)
