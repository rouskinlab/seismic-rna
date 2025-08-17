import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.confusion import (init_confusion_matrix,
                                             calc_confusion_matrix,
                                             calc_confusion_phi,
_count_intersection,
                                             label_significant_pvals)
from seismicrna.core.rel.code import DELET, SUB_A, SUB_C, SUB_G, SUB_T
from seismicrna.core.rel.pattern import RelPattern
from seismicrna.core.seq.region import POS_NAME, BASE_NAME

rng = np.random.default_rng()


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
        x = np.unique(rng.integers(100, size=100))
        y = np.unique(rng.integers(100, size=100))
        count = np.intersect1d(x, y).size
        self.assertEqual(_count_intersection(x, y), count)
        self.assertEqual(_count_intersection(y, x), count)


class TestInitConfusionMatrix(ut.TestCase):

    def test_no_clusters(self):
        pos_index = pd.MultiIndex.from_arrays([[2, 4, 7, 9],
                                               ["C", "C", "A", "A"]],
                                              names=[POS_NAME,
                                                     BASE_NAME])
        expect_index = pd.MultiIndex.from_arrays([[2, 2, 2, 4, 4, 7],
                                                  [4, 7, 9, 7, 9, 9]],
                                                 names=["Position A",
                                                        "Position B"])
        n, a, b, ab = init_confusion_matrix(pos_index)
        for entry in [n, a, b, ab]:
            self.assertIsInstance(entry, pd.Series)
            self.assertTrue(entry.index.equals(expect_index))
            self.assertTrue(entry.index.equal_levels(expect_index))
            self.assertTrue(np.all(entry == 0))


class TestCalcConfusionMatrix(ut.TestCase):

    def test_no_clusters(self):
        """
            Pos.
        Rd. 2479
        0   ====
        2   A=C=
        3   =G=T
        5   .G=.
        8   ..==
        9   T-=.
        """
        pattern = RelPattern.from_counts(count_del=False, count_ins=False)
        read_nums = np.array([0, 2, 3, 5, 8, 9])
        pos_index = pd.MultiIndex.from_arrays([[2, 4, 7, 9],
                                               ["C", "C", "A", "A"]],
                                              names=[POS_NAME, BASE_NAME])
        mutations = {2: {SUB_A: np.array([2]),
                         SUB_T: np.array([9])},
                     4: {SUB_G: np.array([3, 5]),
                         DELET: np.array([9])},
                     7: {SUB_C: np.array([2])},
                     9: {SUB_T: np.array([3])}}
        seg_end5s = np.array([[2, 7],
                              [2, 4],
                              [2, 4],
                              [4, 4],
                              [7, 1],
                              [2, 4]])
        seg_end3s = np.array([[4, 9],
                              [7, 9],
                              [4, 9],
                              [7, 7],
                              [9, 0],
                              [4, 7]])
        seg_ends_mask = seg_end5s > seg_end3s
        expect_index = pd.MultiIndex.from_arrays([[2, 2, 2, 4, 4, 7],
                                                  [4, 7, 9, 7, 9, 9]],
                                                 names=["Position A",
                                                        "Position B"])
        en = pd.Series([4, 4, 3, 5, 3, 4], expect_index)
        ea = pd.Series([2, 2, 1, 2, 1, 1], expect_index)
        eb = pd.Series([1, 1, 1, 1, 1, 1], expect_index)
        eab = pd.Series([0, 1, 0, 0, 1, 0], expect_index)
        for result, expect in zip(
                calc_confusion_matrix(pattern,
                                      mutations,
                                      pos_index,
                                      read_nums,
                                      seg_end5s,
                                      seg_end3s,
                                      seg_ends_mask),
                [en, ea, eb, eab],
                strict=True
        ):
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.index.equals(expect.index))
            self.assertTrue(result.index.equal_levels(expect.index))
            self.assertTrue(result.equals(expect))


class TestCalcConfusionPhi(ut.TestCase):

    def test_series(self):
        n = pd.Series([449491,
                       1000000,
                       8,
                       100,
                       100,
                       100])
        a = pd.Series([17259,
                       200000,
                       4,
                       100,
                       20,
                       50])
        b = pd.Series([5571,
                       100000,
                       2,
                       100,
                       20,
                       50])
        ab = pd.Series([0,
                        10000,
                        1,
                        100,
                        20,
                        0])
        expect = pd.Series([-0.022385334,
                            -1 / 12,
                            0.0,
                            np.nan,
                            1.0,
                            -1.0])
        result = calc_confusion_phi(n, a, b, ab)
        self.assertTrue(np.allclose(result,
                                    expect,
                                    equal_nan=True))


def bh_ref_mask(p, alpha):
    """Reference BH step-up (finite p only), returns boolean mask."""
    p = np.asarray(p, float)
    finite = np.isfinite(p)
    pf = p[finite]
    m = pf.size
    mask = np.zeros_like(p, dtype=bool)
    if m == 0:
        return mask
    order = np.argsort(pf, kind="mergesort")
    pf_sorted = pf[order]
    ranks = np.arange(1, m + 1)
    thresh = (ranks / m) * alpha
    passed = pf_sorted <= thresh
    if not np.any(passed):
        return mask
    k = np.max(np.where(passed)[0])  # 0-based
    cutoff = pf_sorted[k]
    mask[finite] = pf <= cutoff
    return mask


class TestBHLabeling(ut.TestCase):

    def test_basic_hand_calcs(self):
        # Simple, hand-checkable example
        p = np.array([0.001, 0.009, 0.020, 0.051, 0.20])
        alpha = 0.05
        got = label_significant_pvals(p, alpha)
        want = np.array([True, True, True, False, False])
        self.assertTrue(np.array_equal(got, want))

    def test_ties_including_cutoff(self):
        # Ties at cutoff should all be rejected
        p = np.array([0.005, 0.010, 0.010, 0.04, 0.9])
        alpha = 0.05
        got = label_significant_pvals(p, alpha)
        want = np.array([True, True, True, True, False])
        self.assertTrue(np.array_equal(got, want))

    def test_nan_are_ignored_and_false(self):
        p = np.array([0.001, np.nan, 0.02, np.nan, 0.9, np.nan])
        alpha = 0.05
        got = label_significant_pvals(p, alpha)
        # finite subset is [0.001, 0.02, 0.9]; BH rejects first two at 5%
        want = np.array([True, False, True, False, False, False])
        self.assertTrue(np.array_equal(got, want))

    def test_order_invariance(self):
        p = rng.uniform(size=200)
        p[[5, 17, 33]] = np.nan
        alpha = 0.1
        m1 = label_significant_pvals(p, alpha)
        m2 = label_significant_pvals(p[::-1], alpha)[::-1]
        self.assertTrue(np.array_equal(m1, m2))

    def test_alpha_validation(self):
        p = np.array([0.1, 0.2, 0.3])
        with self.assertRaises(ValueError):
            label_significant_pvals(p, 0.0)
        with self.assertRaises(ValueError):
            label_significant_pvals(p, 1.0)
        with self.assertRaises(ValueError):
            label_significant_pvals(p, -0.1)

    def test_dtype_and_shape_preserved(self):
        for dtype in (np.float64, np.float32):
            with self.subTest(dtype=dtype):
                p = np.array([0.01, 0.2, np.nan, 0.03], dtype=dtype)
                mask = label_significant_pvals(p, 0.05)
                self.assertEqual(mask.shape, p.shape)
                self.assertEqual(mask.dtype, bool)

    def test_pandas_series_roundtrip(self):
        s = pd.Series([0.01, 0.2, np.nan, 0.03], index=list("abcd"))
        m = label_significant_pvals(s, 0.05)
        self.assertIsInstance(m, pd.Series)
        self.assertListEqual(m.index.tolist(), list("abcd"))

    def test_monotonicity_property(self):
        # Lowering any p-value cannot reduce the number of rejections
        p = rng.uniform(size=200)
        alpha = 0.05
        base = label_significant_pvals(p, alpha)
        idx = rng.choice(len(p), size=50, replace=False)
        p2 = p.copy()
        p2[idx] = np.maximum(0.0, p2[idx] * 0.5)
        base2 = label_significant_pvals(p2, alpha)
        self.assertGreaterEqual(base2.sum(), base.sum())

    def test_all_or_none(self):
        # None significant
        p = np.array([0.8, 0.6, 0.9, 0.7])
        self.assertFalse(label_significant_pvals(p, 0.05).any())
        # All significant (with a generous alpha)
        p = np.array([1e-10, 1e-9, 1e-8, 1e-7])
        self.assertTrue(label_significant_pvals(p, 0.2).all())

    def test_against_reference_implementation(self):
        p = rng.uniform(size=500)
        p[rng.choice(500, size=20, replace=False)] = np.nan
        alpha = 0.07
        got = label_significant_pvals(p, alpha)
        ref = bh_ref_mask(p, alpha)
        self.assertTrue(np.array_equal(got, ref))

    def test_rejects_values_failing_own_rank_due_to_global_cutoff(self):
        # Construct p-values where p_(1) > (1/m)*alpha, but a larger k passes,
        # so BH rejects ALL p <= p_(k) (including that first one).
        # Sorted: [0.013, 0.026, 0.026, 0.026]
        # Thresholds for alpha=0.05, m=4: [0.0125, 0.025, 0.0375, 0.05]
        # p_(1)=0.013 > 0.0125 (fails own threshold),
        # p_(4)=0.026 <= 0.05 -> k=4, cutoff = 0.026 -> reject all 4.
        p = np.array([0.026, 0.013, 0.051, 0.026, 0.025])  # intentionally unsorted
        alpha = 0.05

        mask = label_significant_pvals(p, alpha)
        self.assertTrue(np.array_equal(
            mask,
            np.array([True, True, False, True, True])
        ))

        # Sanity check: the smallest p actually fails its own rank threshold
        p_sorted = np.sort(p)
        m = p_sorted.size
        t1 = alpha * 1 / m
        self.assertGreater(
            p_sorted[0], t1,
            "Smallest p-value should fail its own (1/m)*alpha threshold in this test."
        )

    def test_invalid_pvalues_raise(self):
        # Values outside [0, 1] must raise
        with self.assertRaises(ValueError):
            label_significant_pvals(np.array([0.1, 2.0, 0.3]), 0.05)
        with self.assertRaises(ValueError):
            label_significant_pvals(np.array([-1.0, 0.2, 0.5]), 0.05)
        with self.assertRaises(ValueError):
            label_significant_pvals(np.array([np.inf, 0.2, 0.5]), 0.05)
        with self.assertRaises(ValueError):
            label_significant_pvals(np.array([-np.inf, 0.2, 0.5]), 0.05)

        # Also check pandas Series input
        with self.assertRaises(ValueError):
            s = pd.Series([0.0, 1.1, 0.2], index=list("abc"))
            label_significant_pvals(s, 0.05)


if __name__ == "__main__":
    ut.main(verbosity=2)
