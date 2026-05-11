import unittest
import pandas as pd
import numpy as np
from seismicrna.core.rna.roc import _compute_fpr_tpr, compute_roc_curve, compute_auc_roc


class TestComputeFprTpr(unittest.TestCase):
    def test_perfect_prediction(self):
        is_paired = pd.Series([True, True, False, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.5, 0.0, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 1.0, 0.5, 0.0])
        ))

    def test_all_false_positives(self):
        is_paired = pd.Series([False, False, True, True], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 1.0, 1.0, 0.5, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 0.5, 0.0, 0.0, 0.0])
        ))

    def test_mixed_prediction(self):
        is_paired = pd.Series([True, False, True, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.5, 0.5, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 0.5, 0.5, 0.0])
        ))

    def test_empty_prediction(self):
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired)
        self.assertTrue(np.array_equal(fpr,
                                       np.array([np.nan]),
                                       equal_nan=True))
        self.assertTrue(np.array_equal(tpr,
                                       np.array([np.nan]),
                                       equal_nan=True))

    def test_invalid_dtype(self):
        is_paired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            _compute_fpr_tpr(is_paired)

    def test_all_positions_unpaired(self):
        is_paired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([np.nan, np.nan, np.nan, np.nan, np.nan]),
            equal_nan=True
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 0.75, 0.5, 0.25, 0.0])
        ))

    def test_mixed_paired_unpaired(self):
        is_paired = pd.Series([True, False, True, False], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, True, False, True], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.5, 0.5, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 0.5, 0.5, 0.0])
        ))

    def test_neither_paired_nor_unpaired(self):
        is_paired = pd.Series([True, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, False, True, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.0, 0.0, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 1.0, 0.0, 0.0])
        ))

    def test_no_positions_unpaired(self):
        is_paired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.75, 0.5, 0.25, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([np.nan, np.nan, np.nan, np.nan, np.nan]),
            equal_nan=True
        ))

    def test_empty_prediction_with_unpaired(self):
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        is_unpaired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        fpr, tpr = _compute_fpr_tpr(is_paired, is_unpaired)
        self.assertTrue(np.array_equal(fpr,
                                       np.array([np.nan]),
                                       equal_nan=True))
        self.assertTrue(np.array_equal(tpr,
                                       np.array([np.nan]),
                                       equal_nan=True))

    def test_invalid_dtype_for_unpaired(self):
        is_paired = pd.Series([True, False, True, False], dtype=bool, index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            _compute_fpr_tpr(is_paired, is_unpaired)


class TestComputeRocCurve(unittest.TestCase):
    def test_perfect_prediction(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, False, True, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.5, 0.0, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 1.0, 0.5, 0.0])
        ))

    def test_all_false_positives(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, True, False, True], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 1.0, 1.0, 0.5, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 0.5, 0.0, 0.0, 0.0])
        ))

    def test_mixed_prediction(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, True, True, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 1.0, 0.5, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 0.5, 0.5, 0.5, 0.0])
        ))

    def test_empty_prediction(self):
        profile = pd.Series([], dtype=float, index=pd.Index([], name='Position'))
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired)
        self.assertTrue(np.array_equal(fpr,
                                       np.array([np.nan]),
                                       equal_nan=True))
        self.assertTrue(np.array_equal(tpr,
                                       np.array([np.nan]),
                                       equal_nan=True))

    def test_invalid_dtype(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            compute_roc_curve(profile, is_paired)

    def test_all_positions_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([np.nan, np.nan, np.nan, np.nan, np.nan]),
            equal_nan=True
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 0.75, 0.5, 0.25, 0.0])
        ))

    def test_neither_paired_nor_unpaired(self):
        profile = pd.Series([0.1, 0.8, 0.9, 0.2], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, False, False, True], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, True, False, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 1.0, 0.0, 0.0, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([1.0, 1.0, 1.0, 0.0, 0.0])
        ))

    def test_no_positions_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired, is_unpaired)
        self.assertTrue(np.array_equal(
            fpr,
            np.array([1.0, 0.75, 0.5, 0.25, 0.0])
        ))
        self.assertTrue(np.array_equal(
            tpr,
            np.array([np.nan, np.nan, np.nan, np.nan, np.nan]),
            equal_nan=True
        ))

    def test_empty_prediction_with_unpaired(self):
        profile = pd.Series([], dtype=float, index=pd.Index([], name='Position'))
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        is_unpaired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        fpr, tpr = compute_roc_curve(profile, is_paired, is_unpaired)
        self.assertTrue(np.array_equal(fpr,
                                       np.array([np.nan]),
                                       equal_nan=True))
        self.assertTrue(np.array_equal(tpr,
                                       np.array([np.nan]),
                                       equal_nan=True))

    def test_invalid_dtype_for_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, False, True, False], dtype=bool, index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            compute_roc_curve(profile, is_paired, is_unpaired)


class TestComputeAucRoc(unittest.TestCase):
    def test_perfect_prediction(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, False, True, False], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired)
        self.assertEqual(auc, 1.0)

    def test_all_false_positives(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, True, False, True], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired)
        self.assertEqual(auc, 0.0)

    def test_mixed_prediction(self):
        profile = pd.Series([0.1, 0.8, 0.2, 0.9], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, True, True, False], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired)
        self.assertEqual(auc, 0.5)

    def test_empty_prediction(self):
        profile = pd.Series([], dtype=float, index=pd.Index([], name='Position'))
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        auc = compute_auc_roc(profile, is_paired)
        self.assertTrue(np.isnan(auc))

    def test_invalid_dtype(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            compute_auc_roc(profile, is_paired)

    def test_all_positions_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired, is_unpaired)
        self.assertTrue(np.isnan(auc))

    def test_neither_paired_nor_unpaired(self):
        profile = pd.Series([0.1, 0.8, 0.9, 0.2], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([False, False, False, True], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, True, False, False], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired, is_unpaired)
        self.assertEqual(auc, 1.0)

    def test_no_positions_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, True, True, True], index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([False, False, False, False], index=pd.Index(range(1, 5), name='Position'))
        auc = compute_auc_roc(profile, is_paired, is_unpaired)
        self.assertTrue(np.isnan(auc))

    def test_empty_prediction_with_unpaired(self):
        profile = pd.Series([], dtype=float, index=pd.Index([], name='Position'))
        is_paired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        is_unpaired = pd.Series([], dtype=bool, index=pd.Index([], name='Position'))
        auc = compute_auc_roc(profile, is_paired, is_unpaired)
        self.assertTrue(np.isnan(auc))

    def test_invalid_dtype_for_unpaired(self):
        profile = pd.Series([0.1, 0.2, 0.3, 0.4], index=pd.Index(range(1, 5), name='Position'))
        is_paired = pd.Series([True, False, True, False], dtype=bool, index=pd.Index(range(1, 5), name='Position'))
        is_unpaired = pd.Series([1, 0, 1, 0], dtype=int, index=pd.Index(range(1, 5), name='Position'))
        with self.assertRaises(AssertionError):
            compute_auc_roc(profile, is_paired, is_unpaired)


if __name__ == "__main__":
    unittest.main()
