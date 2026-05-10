import unittest as ut
import warnings
from itertools import permutations
from types import SimpleNamespace

import numpy as np

from seismicrna.cluster.emk import (assign_clusterings,
                                     calc_mean_arcsine_distance_clusters,
                                     calc_mean_pearson_clusters,
                                     get_common_k,
                                     sort_runs)
from seismicrna.core.array import calc_inverse
from seismicrna.core.error import NoDataError, InconsistentValueError


def _mock_run(k: int, log_like: float):
    """ Minimal stand-in for EMRun with just the attributes used by
    get_common_k / sort_runs. """
    return SimpleNamespace(k=k, log_like=log_like)


class TestAssignClusterings(ut.TestCase):

    def compare_result(self, x: np.ndarray, y: np.ndarray, expect: np.ndarray):
        rows, cols = assign_clusterings(x, y)
        self.assertTupleEqual(rows.shape, cols.shape)
        self.assertTrue(np.array_equal(rows, np.arange(rows.size)))
        self.assertTrue(np.array_equal(cols, expect))

    def test_0_clusters(self):
        for n in range(5):
            x = np.empty((n, 0))
            y = np.empty((n, 0))
            expect = np.empty(0, dtype=int)
            self.compare_result(x, y, expect)

    def test_0_positions(self):
        for n in range(5):
            x = np.empty((0, n))
            y = np.empty((0, n))
            expect = np.arange(n)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                with np.errstate(invalid="ignore"):
                    self.compare_result(x, y, expect)

    def test_1_cluster(self):
        rng = np.random.default_rng(seed=0)
        for n in range(1, 10):
            x = rng.random((n, 1))
            y = rng.random((n, 1))
            expect = np.array([0])
            self.compare_result(x, y, expect)

    def test_more_clusters(self):
        rng = np.random.default_rng(seed=0)
        for ncls in range(2, 5):
            for npos in range(1, 8):
                with self.subTest(ncls=ncls, npos=npos):
                    x = rng.random((npos, ncls))
                    for permutation in permutations(range(ncls)):
                        cols = np.array(permutation)
                        y = x[:, cols]
                        self.compare_result(x, y, calc_inverse(cols))
                        self.compare_result(y, x, cols)


class TestGetCommonK(ut.TestCase):

    def test_empty_raises(self):
        with self.assertRaises(NoDataError):
            get_common_k([])

    def test_single_run(self):
        self.assertEqual(get_common_k([_mock_run(3, -10.)]), 3)

    def test_multiple_runs_same_k(self):
        runs = [_mock_run(2, -5.), _mock_run(2, -7.), _mock_run(2, -6.)]
        self.assertEqual(get_common_k(runs), 2)

    def test_mixed_ks_raises(self):
        runs = [_mock_run(1, -5.), _mock_run(2, -6.)]
        with self.assertRaises(InconsistentValueError):
            get_common_k(runs)

    def test_returns_correct_k(self):
        for k in [1, 2, 5, 10]:
            with self.subTest(k=k):
                runs = [_mock_run(k, float(-i)) for i in range(3)]
                self.assertEqual(get_common_k(runs), k)


class TestSortRuns(ut.TestCase):

    def test_empty_list(self):
        self.assertEqual(sort_runs([]), [])

    def test_single_run(self):
        run = _mock_run(1, -5.)
        result = sort_runs([run])
        self.assertEqual(result, [run])

    def test_already_sorted(self):
        runs = [_mock_run(2, -3.), _mock_run(2, -5.), _mock_run(2, -10.)]
        result = sort_runs(runs)
        log_likes = [r.log_like for r in result]
        self.assertEqual(log_likes, sorted(log_likes, reverse=True))

    def test_reverse_order(self):
        runs = [_mock_run(1, -10.), _mock_run(1, -5.), _mock_run(1, -3.)]
        result = sort_runs(runs)
        log_likes = [r.log_like for r in result]
        self.assertEqual(log_likes, [-3., -5., -10.])

    def test_mixed_ks_raises(self):
        runs = [_mock_run(1, -3.), _mock_run(2, -5.)]
        with self.assertRaises(InconsistentValueError):
            sort_runs(runs)

    def test_best_run_first(self):
        runs = [_mock_run(3, float(-i)) for i in [7, 2, 9, 4]]
        result = sort_runs(runs)
        self.assertEqual(result[0].log_like, -2.)


class TestCalcMeanArcsineDistanceClusters(ut.TestCase):

    def test_identical_clusters_distance_zero(self):
        mus = np.array([[0.1, 0.4],
                        [0.2, 0.5],
                        [0.3, 0.6]])
        result = calc_mean_arcsine_distance_clusters(mus, mus.copy())
        self.assertAlmostEqual(result, 0., places=10)

    def test_single_cluster(self):
        mus = np.array([[0.1], [0.3], [0.5]])
        result = calc_mean_arcsine_distance_clusters(mus, mus.copy())
        self.assertAlmostEqual(result, 0., places=10)

    def test_symmetry(self):
        rng = np.random.default_rng(42)
        mus1 = rng.uniform(0.01, 0.99, (5, 2))
        mus2 = rng.uniform(0.01, 0.99, (5, 2))
        d12 = calc_mean_arcsine_distance_clusters(mus1, mus2)
        d21 = calc_mean_arcsine_distance_clusters(mus2, mus1)
        self.assertAlmostEqual(d12, d21, places=10)

    def test_permuted_columns_same_result(self):
        mus1 = np.array([[0.1, 0.8],
                         [0.2, 0.7]])
        mus2_permuted = mus1[:, [1, 0]]  # swap cluster columns
        result = calc_mean_arcsine_distance_clusters(mus1, mus2_permuted)
        self.assertAlmostEqual(result, 0., places=10)

    def test_nonnegative(self):
        rng = np.random.default_rng(7)
        mus1 = rng.uniform(0.01, 0.99, (6, 3))
        mus2 = rng.uniform(0.01, 0.99, (6, 3))
        result = calc_mean_arcsine_distance_clusters(mus1, mus2)
        self.assertGreaterEqual(result, 0.)


class TestCalcMeanPearsonClusters(ut.TestCase):

    def test_identical_clusters_correlation_one(self):
        mus = np.array([[0.1, 0.4],
                        [0.2, 0.5],
                        [0.3, 0.6]])
        result = calc_mean_pearson_clusters(mus, mus.copy())
        self.assertAlmostEqual(result, 1., places=10)

    def test_single_cluster(self):
        mus = np.array([[0.1], [0.3], [0.5]])
        result = calc_mean_pearson_clusters(mus, mus.copy())
        self.assertAlmostEqual(result, 1., places=10)

    def test_linearly_related_clusters_one(self):
        # Pearson is invariant to positive linear transforms
        mus1 = np.array([[0.1, 0.5], [0.3, 0.3], [0.5, 0.1]])
        # Compress towards centre so all values remain in (0, 1)
        mus2 = 0.5 + 0.5 * (mus1 - 0.5)
        result = calc_mean_pearson_clusters(mus1, mus2)
        self.assertAlmostEqual(result, 1., places=10)

    def test_permuted_columns_same_result(self):
        mus1 = np.array([[0.1, 0.8],
                         [0.2, 0.7],
                         [0.3, 0.6]])
        mus2_permuted = mus1[:, [1, 0]]
        result = calc_mean_pearson_clusters(mus1, mus2_permuted)
        self.assertAlmostEqual(result, 1., places=10)

    def test_returns_float(self):
        mus = np.array([[0.1, 0.4], [0.2, 0.5], [0.3, 0.6]])
        result = calc_mean_pearson_clusters(mus, mus)
        self.assertIsInstance(result, float)


if __name__ == "__main__":
    ut.main(verbosity=2)
