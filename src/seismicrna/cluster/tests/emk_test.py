import unittest as ut
import warnings
from itertools import permutations

import numpy as np

from seismicrna.cluster.emk import assign_clusterings
from seismicrna.core.array import calc_inverse

rng = np.random.default_rng()


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
        for n in range(1, 10):
            x = rng.random((n, 1))
            y = rng.random((n, 1))
            expect = np.array([0])
            self.compare_result(x, y, expect)

    def test_more_clusters(self):
        for ncls in range(2, 5):
            for npos in range(1, 8):
                with self.subTest(ncls=ncls, npos=npos):
                    x = rng.random((npos, ncls))
                    for permutation in permutations(range(ncls)):
                        cols = np.array(permutation)
                        y = x[:, cols]
                        self.compare_result(x, y, calc_inverse(cols))
                        self.compare_result(y, x, cols)


if __name__ == "__main__":
    ut.main(verbosity=2)
