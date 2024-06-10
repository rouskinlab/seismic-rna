import unittest as ut
import warnings
from itertools import permutations

import numpy as np

from seismicrna.cluster.compare import assign_clusterings
from seismicrna.core.array import calc_inverse

rng = np.random.default_rng()


class TestAssignClusterings(ut.TestCase):

    def test_0_clusters(self):
        for n in range(5):
            x = np.empty((n, 0))
            y = np.empty((n, 0))
            self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                           np.empty(0, dtype=int)))

    def test_0_positions(self):
        for n in range(5):
            x = np.empty((0, n))
            y = np.empty((0, n))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                with np.errstate(invalid="ignore"):
                    self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                                   np.arange(n)))

    def test_1_cluster(self):
        for n in range(1, 10):
            x = rng.random((n, 1))
            y = rng.random((n, 1))
            self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                           np.array([0])))

    def test_more_clusters(self):
        for ncls in range(2, 5):
            for npos in range(2, 8):
                with self.subTest(ncls=ncls, npos=npos):
                    x = rng.random((npos, ncls))
                    for permutation in permutations(range(ncls)):
                        cols = np.array(permutation)
                        y = x[:, cols]
                        self.assertTrue(np.array_equal(
                            assign_clusterings(x, y),
                            calc_inverse(cols)
                        ))
                        self.assertTrue(np.array_equal(
                            assign_clusterings(y, x),
                            cols
                        ))


if __name__ == "__main__":
    ut.main()
