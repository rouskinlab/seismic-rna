import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.mu import calc_ranks

rng = np.random.default_rng()


class TestCalcRanks(ut.TestCase):

    def test_array0d(self):
        mus = np.array(0.52)
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_ranks,
                               mus)

    def test_array1d(self):
        mus = np.array([0.49, 0.17, 0.24, 0.90, 0.47, 0.50, 0.22, 0.14, 0.18])
        expect = np.array([6, 1, 4, 8, 5, 7, 3, 0, 2])
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_array2d(self):
        mus = np.array([[0.04, 0.61],
                        [0.60, 0.59],
                        [0.73, 0.81],
                        [0.22, 0.44],
                        [0.88, 0.78],
                        [0.48, 0.58]])
        expect = np.array([[0, 3],
                           [3, 2],
                           [4, 5],
                           [1, 0],
                           [5, 4],
                           [2, 1]])
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_array3d(self):
        mus = np.array([[[0.59, 0.23],
                         [0.67, 0.94],
                         [0.93, 0.73]],
                        [[0.15, 0.06],
                         [0.08, 0.71],
                         [0.89, 0.26]],
                        [[0.53, 0.03],
                         [0.34, 0.76],
                         [0.39, 0.52]],
                        [[0.89, 0.21],
                         [0.43, 0.50],
                         [0.23, 0.66]]])
        expect = np.array([[[2, 3],
                            [3, 3],
                            [3, 3]],
                           [[0, 1],
                            [0, 1],
                            [2, 0]],
                           [[1, 0],
                            [1, 2],
                            [1, 1]],
                           [[3, 2],
                            [2, 0],
                            [0, 2]]])
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_series(self):
        mus = pd.Series([0.38, 0.50, 0.51, 0.64, 0.75, 0.60, 0.18],
                        index=[1, 2, 4, 5, 6, 7, 9])
        expect = pd.Series([1, 2, 3, 5, 6, 4, 0], index=mus.index)
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, pd.Series)
        self.assertTrue(ranks.equals(expect))

    def test_dataframe(self):
        mus = pd.DataFrame([[0.23, 0.16, 0.72],
                            [0.63, 0.10, 0.73],
                            [0.29, 0.19, 0.69],
                            [0.14, 0.09, 0.24],
                            [0.76, 0.55, 0.34]],
                           index=[2, 3, 5, 7, 8],
                           columns=["x", "y", "z"])
        expect = pd.DataFrame([[1, 2, 3],
                               [3, 1, 4],
                               [2, 3, 2],
                               [0, 0, 0],
                               [4, 4, 1]],
                              mus.index,
                              mus.columns)
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, pd.DataFrame)
        self.assertTrue(ranks.equals(expect))

    def test_nan(self):
        mus = np.array([0.49, 0.17, 0.24, 0.90, np.nan, 0.50, 0.22, 0.14, 0.18])
        expect = np.array([5, 1, 4, 7, 6, 3, 0, 2])
        ranks = calc_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))


if __name__ == "__main__":
    ut.main(verbosity=2)
