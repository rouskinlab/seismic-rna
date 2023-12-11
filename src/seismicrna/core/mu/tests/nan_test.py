import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.mu.nan import any_nan, no_nan, removes_nan

rng = np.random.default_rng()


class TestAnyNaN(ut.TestCase):

    def test_array0d(self):
        for mus, expect in [(np.array(np.nan), True),
                            (np.array(0.0), False)]:
            a = any_nan(mus)
            self.assertEqual(a, expect)

    def test_array1d(self):
        mus = np.array([0.0, 0.1, np.nan, 0.3, 0.4])
        expect = np.array([False, False, True, False, False])
        a = any_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_array2d(self):
        mus = np.array([[0.0, 0.5],
                        [0.1, np.nan],
                        [np.nan, np.nan],
                        [0.3, 0.8],
                        [np.nan, 0.9]])
        expect = np.array([False, True, True, False, True])
        a = any_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_array3d(self):
        mus = np.array([[[np.nan, 0.1],
                         [0.2, 0.3]],
                        [[0.0, np.nan],
                         [0.2, 0.3]],
                        [[0.0, 0.1],
                         [0.2, 0.3]],
                        [[0.0, 0.1],
                         [np.nan, 0.3]],
                        [[0.0, 0.1],
                         [0.2, np.nan]]])
        expect = np.array([True, True, False, True, True])
        a = any_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_series(self):
        mus = pd.Series([0.0, 0.1, np.nan, 0.3, 0.4],
                        index=[1, 2, 6, 7, 8])
        expect = pd.Series([False, False, True, False, False],
                           index=[1, 2, 6, 7, 8])
        a = any_nan(mus)
        self.assertIsInstance(a, pd.Series)
        self.assertTrue(a.equals(expect))

    def test_dataframe(self):
        mus = pd.DataFrame([[0.0, 0.5],
                            [0.1, np.nan],
                            [np.nan, np.nan],
                            [0.3, 0.8],
                            [np.nan, 0.9]],
                           index=[1, 2, 6, 7, 8],
                           columns=["x", "y"])
        expect = pd.Series([False, True, True, False, True],
                           index=[1, 2, 6, 7, 8])
        a = any_nan(mus)
        self.assertIsInstance(a, pd.Series)
        self.assertTrue(a.equals(expect))


class TestNoNaN(ut.TestCase):

    def test_array1d(self):
        mus = np.array([0.0, 0.1, np.nan, 0.3, 0.4])
        expect = np.array([True, True, False, True, True])
        a = no_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_array2d(self):
        mus = np.array([[0.0, 0.5],
                        [0.1, np.nan],
                        [np.nan, np.nan],
                        [0.3, 0.8],
                        [np.nan, 0.9]])
        expect = np.array([True, False, False, True, False])
        a = no_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_series(self):
        mus = pd.Series([0.0, 0.1, np.nan, 0.3, 0.4],
                        index=[1, 2, 6, 7, 8])
        expect = pd.Series([True, True, False, True, True],
                           index=[1, 2, 6, 7, 8])
        a = no_nan(mus)
        self.assertIsInstance(a, pd.Series)
        self.assertTrue(a.equals(expect))

    def test_dataframe(self):
        mus = pd.DataFrame([[0.0, 0.5],
                            [0.1, np.nan],
                            [np.nan, np.nan],
                            [0.3, 0.8],
                            [np.nan, 0.9]],
                           index=[1, 2, 6, 7, 8],
                           columns=["x", "y"])
        expect = pd.Series([True, False, False, True, False],
                           index=[1, 2, 6, 7, 8])
        a = no_nan(mus)
        self.assertIsInstance(a, pd.Series)
        self.assertTrue(a.equals(expect))


class TestWithoutNaNs(ut.TestCase):

    def test_array1d(self):
        mus = np.array([0.0, 0.1, np.nan, 0.3, 0.4])
        expect = np.array([0.0, 0.1, 0.3, 0.4])
        a, = removes_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_array2d(self):
        mus = np.array([[0.0, 0.5],
                        [0.1, np.nan],
                        [np.nan, np.nan],
                        [0.3, 0.8],
                        [np.nan, 0.9]])
        expect = np.array([[0.0, 0.5],
                           [0.3, 0.8]])
        a, = removes_nan(mus)
        self.assertIsInstance(a, np.ndarray)
        self.assertTrue(np.array_equal(a, expect))

    def test_series(self):
        mus = pd.Series([0.0, 0.1, np.nan, 0.3, 0.4],
                        index=[1, 2, 6, 7, 8])
        expect = pd.Series([0.0, 0.1, 0.3, 0.4],
                           index=[1, 2, 7, 8])
        a, = removes_nan(mus)
        self.assertIsInstance(a, pd.Series)
        self.assertTrue(a.equals(expect))

    def test_dataframe(self):
        mus = pd.DataFrame([[0.0, 0.5],
                            [0.1, np.nan],
                            [np.nan, np.nan],
                            [0.3, 0.8],
                            [np.nan, 0.9]],
                           index=[1, 2, 6, 7, 8],
                           columns=["x", "y"])
        expect = pd.DataFrame([[0.0, 0.5],
                               [0.3, 0.8]],
                              index=[1, 7],
                              columns=["x", "y"])
        a, = removes_nan(mus)
        self.assertIsInstance(a, pd.DataFrame)
        self.assertTrue(a.equals(expect))

    def test_array_series_dataframe(self):
        n = np.array([0.0, np.nan, 0.2, 0.3, 0.4])
        s = pd.Series([0.5, 0.6, 0.7, np.nan, 0.9],
                      index=[0, 1, 2, 3, 4])
        f = pd.DataFrame([[0.0, 0.5],
                          [0.1, 0.6],
                          [np.nan, 0.7],
                          [0.3, 0.8],
                          [0.4, 0.9]],
                         index=[5, 6, 7, 8, 9],
                         columns=["x", "y"])
        na, sa, fa = removes_nan(n, s, f)
        self.assertIsInstance(na, np.ndarray)
        self.assertTrue(np.array_equal(na, np.array([0.0, 0.4])))
        self.assertIsInstance(sa, pd.Series)
        self.assertTrue(sa.equals(pd.Series([0.5, 0.9],
                                            index=[0, 4])))
        self.assertIsInstance(fa, pd.DataFrame)
        self.assertTrue(fa.equals(pd.DataFrame([[0.0, 0.5],
                                                [0.4, 0.9]],
                                               index=[5, 9],
                                               columns=["x", "y"])))


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
