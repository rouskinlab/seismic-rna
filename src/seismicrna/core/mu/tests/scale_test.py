import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.mu import (calc_quantile,
                                calc_ranks,
                                normalize,
                                winsorize)

rng = np.random.default_rng()


class TestCalcQuantile(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_quantile,
                               np.array(0.5),
                               0.5)

    def test_array1d(self):
        for n in range(1, 5):
            # Create a random order.
            order = np.arange(n)
            rng.shuffle(order)
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., n)
            # Test multiple maximum mutation rates.
            for mu_max in np.linspace(0., 1., n):
                # Create an array of mutation rates from 0 to mu_max and
                # shuffle the values.
                mus = np.linspace(0., mu_max, n)[order]
                # Determine the value associated with each quantile.
                for quantile in quantiles:
                    value = calc_quantile(mus, quantile)
                    # Since mus was generated via np.linspace, the value
                    # of each quantile should be the quantile times the
                    # maximum value of mus.
                    self.assertIsInstance(value, np.ndarray)
                    self.assertEqual(value.shape, ())
                    self.assertTrue(np.isclose(value, quantile * mu_max))

    def test_array1d_some_nan(self):
        for n in range(1, 5):
            for n_nan in range(1, n - 1):
                # Create a random order so that the NaN values are mixed
                # in with the finite values.
                order = np.arange(n + n_nan)
                rng.shuffle(order)
                # Define quantiles of the array to check.
                quantiles = np.linspace(0., 1., n)
                # Test multiple maximum mutation rates.
                for mu_max in np.linspace(0., 1., n):
                    # Create an array of mutation rates from 0 to mu_max
                    # and shuffle the values.
                    mus = np.concatenate([np.linspace(0., mu_max, n),
                                          np.full(n_nan, np.nan)])[order]
                    # Determine the value associated with each quantile.
                    for quantile in quantiles:
                        value = calc_quantile(mus, quantile)
                        # Since mus was generated via np.linspace, the
                        # value of each quantile should be the quantile
                        # times the maximum value of mus.
                        self.assertIsInstance(value, np.ndarray)
                        self.assertEqual(value.shape, ())
                        self.assertTrue(np.isclose(value, quantile * mu_max))

    def test_array1d_all_nan(self):
        for n in range(1, 5):
            quantiles = np.linspace(0., 1., n)
            # Make mus an all-NaN array.
            mus = np.full(n, np.nan)
            # Determine the value associated with each quantile.
            for quantile in quantiles:
                value = calc_quantile(mus, quantile)
                # The value of every quantile should be NaN.
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, ())
                self.assertTrue(np.isnan(value))

    def test_array1d_empty(self):
        # Make mus an empty array.
        mus = np.array([], dtype=float)
        for n in range(1, 5):
            quantiles = np.linspace(0., 1., n)
            # Determine the value associated with each quantile.
            for quantile in quantiles:
                value = calc_quantile(mus, quantile)
                # The value of every quantile should be NaN.
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, ())
                self.assertTrue(np.isnan(value))

    def test_series(self):
        for n in range(1, 5):
            # Create a random order.
            order = np.arange(n)
            rng.shuffle(order)
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., n)
            # Test multiple maximum mutation rates.
            for mu_max in np.linspace(0., 1., n):
                # Create a Series of mutation rates from 0 to mu_max and
                # shuffle the values.
                mus = pd.Series(np.linspace(0., mu_max, n)[order])
                # Determine the value associated with each quantile.
                for quantile in quantiles:
                    value = calc_quantile(mus, quantile)
                    # Since mus was generated via np.linspace, the value
                    # of each quantile should be the quantile times the
                    # maximum value of mus.
                    self.assertIsInstance(value, np.ndarray)
                    self.assertEqual(value.shape, ())
                    self.assertTrue(np.isclose(value, quantile * mu_max))

    def test_array2d(self):
        for nrow in range(1, 5):
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., nrow)
            for ncol in range(4):
                # Create a random order for each column.
                orders = list()
                for _ in range(ncol):
                    order = np.arange(nrow)
                    rng.shuffle(order)
                    orders.append(order)
                # Test multiple maximum mutation rates.
                for mu_max in np.linspace(0., 1., nrow):
                    # Create an array of mutation rates from 0 to mu_max and
                    # shuffle the values.
                    mus = np.empty((nrow, ncol))
                    for col, order in enumerate(orders):
                        mus[:, col] = np.linspace(0., mu_max, nrow)[order]
                    # Determine the value associated with each quantile.
                    for quantile in quantiles:
                        value = calc_quantile(mus, quantile)
                        # Since mus was generated via np.linspace, the value
                        # of each quantile should be the quantile times the
                        # maximum value of mus.
                        self.assertIsInstance(value, np.ndarray)
                        self.assertEqual(value.shape, (ncol,))
                        self.assertTrue(np.allclose(value, quantile * mu_max))

    def test_array2d_one_row_nan(self):
        for nrow in range(1, 5):
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., nrow)
            for ncol in range(4):
                # Test multiple maximum mutation rates.
                for mu_max in np.linspace(0., 1., nrow):
                    for nanrow in range(nrow + 1):
                        # Create an array of mutation rates.
                        mus = np.array(np.broadcast_to(
                            np.linspace(0., mu_max, nrow)[:, np.newaxis],
                            (nrow, ncol)
                        ))
                        # Add one all-NaN row.
                        mus = np.vstack([mus[:nanrow],
                                         np.full((1, ncol), np.nan),
                                         mus[nanrow:]])
                        # The mutation rates should have one extra row
                        # containing only NaN values.
                        self.assertEqual(mus.shape, (nrow + 1, ncol))
                        # Determine the value for each quantile.
                        for quantile in quantiles:
                            value = calc_quantile(mus, quantile)
                            # Since mus was generated via np.linspace,
                            # the value of each quantile should be the
                            # quantile times the maximum value of mus.
                            self.assertIsInstance(value, np.ndarray)
                            self.assertEqual(value.shape, (ncol,))
                            self.assertTrue(np.allclose(value, quantile * mu_max))

    def test_array2d_one_col_nan(self):
        for nrow in range(1, 5):
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., nrow)
            for ncol in range(4):
                # Create a random order for each column.
                orders = list()
                for _ in range(ncol):
                    order = np.arange(nrow)
                    rng.shuffle(order)
                    orders.append(order)
                # Test multiple maximum mutation rates.
                for mu_max in np.linspace(0., 1., nrow):
                    for nancol in range(ncol):
                        # Create an array of mutation rates from 0 to
                        # mu_max and shuffle the values, setting one
                        # column to all NaN values.
                        mus = np.empty((nrow, ncol))
                        for col, order in enumerate(orders):
                            mus[:, col] = (np.linspace(0., mu_max, nrow)[order]
                                           if col != nancol
                                           else np.nan)
                        # Determine the value for each quantile.
                        for quantile in quantiles:
                            value = calc_quantile(mus, quantile)
                            # Since every position (row) has one NaN
                            # value, all positions should be dropped
                            # and thus the result should be all-NaN.
                            self.assertIsInstance(value, np.ndarray)
                            self.assertEqual(value.shape, (ncol,))
                            self.assertTrue(np.all(np.isnan(value)))

    def test_dataframe(self):
        for nrow in range(1, 5):
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., nrow)
            for ncol in range(4):
                # Create a random order for each column.
                orders = list()
                for _ in range(ncol):
                    order = np.arange(nrow)
                    rng.shuffle(order)
                    orders.append(order)
                # Test multiple maximum mutation rates.
                for mu_max in np.linspace(0., 1., nrow):
                    # Create an array of mutation rates from 0 to mu_max and
                    # shuffle the values.
                    mus = np.empty((nrow, ncol))
                    for col, order in enumerate(orders):
                        mus[:, col] = np.linspace(0., mu_max, nrow)[order]
                    mus = pd.DataFrame(mus,
                                       rng.integers(10, size=nrow),
                                       rng.integers(10, size=ncol))
                    # Determine the value associated with each quantile.
                    for quantile in quantiles:
                        value = calc_quantile(mus, quantile)
                        # Since mus was generated via np.linspace, the value
                        # of each quantile should be the quantile times the
                        # maximum value of mus.
                        self.assertIsInstance(value, pd.Series)
                        self.assertEqual(value.shape, (ncol,))
                        self.assertTrue(np.allclose(value, quantile * mu_max))
                        self.assertTrue(value.index.equals(mus.columns))

    def test_invalid_quantiles(self):
        n = 11
        mus = rng.random(n)
        valerr = r"Quantiles must be in the range [\[0, 1\]]"
        # Test that negative quantiles are invalid.
        for quantile in np.linspace(0., -1., n)[1:]:
            self.assertRaisesRegex(ValueError,
                                   valerr,
                                   calc_quantile,
                                   mus,
                                   quantile)
        # Test that quantiles greater than 1 are invalid.
        for quantile in np.linspace(1., 2., n)[1:]:
            self.assertRaisesRegex(ValueError,
                                   valerr,
                                   calc_quantile,
                                   mus,
                                   quantile)
        # Test that NaN is an invalid quantile.
        self.assertRaisesRegex(ValueError,
                               valerr,
                               calc_quantile,
                               mus,
                               np.nan)
        # Test that non-floats are invalid.
        for quantile in [0, 1, np.array([0.5]), np.linspace(0., 1., n)]:
            self.assertRaisesRegex(TypeError,
                                   "Expected quantile to be float, but got",
                                   calc_quantile,
                                   mus,
                                   quantile)


class TestNormalize(ut.TestCase):

    def test_normalize_p0(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 0.0), mus))

    def test_normalize_p50(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 0.5), mus * 20.))

    def test_normalize_p100(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 1.0), mus * 10.))


class TestWinsorize(ut.TestCase):

    def test_winsorize_p0(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 0.0), mus))

    def test_winsorize_p50(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 0.5),
                                        np.where(mus < 0.05, mus * 20., 1.)))

    def test_winsorize_p100(self):
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 1.0), mus * 10.))


class TestCalcRanks(ut.TestCase):

    def test_array0d(self):
        mus = np.array(0.52)
        self.assertRaisesRegex(ValueError,
                               "Cannot rank mutation rates in a 0-D array",
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
        self.assertRaisesRegex(ValueError,
                               "Cannot rank mutation rates with NaN values",
                               calc_ranks,
                               mus)


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
