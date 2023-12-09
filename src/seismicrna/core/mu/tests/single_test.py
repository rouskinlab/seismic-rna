import unittest as ut
from logging import Filter, LogRecord

import numpy as np
import pandas as pd

from seismicrna.core.mu.single import (get_quantile,
                                       get_ranks,
                                       normalize,
                                       winsorize,
                                       logger as single_logger)

rng = np.random.default_rng()


class TestGetRanks(ut.TestCase):
    """ Test the function `get_ranks`. """

    def test_numpy_1d(self):
        """ Test with a 1D NumPy array. """
        mus = np.array([0.49, 0.17, 0.24, 0.90, 0.47, 0.50, 0.22, 0.14, 0.18])
        expect = np.array([6, 1, 4, 8, 5, 7, 3, 0, 2])
        ranks = get_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_numpy_2d(self):
        """ Test with a 2D NumPy array. """
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
        ranks = get_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_numpy_3d(self):
        """ Test with a 3D NumPy array. """
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
        ranks = get_ranks(mus)
        self.assertIsInstance(ranks, np.ndarray)
        self.assertTrue(np.array_equal(ranks, expect))

    def test_series(self):
        """ Test with a Series. """
        mus = pd.Series([0.38, 0.50, 0.51, 0.64, 0.75, 0.60, 0.18],
                        index=[1, 2, 4, 5, 6, 7, 9])
        expect = pd.Series([1, 2, 3, 5, 6, 4, 0], index=mus.index)
        ranks = get_ranks(mus)
        self.assertIsInstance(ranks, pd.Series)
        self.assertTrue(ranks.equals(expect))


class TestGetQuantile(ut.TestCase):
    """ Test the function `get_quantile`. """

    class NanFilter(Filter):
        """ Suppress warnings about NaN quantiles. """

        def filter(self, rec: LogRecord):
            """ Suppress warnings about NaN quantiles. """
            return not rec.msg.startswith("Got NaN quantile")

    nan_filter = NanFilter()

    def test_no_nan(self):
        """ Test with no NaN values. """
        for n in [5, 11, 19]:
            # Create a random order so that the NaN values are mixed in
            # with the finite values.
            order = np.arange(n)
            rng.shuffle(order)
            # Define quantiles of the array to check.
            quantiles = np.linspace(0., 1., n)
            for mu_max in np.linspace(0., 1., n):
                # Create an array of mutation rates from 0 to mu_max and
                # shuffle the values.
                mus = np.linspace(0., mu_max, n)[order]
                # Determine the value associated with each quantile.
                values = np.array([get_quantile(mus, quantile)
                                   for quantile in quantiles])
                # Since the values of mus were obtained via np.linspace,
                # the value for each quantile should be the quantile
                # times the maximum value of mus.
                self.assertTrue(np.allclose(values, quantiles * mu_max))

    def test_some_nan(self):
        """ Test with some (but not all) NaN values. """
        for n in [5, 11, 19]:
            for n_nan in [1, 3, 5]:
                # Create a random order so that the NaN values are mixed
                # in with the finite values.
                order = np.arange(n + n_nan)
                rng.shuffle(order)
                # Define quantiles of the array to check.
                quantiles = np.linspace(0., 1., n)
                # Test different maximum mutation rates.
                for mu_max in np.linspace(0., 1., n):
                    # Create an array of mutation rates from 0 to mu_max
                    # and shuffle the values.
                    mus = np.concatenate([np.linspace(0., mu_max, n),
                                          np.full(n_nan, np.nan)])[order]
                    # Determine the value associated with each quantile.
                    values = np.array([get_quantile(mus, quantile)
                                       for quantile in quantiles])
                    # Because the finite values of mus were obtained via
                    # np.linspace, the value for each quantile should be
                    # the quantile times the maximum value of mus.
                    self.assertTrue(np.allclose(values, quantiles * mu_max))

    def test_all_nan(self):
        """ Test that an all-NaN array returns NaN values. """
        for n in [5, 11, 19]:
            quantiles = np.linspace(0., 1., n)
            # Make mus an all-NaN array.
            mus = np.full(n, np.nan)
            # Temporarily suppress warnings about NaN values.
            single_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                single_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_empty(self):
        """ Test that an empty array always returns NaN values. """
        # Make mus an empty array.
        mus = np.array([], dtype=float)
        for n in [5, 11, 19]:
            quantiles = np.linspace(0., 1., n)
            # Temporarily suppress warnings about NaN values.
            single_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                single_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_invalid_quantiles(self):
        """ Test that invalid quantiles raise errors. """
        n = 11
        mus = rng.random(n)
        errmsg = r"Quantiles must be in the range [\[0, 1\]]"
        # Test that negative quantiles are invalid.
        for quantile in np.linspace(0., -1., n)[1:]:
            self.assertRaisesRegex(ValueError, errmsg,
                                   get_quantile, mus, quantile)
        # Test that quantiles greater than 1 are invalid.
        for quantile in np.linspace(1., 2., n)[1:]:
            self.assertRaisesRegex(ValueError, errmsg,
                                   get_quantile, mus, quantile)
        # Test that NaN is an invalid quantile.
        self.assertRaisesRegex(ValueError, errmsg,
                               get_quantile, mus, np.nan)


class TestNormalize(ut.TestCase):
    """ Test the function `mu.normalize`. """

    def test_normalize_p0(self):
        """ Do not normalize. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 0.0), mus))

    def test_normalize_p50(self):
        """ Normalize to the median. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 0.5), mus * 20.))

    def test_normalize_p100(self):
        """ Normalize to the maximum. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(normalize(mus, 1.0), mus * 10.))


class TestWinsorize(ut.TestCase):
    """ Test the function `mu.winsorize`. """

    def test_winsorize_p0(self):
        """ Do not winsorize. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 0.0), mus))

    def test_winsorize_p50(self):
        """ Winsorize to the median. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 0.5),
                                        np.where(mus < 0.05, mus * 20., 1.)))

    def test_winsorize_p100(self):
        """ Winsorize to the maximum. """
        for n in [5, 12, 19]:
            mus = np.linspace(0.0, 0.1, n)
            self.assertTrue(np.allclose(winsorize(mus, 1.0), mus * 10.))


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
