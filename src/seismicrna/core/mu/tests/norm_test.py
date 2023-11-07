import unittest as ut
from logging import Filter, LogRecord

import numpy as np

from ..norm import (get_quantile,
                    normalize,
                    winsorize,
                    logger as unbias_logger)
from ...rand import rng


class TestGetQuantile(ut.TestCase):
    """ Test the function `mu.get_quantile`. """

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
            unbias_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                unbias_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_empty(self):
        """ Test that an empty array always returns NaN values. """
        # Make mus an empty array.
        mus = np.array([], dtype=float)
        for n in [5, 11, 19]:
            quantiles = np.linspace(0., 1., n)
            # Temporarily suppress warnings about NaN values.
            unbias_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                unbias_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_invalid_quantiles(self):
        """ Test that invalid quantiles raise errors. """
        n = 11
        mus = rng.random(n)
        errmsg = "Quantiles must be in the range [\[0, 1\]]"
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
