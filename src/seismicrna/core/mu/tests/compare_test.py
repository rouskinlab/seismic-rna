import unittest as ut
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from seismicrna.core.mu.compare import (calc_coeff_determ,
                                        calc_pearson,
                                        calc_rmsd,
                                        calc_spearman)

rng = np.random.default_rng()


class TestCalcRMSD(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot count positions in 0-D arrays",
                               calc_rmsd,
                               np.array(0.),
                               np.array(0.))

    def test_array1d_allzero(self):
        for n in range(10):
            x = np.zeros(n, dtype=float)
            y = np.zeros(n, dtype=float)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                rmsd = calc_rmsd(x, y)
            self.assertIsInstance(rmsd, float)
            self.assertTrue(np.isnan(rmsd))

    def test_array1d_extremes(self):
        for fx in [0.0001, 0.01, 1.0]:
            for fy in [0.0001, 0.01, 1.0]:
                x = fx * np.array([0., 1.])
                y = fy * np.array([0., 1.])
                self.assertTrue(np.isclose(calc_rmsd(x, y), 0.))
                x = fx * np.array([0., 1.])
                y = fy * np.array([1., 0.])
                self.assertTrue(np.isclose(calc_rmsd(x, y), 1.))
                self.assertTrue(np.isclose(calc_rmsd(y, x), 1.))

    def test_array1d_examples(self):
        for fx in [0.0001, 0.01, 1.0]:
            for fy in [0.0001, 0.01, 1.0]:
                x = fx * np.array([1., 0., 0., 0.])
                y = fy * np.array([0., 0., 1., 0.])
                self.assertTrue(np.isclose(calc_rmsd(x, y), np.sqrt(1 / 2)))
                self.assertTrue(np.isclose(calc_rmsd(y, x), np.sqrt(1 / 2)))
                x = fx * np.array([0.4, 0.1, 0.8])
                y = fy * np.array([0.3, 0.2, 0.6])
                self.assertTrue(np.isclose(calc_rmsd(x, y), np.sqrt(1 / 72)))
                self.assertTrue(np.isclose(calc_rmsd(y, x), np.sqrt(1 / 72)))


class TestCalcPearson(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return pearsonr(x, y).statistic

    def test_array1d(self):
        # Vary number of rows.
        for nr in range(2, 10):
            x = rng.random(nr)
            y = rng.random(nr)
            s = calc_pearson(x, y)
            self.assertIsInstance(s, float)
            self.assertTrue(np.isclose(s, self.calc_true(x, y)))

    def test_array2d(self):
        # Vary number of columns.
        for nc in range(1, 3):
            # Vary number of rows.
            for nr in range(2, 10):
                x = rng.random((nr, nc))
                y = rng.random((nr, nc))
                s = calc_pearson(x, y)
                self.assertIsInstance(s, np.ndarray)
                self.assertEqual(s.shape, (nc,))
                # Compare the correlation for each column.
                for ic, sc in enumerate(s):
                    self.assertTrue(np.isclose(sc, self.calc_true(x[:, ic],
                                                                  y[:, ic])))

    def test_series(self):
        # Vary number of rows.
        for nr in range(2, 10):
            x = pd.Series(rng.random(nr))
            y = pd.Series(rng.random(nr))
            s = calc_pearson(x, y)
            self.assertIsInstance(s, float)
            self.assertTrue(np.isclose(s, self.calc_true(x, y)))

    def test_dataframe(self):
        # Vary number of columns.
        for nc in range(1, 3):
            # Vary number of rows.
            for nr in range(2, 10):
                x = pd.DataFrame(rng.random((nr, nc)))
                y = pd.DataFrame(rng.random((nr, nc)))
                s = calc_pearson(x, y)
                self.assertIsInstance(s, pd.Series)
                self.assertEqual(s.shape, (nc,))
                # Compare the correlation for each column.
                for ic, sc in enumerate(s):
                    self.assertTrue(np.isclose(
                        sc, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


class TestCalcCoeffDeterm(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return pearsonr(x, y).statistic ** 2

    def test_dataframe(self):
        # Vary number of columns.
        for nc in range(1, 3):
            # Vary number of rows.
            for nr in range(2, 10):
                x = pd.DataFrame(rng.random((nr, nc)))
                y = pd.DataFrame(rng.random((nr, nc)))
                s = calc_coeff_determ(x, y)
                self.assertIsInstance(s, pd.Series)
                self.assertEqual(s.shape, (nc,))
                # Compare the correlation for each column.
                for ic, sc in enumerate(s):
                    self.assertTrue(np.isclose(
                        sc, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


class TestCalcSpearman(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return spearmanr(x, y).statistic

    def test_array1d(self):
        # Vary number of rows.
        for nr in range(2, 10):
            x = rng.random(nr)
            y = rng.random(nr)
            s = calc_spearman(x, y)
            self.assertIsInstance(s, float)
            self.assertTrue(np.isclose(s, self.calc_true(x, y)))

    def test_array2d(self):
        # Vary number of columns.
        for nc in range(1, 3):
            # Vary number of rows.
            for nr in range(2, 10):
                x = rng.random((nr, nc))
                y = rng.random((nr, nc))
                s = calc_spearman(x, y)
                self.assertIsInstance(s, np.ndarray)
                self.assertEqual(s.shape, (nc,))
                # Compare the correlation for each column.
                for ic, sc in enumerate(s):
                    self.assertTrue(np.isclose(sc, self.calc_true(x[:, ic],
                                                                  y[:, ic])))

    def test_series(self):
        # Vary number of rows.
        for nr in range(2, 10):
            x = pd.Series(rng.random(nr))
            y = pd.Series(rng.random(nr))
            s = calc_spearman(x, y)
            self.assertIsInstance(s, float)
            self.assertTrue(np.isclose(s, self.calc_true(x, y)))

    def test_dataframe(self):
        # Vary number of columns.
        for nc in range(1, 3):
            # Vary number of rows.
            for nr in range(2, 10):
                x = pd.DataFrame(rng.random((nr, nc)))
                y = pd.DataFrame(rng.random((nr, nc)))
                s = calc_spearman(x, y)
                self.assertIsInstance(s, pd.Series)
                self.assertEqual(s.shape, (nc,))
                # Compare the correlation for each column.
                for ic, sc in enumerate(s):
                    self.assertTrue(np.isclose(
                        sc, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


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
