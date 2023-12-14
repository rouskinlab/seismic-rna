import unittest as ut
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from seismicrna.core.mu import (calc_coeff_determ,
                                calc_pearson,
                                calc_rmsd,
                                calc_spearman,
                                compare_windows,
                                get_comp_abbr,
                                get_comp_func,
                                get_comp_name)
from seismicrna.core.seq import DNA, seq_pos_to_index

rng = np.random.default_rng()


class TestCalcRMSD(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_rmsd,
                               rng.random(()),
                               rng.random(()))

    def test_array1d_allzero(self):
        for n in range(5):
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
                self.assertTrue(np.isclose(calc_rmsd(x, y), 2. ** -0.5))
                self.assertTrue(np.isclose(calc_rmsd(y, x), 2. ** -0.5))
                x = fx * np.array([0.4, 0.1, 0.8])
                y = fy * np.array([0.3, 0.2, 0.6])
                self.assertTrue(np.isclose(calc_rmsd(x, y), 72. ** -0.5))
                self.assertTrue(np.isclose(calc_rmsd(y, x), 72. ** -0.5))
                x = fx * np.array([np.nan, 0.4, 0.1, 0.3, 0.8])
                y = fy * np.array([0.5, 0.3, 0.2, np.nan, 0.6])
                self.assertTrue(np.isclose(calc_rmsd(x, y), 72. ** -0.5))
                self.assertTrue(np.isclose(calc_rmsd(y, x), 72. ** -0.5))

    def test_array2d(self):
        x = np.array([[0.8, 0.0],
                      [0.9, np.nan],
                      [0.4, 0.4],
                      [0.5, 0.6],
                      [0.1, 0.0]])
        y = np.array([[0.6, 0.9],
                      [1.0, 0.2],
                      [0.3, 0.0],
                      [np.nan, 0.3],
                      [0.2, 0.0]])
        rmsd = calc_rmsd(x, y)
        self.assertIsInstance(rmsd, np.ndarray)
        self.assertEqual(rmsd.shape, (2,))
        self.assertTrue(np.allclose(rmsd, [72. ** -0.5,
                                           (2. / 3.) ** 0.5]))

    def test_series(self):
        x = pd.Series([np.nan, 0.4, 0.1, 0.3, 0.8])
        y = pd.Series([0.5, 0.3, 0.2, np.nan, 0.6])
        self.assertTrue(np.isclose(calc_rmsd(x, y), 72. ** -0.5))
        self.assertTrue(np.isclose(calc_rmsd(y, x), 72. ** -0.5))

    def test_dataframe(self):
        index = pd.Index([2, 4, 5, 7, 9])
        x = pd.DataFrame([[0.8, 0.0],
                          [0.9, np.nan],
                          [0.4, 0.4],
                          [0.5, 0.6],
                          [0.1, 0.0]],
                         index=index,
                         columns=["i", "j"])
        y = pd.DataFrame([[0.9, 0.6],
                          [0.2, 1.0],
                          [0.0, 0.3],
                          [0.3, np.nan],
                          [0.0, 0.2]],
                         index=index,
                         columns=["j", "i"])
        rmsd = calc_rmsd(x, y)
        self.assertIsInstance(rmsd, pd.Series)
        self.assertEqual(rmsd.shape, (2,))
        self.assertTrue(np.allclose(rmsd, [72. ** -0.5,
                                           (2. / 3.) ** 0.5]))
        self.assertTrue(rmsd.index.equals(x.columns))


class TestCalcPearson(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return pearsonr(x, y).statistic

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_pearson,
                               rng.random(()),
                               rng.random(()))

    def test_array1d(self):
        for nrow in range(10):
            x = rng.random(nrow)
            y = rng.random(nrow)
            if nrow >= 2:
                r = calc_pearson(x, y)
                self.assertIsInstance(r, float)
                self.assertTrue(np.isclose(r, self.calc_true(x, y)))
            else:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",
                                            "Mean of empty slice",
                                            RuntimeWarning)
                    with np.errstate(invalid="ignore"):
                        r = calc_pearson(x, y)
                        self.assertIsInstance(r, float)
                        self.assertTrue(np.isnan(r))

    def test_array2d(self):
        for ncol in range(1, 3):
            for nrow in range(2, 10):
                x = rng.random((nrow, ncol))
                y = rng.random((nrow, ncol))
                r = calc_pearson(x, y)
                self.assertIsInstance(r, np.ndarray)
                self.assertEqual(r.shape, (ncol,))
                # Compare the correlation for each column.
                for ic, rc in enumerate(r):
                    self.assertTrue(np.isclose(rc, self.calc_true(x[:, ic],
                                                                  y[:, ic])))

    def test_series(self):
        for nrow in range(2, 10):
            x = pd.Series(rng.random(nrow))
            y = pd.Series(rng.random(nrow))
            r = calc_pearson(x, y)
            self.assertIsInstance(r, float)
            self.assertTrue(np.isclose(r, self.calc_true(x, y)))

    def test_dataframe(self):
        for ncol in range(1, 3):
            for nrow in range(2, 10):
                x = pd.DataFrame(rng.random((nrow, ncol)))
                y = pd.DataFrame(rng.random((nrow, ncol)))
                r = calc_pearson(x, y)
                self.assertIsInstance(r, pd.Series)
                self.assertEqual(r.shape, (ncol,))
                # Compare the correlation for each column.
                for ic, rc in enumerate(r):
                    self.assertTrue(np.isclose(
                        rc, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


class TestCalcCoeffDeterm(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return pearsonr(x, y).statistic ** 2

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_coeff_determ,
                               rng.random(()),
                               rng.random(()))

    def test_array1d(self):
        for nrow in range(10):
            x = rng.random(nrow)
            y = rng.random(nrow)
            if nrow >= 2:
                r2 = calc_coeff_determ(x, y)
                self.assertIsInstance(r2, float)
                self.assertTrue(np.isclose(r2, self.calc_true(x, y)))
            else:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",
                                            "Mean of empty slice",
                                            RuntimeWarning)
                    with np.errstate(invalid="ignore"):
                        r2 = calc_coeff_determ(x, y)
                        self.assertIsInstance(r2, float)
                        self.assertTrue(np.isnan(r2))

    def test_dataframe(self):
        for ncol in range(1, 3):
            for nrow in range(2, 10):
                x = pd.DataFrame(rng.random((nrow, ncol)))
                y = pd.DataFrame(rng.random((nrow, ncol)))
                r2 = calc_coeff_determ(x, y)
                self.assertIsInstance(r2, pd.Series)
                self.assertEqual(r2.shape, (ncol,))
                # Compare the correlation for each column.
                for ic, r2c in enumerate(r2):
                    self.assertTrue(np.isclose(
                        r2c, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


class TestCalcSpearman(ut.TestCase):

    @classmethod
    def calc_true(cls, x, y):
        """ Calculate the "true" coefficient using a trusted method. """
        return spearmanr(x, y, nan_policy="omit").statistic

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_spearman,
                               rng.random(()),
                               rng.random(()))

    def test_array1d(self):
        for nrow in range(10):
            x = rng.random(nrow)
            y = rng.random(nrow)
            if nrow >= 2:
                rho = calc_spearman(x, y)
                self.assertIsInstance(rho, float)
                self.assertTrue(np.isclose(rho, self.calc_true(x, y)))
            else:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",
                                            "Mean of empty slice",
                                            RuntimeWarning)
                    with np.errstate(invalid="ignore"):
                        rho = calc_spearman(x, y)
                        self.assertIsInstance(rho, float)
                        self.assertTrue(np.isnan(rho))

    def test_array1d_nan(self):
        for nrow in range(10):
            for num_nan in range(nrow + 1):
                x = rng.random(nrow)
                x[:num_nan] = np.nan
                self.assertEqual(np.count_nonzero(np.isnan(x)), num_nan)
                y = rng.random(nrow)
                if nrow - num_nan >= 2:
                    rho = calc_spearman(x, y)
                    self.assertIsInstance(rho, float)
                    errstate = "warn" if nrow - num_nan >= 3 else "ignore"
                    with np.errstate(invalid=errstate):
                        self.assertTrue(np.isclose(rho, self.calc_true(x, y)))
                else:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore",
                                                "Mean of empty slice",
                                                RuntimeWarning)
                        with np.errstate(invalid="ignore"):
                            rho = calc_spearman(x, y)
                            self.assertIsInstance(rho, float)
                            self.assertTrue(np.isnan(rho))

    def test_array2d(self):
        for ncol in range(1, 3):
            for nrow in range(2, 10):
                x = rng.random((nrow, ncol))
                y = rng.random((nrow, ncol))
                rho = calc_spearman(x, y)
                self.assertIsInstance(rho, np.ndarray)
                self.assertEqual(rho.shape, (ncol,))
                # Compare the correlation for each column.
                for ic, rhoc in enumerate(rho):
                    self.assertTrue(np.isclose(rhoc, self.calc_true(x[:, ic],
                                                                    y[:, ic])))

    def test_series(self):
        for nrow in range(2, 10):
            x = pd.Series(rng.random(nrow))
            y = pd.Series(rng.random(nrow))
            rho = calc_spearman(x, y)
            self.assertIsInstance(rho, float)
            self.assertTrue(np.isclose(rho, self.calc_true(x, y)))

    def test_dataframe(self):
        for ncol in range(1, 3):
            for nrow in range(2, 10):
                x = pd.DataFrame(rng.random((nrow, ncol)))
                y = pd.DataFrame(rng.random((nrow, ncol)))
                rho = calc_spearman(x, y)
                self.assertIsInstance(rho, pd.Series)
                self.assertEqual(rho.shape, (ncol,))
                # Compare the correlation for each column.
                for ic, rhoc in enumerate(rho):
                    self.assertTrue(np.isclose(
                        rhoc, self.calc_true(x.iloc[:, ic], y.iloc[:, ic])
                    ))


class TestGetComp(ut.TestCase):

    def test_comps(self):
        self.assertIs(get_comp_func("rmsd"), calc_rmsd)
        self.assertEqual(get_comp_abbr("rmsd"), "RMSD")
        self.assertEqual(get_comp_name("rmsd"),
                         "Root-Mean-Square Deviation")
        self.assertIs(get_comp_func("r"), calc_pearson)
        self.assertEqual(get_comp_abbr("r"), "PCC")
        self.assertEqual(get_comp_name("r"),
                         "Pearson Correlation Coefficient")
        self.assertIs(get_comp_func("r2"), calc_coeff_determ)
        self.assertEqual(get_comp_abbr("r2"), "R-Squared")
        self.assertEqual(get_comp_name("r2"),
                         "Coefficient of Determination")
        self.assertIs(get_comp_func("rho"), calc_spearman)
        self.assertEqual(get_comp_abbr("rho"), "SCC")
        self.assertEqual(get_comp_name("rho"),
                         "Spearman Correlation Coefficient")
        for get_comp in [get_comp_func, get_comp_name, get_comp_abbr]:
            self.assertRaisesRegex(ValueError,
                                   "Invalid method of comparison: 'other'",
                                   get_comp,
                                   "other")


class TestCompareWindows(ut.TestCase):

    def test_contiguous(self):
        for seqlen in [0, 10, 20]:
            seq = DNA.random(seqlen)
            index = seq_pos_to_index(seq, np.arange(1, seqlen + 1), start=1)
            mus1 = pd.Series(rng.random(seqlen), index)
            mus2 = pd.Series(rng.random(seqlen), index)
            for winlen in range(2, 6):
                for mucomp in [calc_rmsd, calc_pearson, calc_spearman]:
                    values = compare_windows(mus1, mus2, mucomp, winlen)
                    self.assertIsInstance(values, pd.Series)
                    self.assertTrue(values.index.equals(index))


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
