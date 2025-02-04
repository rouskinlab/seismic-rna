import unittest as ut
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from seismicrna.core.mu import (calc_arcsine_distance,
                                calc_sum_arcsine_distance,
                                calc_mean_arcsine_distance,
                                calc_coeff_determ,
                                calc_pearson,
                                calc_spearman,
                                compare_windows,
                                get_comp_func,
                                get_comp_name)
from seismicrna.core.seq import DNA, seq_pos_to_index

rng = np.random.default_rng()


class TestCalcArcsineDistance(ut.TestCase):

    def test_array0d(self):
        mus1 = np.array(0.1)
        mus2 = np.array(0.6)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(0.35926145215)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_array1d(self):
        mus1 = np.array([0.6])
        mus2 = np.array([0.1])
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array([0.35926145215])
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_array2d(self):
        mus1 = np.array([[0.2], [0.4]])
        mus2 = np.array([[0.3], [0.4]])
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array([[0.0738428842647], [0.0]])
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_series(self):
        index = pd.Index([4])
        mus1 = pd.Series([0.6], index)
        mus2 = pd.Series([0.1], index)
        result = calc_arcsine_distance(mus1, mus2)
        expect = pd.Series([0.35926145215], index)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))
        self.assertTrue(result.index.equals(index))

    def test_dataframe(self):
        index = pd.Index([4, 6])
        columns = pd.Index([8])
        mus1 = pd.DataFrame([[0.2], [0.4]], index, columns)
        mus2 = pd.DataFrame([[0.3], [0.4]], index, columns)
        result = calc_arcsine_distance(mus1, mus2)
        expect = pd.DataFrame([[0.0738428842647], [0.0]], index, columns)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))
        self.assertTrue(result.index.equals(index))
        self.assertTrue(result.columns.equals(columns))

    def test_both_zeros(self):
        mus1 = np.array(0.0)
        mus2 = np.array(0.0)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(0.0)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_both_ones(self):
        mus1 = np.array(1.0)
        mus2 = np.array(1.0)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(0.0)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_zero_half(self):
        mus1 = np.array(0.0)
        mus2 = np.array(0.5)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(0.5)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_one_half(self):
        mus1 = np.array(1.0)
        mus2 = np.array(0.5)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(0.5)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_zero_one(self):
        mus1 = np.array(0.0)
        mus2 = np.array(1.0)
        result = calc_arcsine_distance(mus1, mus2)
        expect = np.array(1.0)
        self.assertTupleEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcSumAbsDiffLogOdds(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_sum_arcsine_distance,
                               rng.random(()),
                               rng.random(()))

    def test_array1d(self):
        mus1 = np.array([0.4, 0.6, 0.3])
        mus2 = np.array([0.5, 0.1, 0.8])
        result = calc_sum_arcsine_distance(mus1, mus2)
        expect = np.sum([0.064094216849, 0.35926145215, 0.335822645134])
        self.assertIsInstance(result, float)
        self.assertTrue(np.isclose(result, expect))


class TestCalcMeanAbsFoldChangeOdds(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               calc_mean_arcsine_distance,
                               rng.random(()),
                               rng.random(()))

    def test_array1d(self):
        mus1 = np.array([0.4, 0.6, 0.3])
        mus2 = np.array([0.5, 0.1, 0.8])
        result = calc_mean_arcsine_distance(mus1, mus2)
        expect = np.mean([0.064094216849, 0.35926145215, 0.335822645134])
        self.assertIsInstance(result, float)
        self.assertTrue(np.isclose(result, expect))


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
        for key in ["MARCD", "marcd"]:
            self.assertIs(get_comp_func(key), calc_mean_arcsine_distance)
            self.assertEqual(get_comp_name(key),
                             "Mean Arcsine Distance")
        for key in ["PCC", "pcc"]:
            self.assertIs(get_comp_func(key), calc_pearson)
            self.assertEqual(get_comp_name(key),
                             "Pearson Correlation Coefficient")
        for key in ["SCC", "scc"]:
            self.assertIs(get_comp_func(key), calc_spearman)
            self.assertEqual(get_comp_name(key),
                             "Spearman Correlation Coefficient")
        for key in ["R2", "r2"]:
            self.assertIs(get_comp_func(key), calc_coeff_determ)
            self.assertEqual(get_comp_name(key),
                             "Coefficient of Determination")
        for get_comp in [get_comp_func, get_comp_name]:
            self.assertRaisesRegex(ValueError,
                                   "Invalid method of comparison: 'other'",
                                   get_comp,
                                   "other")


class TestCompareWindows(ut.TestCase):

    def test_contiguous(self):
        for seqlen in [0, 10, 20]:
            seq = DNA.random(seqlen)
            index = seq_pos_to_index(seq, np.arange(1, seqlen + 1), start=1)
            mus = pd.Series(rng.random(seqlen), index)
            for size in range(2, 6):
                for mc in range(2, size + 1):
                    nan5 = min(seqlen, max(0, mc - (1 + size // 2)))
                    nan3 = min(seqlen, max(0, mc - (1 + size) // 2))
                    nval = seqlen - (nan5 + nan3)
                    for method in ["MARCD", "PCC", "SCC", "R2"]:
                        result = compare_windows(mus, mus, method, size, mc)
                        self.assertIsInstance(result, pd.Series)
                        self.assertTrue(result.index.equals(index))
                        expect = np.hstack(
                            [np.full(nan5, np.nan),
                             np.full(nval, 0. if method == "MARCD" else 1.),
                             np.full(nan3, np.nan)]
                        )
                        self.assertTrue(np.allclose(result,
                                                    expect,
                                                    equal_nan=True))


if __name__ == "__main__":
    ut.main(verbosity=2)
