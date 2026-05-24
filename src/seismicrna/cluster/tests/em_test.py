import unittest as ut

import numpy as np

from seismicrna.cluster.em import _calc_bic, _calc_log_like


class TestCalcBic(ut.TestCase):
    def test_formula(self):
        # BIC = n_params * log(n_data) - 2 * log_like
        result = _calc_bic(n_params=3, n_data=100, log_like=-50.0)
        expected = 3 * np.log(100) - 2 * (-50.0)
        self.assertAlmostEqual(result, expected)

    def test_zero_params(self):
        result = _calc_bic(n_params=0, n_data=100, log_like=-10.0)
        self.assertAlmostEqual(result, -2 * (-10.0))

    def test_large_n_data(self):
        result = _calc_bic(n_params=5, n_data=10000, log_like=-200.0)
        expected = 5 * np.log(10000) - 2 * (-200.0)
        self.assertAlmostEqual(result, expected)

    def test_zero_log_like(self):
        # log_like = 0 is allowed (boundary)
        result = _calc_bic(n_params=2, n_data=50, log_like=0.0)
        expected = 2 * np.log(50)
        self.assertAlmostEqual(result, expected)

    def test_positive_log_like_raises(self):
        with self.assertRaises(ValueError):
            _calc_bic(n_params=2, n_data=50, log_like=0.1)

    def test_negative_n_params_raises(self):
        with self.assertRaises(ValueError):
            _calc_bic(n_params=-1, n_data=50, log_like=-10.0)

    def test_negative_n_data_raises(self):
        with self.assertRaises(ValueError):
            _calc_bic(n_params=2, n_data=-1, log_like=-10.0)

    def test_zero_n_data_returns_nan(self):
        # n_data=0 → n_params * log(0) = -inf; verify it doesn't raise
        result = _calc_bic(n_params=2, n_data=0, log_like=0.0)
        self.assertTrue(np.isneginf(result))

    def test_returns_float(self):
        result = _calc_bic(n_params=1, n_data=10, log_like=-5.0)
        self.assertIsInstance(result, float)


class TestCalcLogLike(ut.TestCase):
    def test_single_read(self):
        logp = np.array([-2.0])
        counts = np.array([1])
        self.assertAlmostEqual(_calc_log_like(logp, counts), -2.0)

    def test_single_uniq_multiple_counts(self):
        logp = np.array([-1.5])
        counts = np.array([4])
        self.assertAlmostEqual(_calc_log_like(logp, counts), -6.0)

    def test_multiple_uniqs(self):
        logp = np.array([-1.0, -2.0, -3.0])
        counts = np.array([2, 3, 1])
        # vdot = 2*(-1) + 3*(-2) + 1*(-3) = -2 - 6 - 3 = -11
        self.assertAlmostEqual(_calc_log_like(logp, counts), -11.0)

    def test_zero_counts(self):
        logp = np.array([-1.0, -2.0])
        counts = np.array([0, 0])
        self.assertAlmostEqual(_calc_log_like(logp, counts), 0.0)

    def test_returns_float(self):
        result = _calc_log_like(np.array([-1.0]), np.array([1]))
        self.assertIsInstance(result, float)


if __name__ == "__main__":
    ut.main(verbosity=2)
