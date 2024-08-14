import unittest as ut

import numpy as np

from seismicrna.cluster.jackpot import calc_jackpot_g_stat
from seismicrna.core.array import get_length

rng = np.random.default_rng(0)


class TestCalcJackpottingGStat(ut.TestCase):

    def test_diff_lengths(self):
        num_obs = np.array([1])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               (r"Lengths differ between observed \(1\) "
                                r"and expected \(2\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_neg_num_obs(self):
        num_obs = np.array([-1, 2])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               r"All num_obs must be ≥ 0, but got \[-1\]",
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_neg_min_exp(self):
        num_obs = np.array([], dtype=int)
        num_exp = np.array([])
        min_exp = -0.1
        self.assertRaisesRegex(ValueError,
                               r"min_exp must be ≥ 0, but got -0.1",
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp),
                               min_exp)

    def test_each_exp_below_min_exp_all_exp_equal_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([2., 3., 1.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_equal_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_less_than_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 2.5])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 1)

    def test_each_exp_below_min_exp_sum_exp_greater_than_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.5])
        self.assertRaisesRegex(ValueError,
                               (r"Total observed reads \(6\) is less than "
                                r"total expected reads \(6.5\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_each_exp_at_least_min_exp_all_exp_equal_obs(self):
        min_obs = 4
        for n in range(5):
            num_obs = np.arange(min_obs, min_obs + n)
            num_exp = np.asarray(num_obs, dtype=float)
            g_stat, df = calc_jackpot_g_stat(num_obs,
                                             np.log(num_exp),
                                             float(min_obs))
            self.assertIsInstance(g_stat, float)
            self.assertEqual(g_stat, 0.)
            self.assertIsInstance(df, int)
            self.assertEqual(df, max(n - 1, 0))

    def test_each_exp_at_least_min_exp_sum_exp_equal_sum_obs(self):
        num_obs = np.array([5, 4, 10])
        num_exp = np.array([4., 5., 10.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   sum([10. * np.log(1.25),
                                        8. * np.log(0.8)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 2)

    def test_each_exp_at_least_min_exp_sum_exp_less_than_sum_obs(self):
        num_obs = np.array([5, 4, 11])
        num_exp = np.array([4., 5., 10.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([5 * np.log(5 / 4),
                                             4 * np.log(4 / 5),
                                             11 * np.log(11 / 10)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_each_exp_at_least_min_exp_sum_exp_greater_than_sum_obs(self):
        num_obs = np.array([5, 4, 10])
        num_exp = np.array([4., 5., 11.])
        self.assertRaisesRegex(ValueError,
                               (r"Total observed reads \(19\) is less than "
                                r"total expected reads \(20.0\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_mixed_all_exp_equal_obs(self):
        num_obs = np.array([2, 4, 5, 8, 3])
        num_exp = np.asarray(num_obs, dtype=float)
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_mixed_sum_exp_equal_sum_obs(self):
        num_obs = np.array([0, 3, 7, 6, 6])
        num_exp = np.array([2., 4., 5., 8., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             7 * np.log(7 / 5),
                                             6 * np.log(6 / 8),
                                             6 * np.log(6 / 5)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_mixed_sum_exp_under_min_exp_less_than_sum_obs(self):
        num_obs = np.array([0, 3, 7, 6, 6])
        num_exp = np.array([1., 4., 5., 8., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             7 * np.log(7 / 5),
                                             6 * np.log(6 / 8),
                                             6 * np.log(6 / 5)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 4)

    def test_mixed_sum_exp_over_min_exp_less_than_sum_obs(self):
        num_obs = np.array([0, 3, 5, 6, 8])
        num_exp = np.array([2., 4., 5., 7., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             5 * np.log(5 / 5),
                                             6 * np.log(6 / 7),
                                             8 * np.log(8 / 6)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 4)


def sim_reads(n: int, mus: np.ndarray):
    return rng.random((n, get_length(mus))) < mus


def count_uniq_reads(reads: np.ndarray):
    uniq_reads, num_obs = np.unique(reads, axis=0, return_counts=True)
    return uniq_reads, num_obs


def calc_log_exp(n: int, mus: np.ndarray, uniq_reads: np.ndarray):
    return np.log(n) + np.sum(np.where(uniq_reads,
                                       np.log(mus),
                                       np.log(1. - mus)),
                              axis=1)


def sim_obs_exp_mixture(n: int,
                        mus1: np.ndarray,
                        mus2: np.ndarray,
                        pi2: float):
    """ Simulate observed and expected counts from a mixture of two
    sets of mutation rates. """
    n2 = round(n * pi2)
    n1 = n - n2
    reads = np.vstack([sim_reads(n1, mus1),
                       sim_reads(n2, mus2)])
    uniq_reads, num_obs = count_uniq_reads(reads)
    # Calculate the expected counts using the average of mus1 and mus2.
    # The discrepancy between this average and the actual underlying
    # structure of the simulated reads is what produces the appearance
    # of jackpotting.
    mus = (n1 / n) * mus1 + (n2 / n) * mus2
    log_exp = calc_log_exp(n, mus, uniq_reads)
    return num_obs, log_exp


if __name__ == "__main__":
    ut.main()
