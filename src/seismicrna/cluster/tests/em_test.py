import unittest as ut

import numpy as np
from scipy.stats import kstest

from seismicrna.cluster.em import (_expectation,
                                   _calc_jackpotting_g_stat,
                                   _calc_jackpotting_p_value)
from seismicrna.core.array import get_length

rng = np.random.default_rng(0)


class TestExpectation(ut.TestCase):

    def compare(self,
                p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                end5s: np.ndarray,
                end3s: np.ndarray,
                unmasked: np.ndarray,
                muts_per_pos: list[np.ndarray],
                min_mut_gap: int,
                expect_log_marginals: np.ndarray,
                expect_memberships: np.ndarray):
        result_log_marginals, result_memberships = _expectation(p_mut,
                                                                p_ends,
                                                                p_clust,
                                                                end5s,
                                                                end3s,
                                                                unmasked,
                                                                muts_per_pos,
                                                                min_mut_gap)
        self.assertEqual(result_log_marginals.shape, expect_log_marginals.shape)
        self.assertTrue(np.allclose(result_log_marginals, expect_log_marginals))
        self.assertEqual(expect_memberships.shape, result_memberships.shape)
        self.assertTrue(np.allclose(expect_memberships, result_memberships))

    def test_1pos(self):
        p_mut = np.array([[0.1]])
        p_ends = np.array([[1.]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0])
        end3s = np.array([0, 0])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1])]
        min_mut_gap = 0
        expect_log_marginals = np.log([0.9, 0.1])
        expect_memberships = np.ones((2, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

    def test_2pos_gap0(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3, 5]), np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.36, 0.04, 0.09, 0.01, 0.24, 0.06]
        )
        expect_memberships = np.ones((8, 1))
        min_mut_gap = 0
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

    def test_2pos_gap1(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3]), np.array([4, 6])]
        expect_log_marginals = np.log([0.18 / 0.99,
                                       0.02 / 0.99,
                                       0.36 / 0.99,
                                       0.04 / 0.99,
                                       0.09 / 0.99,
                                       0.24 / 0.99,
                                       0.06 / 0.99])
        expect_memberships = np.ones((7, 1))
        min_mut_gap = 1
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

    def test_2pos_masked0(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([1])
        muts_per_pos = [np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.20, 0.20, 0.40, 0.40, 0.10, 0.10, 0.24, 0.06]
        )
        expect_memberships = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         expect_log_marginals,
                         expect_memberships)

    def test_2pos_masked1(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1, 3, 5])]
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.45, 0.05, 0.45, 0.05, 0.30, 0.30]
        )
        expect_memberships = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         expect_log_marginals,
                         expect_memberships)

    def test_1pos_2clusters(self):
        p_mut = np.array([[0.1, 0.2]])
        p_ends = np.array([[1.]])
        p_clust = np.array([0.4, 0.6])
        end5s = np.array([0, 0])
        end3s = np.array([0, 0])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1])]
        min_mut_gap = 0
        expect_log_marginals = np.log([0.84, 0.16])
        expect_memberships = np.array([[3. / 7., 4. / 7.],
                                       [1. / 4., 3. / 4.]])
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

    def test_2pos_gap0_2clusters(self):
        p_mut = np.array([[0.1, 0.4],
                          [0.2, 0.3]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([0.4, 0.6])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3, 5]), np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.144, 0.056, 0.270, 0.100, 0.090, 0.040, 0.222, 0.078]
        )
        expect_memberships = np.array([[1. / 2., 1. / 2.],
                                       [1. / 7., 6. / 7.],
                                       [8. / 15., 7. / 15.],
                                       [4. / 25., 21. / 25.],
                                       [2. / 5., 3. / 5.],
                                       [1. / 10., 9. / 10.],
                                       [16. / 37., 21. / 37.],
                                       [4. / 13., 9. / 13.]])
        min_mut_gap = 0
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)


class TestCalcJackpottingGStat(ut.TestCase):

    def test_diff_lengths(self):
        num_obs = np.array([1])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               (r"Lengths differ between observed \(1\) "
                                r"and expected \(2\)"),
                               _calc_jackpotting_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_neg_num_obs(self):
        num_obs = np.array([-1, 2])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               r"All num_obs must be ≥ 0, but got \[-1\]",
                               _calc_jackpotting_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_zero_min_exp(self):
        num_obs = np.array([], dtype=int)
        num_exp = np.array([])
        min_exp = 0.
        self.assertRaisesRegex(ValueError,
                               r"min_exp must be > 0, but got 0.0",
                               _calc_jackpotting_g_stat,
                               num_obs,
                               np.log(num_exp),
                               min_exp)

    def test_each_exp_below_min_exp_all_exp_equal_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([2., 3., 1.])
        min_exp = 4.
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_equal_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.])
        min_exp = 4.
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_less_than_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 2.5])
        min_exp = 4.
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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
                               _calc_jackpotting_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_each_exp_at_least_min_exp_all_exp_equal_obs(self):
        min_obs = 4
        for n in range(5):
            num_obs = np.arange(min_obs, min_obs + n)
            num_exp = np.asarray(num_obs, dtype=float)
            g_stat, df = _calc_jackpotting_g_stat(num_obs,
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
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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
                               _calc_jackpotting_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_mixed_all_exp_equal_obs(self):
        num_obs = np.array([2, 4, 5, 8, 3])
        num_exp = np.asarray(num_obs, dtype=float)
        min_exp = 4.
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_mixed_sum_exp_equal_sum_obs(self):
        num_obs = np.array([0, 3, 7, 6, 6])
        num_exp = np.array([2., 4., 5., 8., 3.])
        min_exp = 4.
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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
        g_stat, df = _calc_jackpotting_g_stat(num_obs, np.log(num_exp), min_exp)
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


class TestJackpottingPValue(ut.TestCase):

    @staticmethod
    def simulate(pi2: float,
                 n_trials: int = 1024,
                 n_reads: int = 10000,
                 read_length: int = 24,
                 beta_a: float = 2.,
                 beta_b: float = 98.):
        p_values = list()
        for _ in range(n_trials):
            mus1 = rng.beta(beta_a, beta_b, size=read_length)
            mus2 = rng.beta(beta_a, beta_b, size=read_length)
            num_obs, log_exp = sim_obs_exp_mixture(n_reads, mus1, mus2, pi2)
            g_stat, df = _calc_jackpotting_g_stat(num_obs, log_exp)
            p_value = _calc_jackpotting_p_value(g_stat, df)
            p_values.append(p_value)
        return np.array(p_values)

    def test_without_jackpotting(self):
        """ Test that the P-value distribution is uniform for samples
        with no jackpotting. """
        pi2 = 0.0
        alpha = 0.01
        p_values = self.simulate(pi2)
        result = kstest(p_values, "uniform", nan_policy="raise")
        self.assertGreaterEqual(result.pvalue, alpha)

    def test_with_jackpotting(self):
        pi2 = 0.15
        alpha = 0.01
        p_values = self.simulate(pi2)
        result = kstest(p_values, "uniform", nan_policy="raise")
        self.assertLess(result.pvalue, alpha)


if __name__ == "__main__":
    ut.main()
