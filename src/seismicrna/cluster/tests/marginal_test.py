import unittest as ut
from collections import Counter

import numpy as np

from seismicrna.cluster.marginal import calc_marginal, calc_marginal_resps
from seismicrna.core.tests import unbias_test
from seismicrna.core.tests.unbias_test import (label_no_close_muts,
                                               merge_mutations_right_to_left,
                                               simulate_reads)


class TestMarginalResps(ut.TestCase):

    def compare(self,
                p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                end5s: np.ndarray,
                end3s: np.ndarray,
                unmasked: np.ndarray,
                muts_per_pos: list[np.ndarray],
                min_mut_gap: int,
                mut_collisions: str,
                expect_log_marginals: np.ndarray,
                expect_resps: np.ndarray):
        result_log_marginals, result_resps = calc_marginal_resps(p_mut,
                                                                 p_ends,
                                                                 p_clust,
                                                                 end5s,
                                                                 end3s,
                                                                 unmasked,
                                                                 muts_per_pos,
                                                                 min_mut_gap,
                                                                 mut_collisions)
        self.assertEqual(result_log_marginals.shape, expect_log_marginals.shape)
        self.assertTrue(np.allclose(result_log_marginals, expect_log_marginals))
        self.assertEqual(expect_resps.shape, result_resps.shape)
        self.assertTrue(np.allclose(expect_resps, result_resps))

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
        expect_resps = np.ones((2, 1))
        for mut_collisions in ["drop", "merge"]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         mut_collisions,
                         expect_log_marginals,
                         expect_resps)

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
        expect_resps = np.ones((8, 1))
        min_mut_gap = 0
        for mut_collisions in ["drop", "merge"]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         mut_collisions,
                         expect_log_marginals,
                         expect_resps)

    def test_2pos_gap1_drop(self):
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
        expect_resps = np.ones((7, 1))
        min_mut_gap = 1
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     "drop",
                     expect_log_marginals,
                     expect_resps)

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
        expect_resps = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            for mut_collisions in ["drop", "merge"]:
                self.compare(p_mut,
                             p_ends,
                             p_clust,
                             end5s,
                             end3s,
                             unmasked,
                             muts_per_pos,
                             min_mut_gap,
                             mut_collisions,
                             expect_log_marginals,
                             expect_resps)

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
        expect_resps = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            for mut_collisions in ["drop", "merge"]:
                self.compare(p_mut,
                             p_ends,
                             p_clust,
                             end5s,
                             end3s,
                             unmasked,
                             muts_per_pos,
                             min_mut_gap,
                             mut_collisions,
                             expect_log_marginals,
                             expect_resps)

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
        expect_resps = np.array([[3. / 7., 4. / 7.],
                                 [1. / 4., 3. / 4.]])
        for mut_collisions in ["drop", "merge"]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         mut_collisions,
                         expect_log_marginals,
                         expect_resps)

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
        expect_resps = np.array([[1. / 2., 1. / 2.],
                                 [1. / 7., 6. / 7.],
                                 [8. / 15., 7. / 15.],
                                 [4. / 25., 21. / 25.],
                                 [2. / 5., 3. / 5.],
                                 [1. / 10., 9. / 10.],
                                 [16. / 37., 21. / 37.],
                                 [4. / 13., 9. / 13.]])
        min_mut_gap = 0
        for mut_collisions in ["drop", "merge"]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         mut_collisions,
                         expect_log_marginals,
                         expect_resps)

    def test_2pos_gap1_merge(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3]), np.array([4, 6])]
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.36, 0.04, 0.10, 0.24, 0.06]
        )
        expect_resps = np.ones((7, 1))
        min_mut_gap = 1
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     "merge",
                     expect_log_marginals,
                     expect_resps)

    def _setup_6pos_gap2_2clusters(self):
        """ Shared inputs for the 6 positions, 2 clusters, min_mut_gap=2
        tests where "drop" and "merge" produce different results. """
        p_mut = np.array([[0.05, 0.10],
                          [0.10, 0.20],
                          [0.15, 0.05],
                          [0.05, 0.15],
                          [0.20, 0.10],
                          [0.10, 0.05]])
        p_ends = np.zeros((6, 6))
        p_ends[0, 5] = 0.4
        p_ends[0, 4] = 0.2
        p_ends[1, 5] = 0.3
        p_ends[2, 5] = 0.1
        p_clust = np.array([0.4, 0.6])
        end5s = np.array([0, 0, 0, 0, 0, 1, 2, 0])
        end3s = np.array([5, 5, 5, 5, 4, 5, 5, 5])
        unmasked = np.array([0, 1, 2, 3, 4, 5])
        muts_per_pos = [np.array([7]),
                        np.array([2, 3]),
                        np.array([1, 6]),
                        np.array([3]),
                        np.array([2, 5]),
                        np.array([6])]
        min_mut_gap = 2
        return (p_mut, p_ends, p_clust, end5s, end3s,
                unmasked, muts_per_pos, min_mut_gap)

    def test_6pos_gap2_2clusters_drop(self):
        (p_mut, p_ends, p_clust, end5s, end3s,
         unmasked, muts_per_pos, min_mut_gap) = (
            self._setup_6pos_gap2_2clusters()
        )
        expect_log_marginals = np.log([0.2152940819,
                                       0.0219959898,
                                       0.0059803912,
                                       0.0062025729,
                                       0.1158307341,
                                       0.0289576835,
                                       0.0006179835,
                                       0.0188854458])
        expect_resps = np.array([[0.4000000000, 0.6000000000],
                                 [0.6909090909, 0.3090909091],
                                 [0.4000000000, 0.6000000000],
                                 [0.0811940299, 0.9188059701],
                                 [0.4130434783, 0.5869565217],
                                 [0.5869565217, 0.4130434783],
                                 [0.7989487516, 0.2010512484],
                                 [0.2400000000, 0.7600000000]])
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     "drop",
                     expect_log_marginals,
                     expect_resps)

    def test_6pos_gap2_2clusters_merge(self):
        (p_mut, p_ends, p_clust, end5s, end3s,
         unmasked, muts_per_pos, min_mut_gap) = (
            self._setup_6pos_gap2_2clusters()
        )
        expect_log_marginals = np.log([0.1988388,
                                       0.0251370,
                                       0.0074400,
                                       0.0083350,
                                       0.1069776,
                                       0.0331200,
                                       0.0007500,
                                       0.0174420])
        expect_resps = np.array([[0.4000000000, 0.6000000000],
                                 [0.6530612245, 0.3469387755],
                                 [0.3870967742, 0.6129032258],
                                 [0.0767846431, 0.9232153569],
                                 [0.4130434783, 0.5869565217],
                                 [0.5869565217, 0.4130434783],
                                 [0.8000000000, 0.2000000000],
                                 [0.2400000000, 0.7600000000]])
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     "merge",
                     expect_log_marginals,
                     expect_resps)


class TestCalcMarginalSimulated(ut.TestCase):
    """ Verify that calc_marginal returns log probabilities that match
    the empirical frequencies of simulated reads to within statistical
    tolerance, for both mut_collisions="drop" and "merge". """

    @staticmethod
    def _empirical_counts(muts: np.ndarray,
                          end5s: np.ndarray,
                          end3s: np.ndarray):
        """ Return (unique_end5s, unique_end3s, unique_muts, counts)
        for the unique (end5, end3, mutation pattern) tuples. """
        keys = [(int(a), int(b), tuple(row.tolist()))
                for a, b, row in zip(end5s, end3s, muts, strict=True)]
        counts = Counter(keys)
        unique = list(counts.keys())
        u_end5s = np.array([k[0] for k in unique])
        u_end3s = np.array([k[1] for k in unique])
        u_muts = np.array([k[2] for k in unique], dtype=bool)
        u_counts = np.array([counts[k] for k in unique])
        return u_end5s, u_end3s, u_muts, u_counts

    def _check(self,
               mut_collisions: str,
               p_mut_1d: np.ndarray,
               p_ends: np.ndarray,
               min_mut_gap: int,
               n_reads: int,
               z_threshold: float,
               seed: int):
        # Seed the simulator's RNG for reproducibility.
        unbias_test.rng = np.random.default_rng(seed=seed)
        n_pos = p_mut_1d.size
        muts, _, end5s, end3s = simulate_reads(n_reads, p_mut_1d, p_ends)
        # Process the simulated reads according to the collision policy.
        if mut_collisions == "drop":
            keep = label_no_close_muts(muts, min_mut_gap)
            proc_muts = muts[keep]
            proc_end5s = end5s[keep]
            proc_end3s = end3s[keep]
        else:
            proc_muts = merge_mutations_right_to_left(muts, min_mut_gap)
            proc_end5s = end5s
            proc_end3s = end3s
        n_total = proc_end5s.size
        self.assertGreater(n_total, 0)
        # Aggregate identical processed reads.
        u_end5s, u_end3s, u_muts, u_counts = self._empirical_counts(
            proc_muts, proc_end5s, proc_end3s
        )
        unmasked = np.arange(n_pos)
        muts_per_pos = [np.flatnonzero(u_muts[:, j]) for j in range(n_pos)]
        # Compute the model probability of each unique read.
        p_mut = p_mut_1d[:, np.newaxis]
        p_clust = np.array([1.0])
        log_marginals = calc_marginal(p_mut,
                                      p_ends,
                                      p_clust,
                                      u_end5s,
                                      u_end3s,
                                      unmasked,
                                      muts_per_pos,
                                      min_mut_gap,
                                      mut_collisions)
        expected = np.exp(log_marginals)
        # The probabilities of all possible reads must sum to 1.
        self.assertAlmostEqual(float(expected.sum()), 1.0, places=6)
        # Compare to empirical frequencies via a normal approximation.
        # Restrict to bins with enough expected counts for the normal
        # approximation to be valid.
        empirical = u_counts / n_total
        check = expected * n_total >= 10
        self.assertTrue(check.any(),
                        "no bins have enough counts to test")
        se = np.sqrt(expected[check] * (1.0 - expected[check]) / n_total)
        z = (empirical[check] - expected[check]) / se
        self.assertLess(float(np.max(np.abs(z))), z_threshold)

    def test_5pos_gap2_drop(self):
        p_mut_1d = np.array([0.10, 0.05, 0.20, 0.05, 0.15])
        p_ends = np.zeros((5, 5))
        p_ends[0, 4] = 0.40
        p_ends[0, 3] = 0.20
        p_ends[1, 4] = 0.25
        p_ends[2, 4] = 0.10
        p_ends[1, 3] = 0.05
        self._check(mut_collisions="drop",
                    p_mut_1d=p_mut_1d,
                    p_ends=p_ends,
                    min_mut_gap=2,
                    n_reads=1_000_000,
                    z_threshold=5.0,
                    seed=42)

    def test_5pos_gap2_merge(self):
        p_mut_1d = np.array([0.10, 0.05, 0.20, 0.05, 0.15])
        p_ends = np.zeros((5, 5))
        p_ends[0, 4] = 0.40
        p_ends[0, 3] = 0.20
        p_ends[1, 4] = 0.25
        p_ends[2, 4] = 0.10
        p_ends[1, 3] = 0.05
        self._check(mut_collisions="merge",
                    p_mut_1d=p_mut_1d,
                    p_ends=p_ends,
                    min_mut_gap=2,
                    n_reads=1_000_000,
                    z_threshold=5.0,
                    seed=42)


if __name__ == "__main__":
    ut.main(verbosity=2)
