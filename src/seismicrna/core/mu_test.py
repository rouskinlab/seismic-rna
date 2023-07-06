"""
Core -- Mutation Module
========================================================================
Auth: Matty

Unit tests for `core.mu`.
"""

import unittest as ut

import numpy as np

from .mu import calc_mu_adj, calc_mu_obs


def has_close_muts(bitvec: np.ndarray, min_gap: int):
    """ Return True if the bit vector has two mutations separated by
    fewer than `min_gap` non-mutated bits, otherwise False. """
    if bitvec.ndim != 1:
        raise ValueError(f"bitvec must have 1 dimension, but got {bitvec.ndim}")
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    if min_gap == 0:
        # Two mutations cannot be separated by fewer than 0 positions.
        return False
    # Close mutations are separated from each other by less than min_gap
    # non-mutated bits. Equivalently, the difference in their positions
    # (their distance) is < min_gap + 1, or ≤ min_gap. These distances
    # are computed as the differences (using np.diff) in the positions
    # of consecutive mutations (using np.flatnonzero).
    dists = np.diff(np.flatnonzero(bitvec))
    return dists.size > 0 and np.min(dists) <= min_gap


def label_close_muts(bitvecs: np.ndarray, min_gap: int):
    """ Return a 1D vector that is True for every row in `bitvecs` that
    has two mutations that are too close, and otherwise False. """
    if bitvecs.ndim != 2:
        raise ValueError(
            f"bitvects must have 2 dimensions, but got {bitvecs.ndim}")
    return np.array([has_close_muts(bitvec, min_gap) for bitvec in bitvecs])


def drop_close_muts(bitvecs: np.ndarray, min_gap: int):
    """ Return a new array without every row in `bitvecs` that has two
    mutations that are too close. """
    return bitvecs[np.logical_not(label_close_muts(bitvecs, min_gap))]


class TestHasCloseMuts(ut.TestCase):
    """ Test function `has_close_muts`. """

    def test_no_muts(self):
        """ Test that a bit vector with no mutations returns False. """
        bitvec = np.zeros(5, dtype=bool)
        for g in range(bitvec.size):
            self.assertFalse(has_close_muts(bitvec, g))

    def test_one_mut(self):
        """ Test that a bit vector with one mutation returns False. """
        for i in range(5):
            bitvec = np.zeros(5, dtype=bool)
            bitvec[i] = 1
            for g in range(bitvec.size):
                self.assertFalse(has_close_muts(bitvec, g))

    def test_more_muts(self):
        """ Test every bit vector with 2 - 5 mutations. """
        bitvecs_gaps = {
            # 2 mutations (n = 10)
            (0, 0, 0, 1, 1): 0,
            (0, 0, 1, 0, 1): 1,
            (0, 0, 1, 1, 0): 0,
            (0, 1, 0, 0, 1): 2,
            (0, 1, 0, 1, 0): 1,
            (0, 1, 1, 0, 0): 0,
            (1, 0, 0, 0, 1): 3,
            (1, 0, 0, 1, 0): 2,
            (1, 0, 1, 0, 0): 1,
            (1, 1, 0, 0, 0): 0,
            # 3 mutations (n = 10)
            (0, 0, 1, 1, 1): 0,
            (0, 1, 0, 1, 1): 0,
            (0, 1, 1, 0, 1): 0,
            (0, 1, 1, 1, 0): 0,
            (1, 0, 0, 1, 1): 0,
            (1, 0, 1, 0, 1): 1,
            (1, 0, 1, 1, 0): 0,
            (1, 1, 0, 0, 1): 0,
            (1, 1, 0, 1, 0): 0,
            (1, 1, 1, 0, 0): 0,
            # 4 mutations (n = 5)
            (0, 1, 1, 1, 1): 0,
            (1, 0, 1, 1, 1): 0,
            (1, 1, 0, 1, 1): 0,
            (1, 1, 1, 0, 1): 0,
            (1, 1, 1, 1, 0): 0,
            # 5 mutations (n = 1)
            (1, 1, 1, 1, 1): 0,
        }
        for bitvec, bit_gap in bitvecs_gaps.items():
            for min_gap in range(len(bitvec)):
                self.assertEqual(has_close_muts(np.array(bitvec), min_gap),
                                 bit_gap < min_gap)


class TestCalcMuAdj(ut.TestCase):
    """ Test function `mu.calc_mu_adj`. """

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random observed mutation rates.
                    mus_obs = np.random.default_rng().random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_adj_sim = calc_mu_adj(mus_obs, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_adj_sep = np.empty_like(mus_obs)
                    for i in range(k):
                        mus_obs_i = mus_obs[:, i].reshape((n_pos, 1))
                        mus_adj_i = calc_mu_adj(mus_obs_i, g).reshape(n_pos)
                        mus_adj_sep[:, i] = mus_adj_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_adj_sim, mus_adj_sep))

    def test_inv_calc_mu_obs(self):
        """ Test that this function inverts `mu.calc_mu_obs`. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Generate random real mutation rates.
            mus = np.random.default_rng().random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the observed mutation rates.
                    mus_obs = calc_mu_obs(mus, g)
                    # Adjust the observed mutation rates.
                    mus_adj = calc_mu_adj(mus_obs, g)
                    # Test if adjusted and initial mutation rates match.
                    self.assertTrue(np.allclose(mus_adj, mus))


class TestCalcMuObs(ut.TestCase):
    """ Test function `mu.calc_mu_obs`. """

    @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_obs_empirical(self):
        """ Test that this function accurately predicts the mutation
        rates that are actually observed when simulated bit vectors are
        filtered to remove mutations that are too close. """
        n_pos = 10
        min_g, max_g = 0, 4
        min_m, max_m = 0.01, 0.1
        # Choose the number of vectors to simulate as follows:
        # The number of mutations at each position in the simulated bit
        # vectors follows a binomial distribution, whose std. dev. is
        # sqrt(p * (1 - p) * n).
        # The proportion of mutations is this quantity divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of vectors) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.allclose
        n_vec = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Generate random real mutation rates.
        mus = min_m + np.random.default_rng().random(n_pos) * (max_m - min_m)
        # Generate random bit vectors with the expected mutation rates.
        bvecs = None
        while bvecs is None or not np.allclose(np.mean(bvecs, axis=0), mus,
                                               atol=tol, rtol=0.):
            bvecs = np.less(np.random.default_rng().random((n_vec, n_pos)), mus)
        # Test each minimum gap between mutations (g).
        for g in range(min_g, max_g + 1):
            with self.subTest(g=g):
                # Drop bit vectors with mutations too close.
                bvecs_g = drop_close_muts(bvecs, g)
                # Compute the empirically observed mutation rates.
                mus_obs_emp = np.mean(bvecs_g, axis=0)
                # Predict the observed mutation rates with calc_mu_obs.
                mus_obs_prd = calc_mu_obs(mus.reshape((-1, 1)), g).reshape(-1)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.allclose(mus_obs_emp, mus_obs_prd,
                                            atol=tol, rtol=0.))

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random real mutation rates.
                    mus = np.random.default_rng().random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_obs_sim = calc_mu_obs(mus, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_obs_sep = np.empty_like(mus)
                    for i in range(k):
                        mus_i = mus[:, i].reshape((n_pos, 1))
                        mus_obs_i = calc_mu_obs(mus_i, g).reshape(n_pos)
                        mus_obs_sep[:, i] = mus_obs_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_obs_sim, mus_obs_sep))

    def test_inv_calc_mu_adj(self):
        """ Test that this function inverts `mu.calc_mu_adj`. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Generate random observed mutation rates.
            mus_obs = np.random.default_rng().random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the adjusted mutation rates.
                    mus_adj = calc_mu_adj(mus_obs, g)
                    # Recompute the observed mutation rates.
                    mus_reobs = calc_mu_obs(mus_adj, g)
                    # Compare observed and reobserved mutation rates.
                    self.assertTrue(np.allclose(mus_obs, mus_reobs))
