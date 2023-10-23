"""

Tests for Mutation Rate Core Module

========================================================================

"""

import unittest as ut
from logging import Filter, LogRecord

import numpy as np

from ..unbias import (MAX_MU,
                      _calc_mu_obs,
                      calc_mu_adj_numpy,
                      calc_f_obs_numpy,
                      clip,
                      logger as unbias_logger)
from ...rand import rng


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
        n_pos = 5
        bitvec = np.zeros(n_pos, dtype=bool)
        for g in range(bitvec.size):
            self.assertFalse(has_close_muts(bitvec, g))

    def test_one_mut(self):
        """ Test that a bit vector with one mutation returns False. """
        n_pos = 5
        for i in range(n_pos):
            bitvec = np.zeros(n_pos, dtype=bool)
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


class TestClip(ut.TestCase):
    """ Test function `mu.clip`. """

    def test_with_clip(self):
        """ Test that values outside [0, MAX_MU] are clipped. """
        n_pos = 64
        n_nan = 16
        min_scale = 1
        max_scale = 10
        nan_indexes = rng.choice(n_pos, n_nan, replace=False)

        class ClipFilter(Filter):
            """ Suppress warnings about invalid mutation rates. """

            def filter(self, rec: LogRecord):
                """ Suppress warnings about invalid mutation rates. """
                msg = f"Mutation rates outside [0, {MAX_MU}]"
                return not rec.msg.startswith(msg)

        for scale in range(min_scale, max_scale + 1):
            # Generate random mutation rates, some of which may not be
            # in the range [0, 1].
            mus = scale * rng.random(n_pos, dtype=float) - (scale - 1) / 2.
            # Set several values to NaN to ensure that clip can also
            # replace missing values with 0.
            mus[nan_indexes] = np.nan
            # Confirm that the values are NaN.
            self.assertEqual(np.count_nonzero(np.isnan(mus)), n_nan)
            # Suppress the warnings that mu.clip() issues for mutation
            # rates being outside the bounds, since in this case they
            # are out of bounds deliberately.
            unbias_logger.addFilter(clip_filter := ClipFilter())
            try:
                # Clip the mutation rates to the bounds.
                clipped = clip(mus)
            finally:
                # Re-enable the warnings for mu.clip().
                unbias_logger.removeFilter(clip_filter)
            # Test that all clipped mutation rates are in [0, MAX_MU].
            self.assertTrue(np.all(clipped >= 0.) and np.all(clipped <= MAX_MU))
            # Test that NaN values in mus become 0 values in clipped.
            self.assertEqual(np.count_nonzero(clipped[nan_indexes]), 0)
            self.assertFalse(np.any(np.isnan(clipped)))
            self.assertTrue(np.all(np.isnan(mus[nan_indexes])))

    def test_without_clip(self):
        """ Test that values inside [0, MAX_MU] are not clipped. """
        mus = rng.random(64, dtype=float) * MAX_MU
        self.assertTrue(np.allclose(mus, clip(mus)))


class TestCalcFObsNumpy(ut.TestCase):
    """ Test function `mu.calc_f_obs_numpy`. """

    @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_obs_empirical(self):
        """ Test that this function accurately predicts the fraction of
        bit vectors without mutations that are too close. """
        n_pos = 16
        min_m, max_m = 0.01, 0.1
        # Choose the number of vectors to simulate as follows:
        # The number of bit vectors with mutations that are too close
        # follows a binomial distribution, whose std. dev. is
        # sqrt(p * (1 - p) * n).
        # The proportion of bit vectors is the above divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of vectors) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.isclose
        n_vec = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Choose random mutation rates.
        mus = min_m + rng.random(n_pos) * (max_m - min_m)
        # Test each minimum gap between mutations (g).
        for g in [0, 3, n_pos]:
            with self.subTest(g=g):
                # Generate random bit vectors with the expected mutation
                # rates and determine how many have mutations too close.
                n_close = sum(has_close_muts(np.less(rng.random(n_pos), mus), g)
                              for _ in range(n_vec))
                # Find the simulated observed fraction of bit vectors.
                f_obs_sim = 1. - n_close / n_vec
                # Predict the observed fraction.
                f_obs_prd = calc_f_obs_numpy(mus, g)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.isclose(f_obs_sim, f_obs_prd,
                                           atol=tol, rtol=0.))

    def test_f_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        max_k = 5
        max_g = 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random real mutation rates.
                    mus = rng.random((n_pos, k)) * max_m
                    # Compute the observed fractions simultaneously.
                    f_obs_sim = calc_f_obs_numpy(mus, g)
                    # Compute the fractions separately.
                    f_obs_sep = np.empty_like(f_obs_sim)
                    for i in range(k):
                        f_obs_sep[i] = calc_f_obs_numpy(mus[:, i], g)
                    # Compare the results.
                    self.assertTrue(np.allclose(f_obs_sim, f_obs_sep))

    def test_1_dim(self):
        """ Test that giving a 1D array returns a float. """
        max_n = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for gap in range(max_g + 1):
                self.assertIsInstance(calc_f_obs_numpy(np.zeros((n_pos,)), gap),
                                      float)

    def test_2_dim(self):
        """ Test that giving a 2D array returns a 1D array with a value
        for each cluster. """
        max_n = 5
        max_c = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for n_clust in range(max_c + 1):
                for gap in range(max_g + 1):
                    f_obs = calc_f_obs_numpy(np.zeros((n_pos, n_clust)), gap)
                    self.assertIsInstance(f_obs, np.ndarray)
                    self.assertEqual(f_obs.ndim, 1)
                    self.assertEqual(f_obs.shape, (n_clust,))

    def test_invalid_dim(self):
        """ Test that other dimensionalities raise an error. """
        for n_dim in range(5):
            if n_dim == 1 or n_dim == 2:
                # Skip the dimensions that are valid.
                continue
            err_msg = f"Expected 1 or 2 dimensions, but got {n_dim}"
            for size in range(5):
                dims = (size,) * n_dim
                for gap in range(5):
                    self.assertRaisesRegex(ValueError, err_msg,
                                           calc_f_obs_numpy,
                                           np.zeros(dims), gap)


class TestCalcMuObs(ut.TestCase):
    """ Test function `mu._calc_mu_obs`. """

    @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_obs_empirical(self):
        """ Test that this function accurately predicts the mutation
        rates that are actually observed when simulated bit vectors are
        filtered to remove mutations that are too close. """
        n_pos = 10
        min_m, max_m = 0.01, 0.1
        # Choose the number of vectors to simulate as follows:
        # The number of mutations at each position in the simulated bit
        # vectors follows a binomial distribution, whose std. dev. is
        # sqrt(p * (1 - p) * n).
        # The proportion of mutations is the above divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of vectors) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.allclose
        n_vec = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Generate random real mutation rates.
        mus = min_m + rng.random(n_pos) * (max_m - min_m)
        # Generate random bit vectors with the expected mutation rates.
        bvecs = None
        while bvecs is None or not np.allclose(np.mean(bvecs, axis=0), mus,
                                               atol=tol, rtol=0.):
            bvecs = np.less(rng.random((n_vec, n_pos)), mus)
        # Test each minimum gap between mutations (g).
        for g in [0, 3, n_pos]:
            with self.subTest(g=g):
                # Drop bit vectors with mutations too close.
                bvecs_g = drop_close_muts(bvecs, g)
                # Compute the empirically observed mutation rates.
                mus_obs_emp = np.mean(bvecs_g, axis=0)
                # Predict the observed mutation rates with calc_mu_obs.
                mus_obs_prd = _calc_mu_obs(mus.reshape((-1, 1)), g).reshape(-1)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.allclose(mus_obs_emp, mus_obs_prd,
                                            atol=tol, rtol=0.))

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        max_k = 5
        max_g = 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random real mutation rates.
                    mus = rng.random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_obs_sim = _calc_mu_obs(mus, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_obs_sep = np.empty_like(mus_obs_sim)
                    for i in range(k):
                        mus_i = mus[:, i].reshape((n_pos, 1))
                        mus_obs_i = _calc_mu_obs(mus_i, g).reshape(n_pos)
                        mus_obs_sep[:, i] = mus_obs_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_obs_sim, mus_obs_sep))

    def test_inv_calc_mu_adj(self):
        """ Test that this function inverts `mu.calc_mu_adj`. """
        n_pos = 16
        max_k = 5
        max_g = 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(max_k + 1):
            # Generate random observed mutation rates.
            mus_obs = rng.random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the adjusted mutation rates.
                    mus_adj = calc_mu_adj_numpy(mus_obs, g)
                    # Recompute the observed mutation rates.
                    mus_reobs = _calc_mu_obs(mus_adj, g)
                    # Compare observed and reobserved mutation rates.
                    self.assertTrue(np.allclose(mus_obs, mus_reobs))


class TestCalcMuAdjNumpy(ut.TestCase):
    """ Test function `mu.calc_mu_adj_numpy`. """

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        max_k = 5
        max_g = 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random observed mutation rates.
                    mus_obs = rng.random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_adj_sim = calc_mu_adj_numpy(mus_obs, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_adj_sep = np.empty_like(mus_obs)
                    for i in range(k):
                        obs_i = mus_obs[:, i].reshape((n_pos, 1))
                        adj_i = calc_mu_adj_numpy(obs_i, g).reshape(n_pos)
                        mus_adj_sep[:, i] = adj_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_adj_sim, mus_adj_sep))

    def test_inv_calc_mu_obs(self):
        """ Test that this function inverts `mu._calc_mu_obs`. """
        n_pos = 16
        max_k = 5
        max_g = 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(max_k + 1):
            # Generate random real mutation rates.
            mus = rng.random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the observed mutation rates.
                    mus_obs = _calc_mu_obs(mus, g)
                    # Adjust the observed mutation rates.
                    mus_adj = calc_mu_adj_numpy(mus_obs, g)
                    # Test if adjusted and initial mutation rates match.
                    self.assertTrue(np.allclose(mus_adj, mus))

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
