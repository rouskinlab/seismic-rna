"""

Tests for Mutation Rate Core Module

========================================================================

"""

import unittest as ut
from logging import getLogger, Filter, LogRecord

import numpy as np
import pandas as pd

from .. import mu
from ..mu import (MAX_MU, clip, _calc_mu_obs,
                  calc_mu_adj_numpy, calc_mu_adj_df, calc_mu_adj_series,
                  calc_f_obs_numpy, calc_f_obs_df, calc_f_obs_series,
                  get_mu_quantile, normalize, winsorize)
from ..sim import rng
from ..sect import seq_pos_to_index, Section
from ..seq import DNA

mu_logger = getLogger(mu.__name__)


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

        clip_filter = ClipFilter()
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
            mu_logger.addFilter(clip_filter)
            try:
                # Clip the mutation rates to the bounds.
                clipped = clip(mus)
            finally:
                # Re-enable the warnings for mu.clip().
                mu_logger.removeFilter(clip_filter)
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


class TestCalcDataFrame(ut.TestCase):
    """ Test functions `mu.calc_mu_adj_df` and mu.calc_f_obs_df. """

    def test_equals_numpy(self):
        """ Check if the output of `calc_mu_adj_df` equals that of
        `calc_mu_adj_numpy` """
        max_mu = 0.1
        start = 1
        gaps = [0, 3]
        for length in range(1, 10):
            # Generate a random reference sequence.
            refseq = DNA.random(length)
            # Make a section for the sequence.
            section = Section("myref", refseq)
            for n_pos in range(length):
                # Choose a random set of positions, and sort them.
                pos = np.sort(rng.choice(length, n_pos, replace=False))
                # Make an index from those positions.
                index = seq_pos_to_index(refseq, pos + start, start)
                for n_clust in range(5):
                    clusters = pd.Index([f"Cluster-{i}"
                                         for i in range(1, n_clust + 1)])
                    # Generate random mutation rates.
                    mus_obs_values = max_mu * rng.random((n_pos, n_clust))
                    mus_obs_df = pd.DataFrame(mus_obs_values,
                                              index=index,
                                              columns=clusters)
                    # To run calc_mu_adj_numpy, create an array of the
                    # mutation rates where values in missing positions
                    # are set to 0.
                    mus_obs_np = np.zeros((length, n_clust))
                    for i_value, i_numpy in enumerate(pos):
                        mus_obs_np[i_numpy] = mus_obs_values[i_value]
                    for gap in gaps:
                        # Run calc_mu_adj_df.
                        mus_adj_df = calc_mu_adj_df(mus_obs_df, section, gap)
                        # Run calc_mu_adj_numpy.
                        mus_adj_np = calc_mu_adj_numpy(mus_obs_np, gap)
                        # Compare the results.
                        self.assertIsInstance(mus_adj_df, pd.DataFrame)
                        self.assertTrue(np.allclose(mus_adj_df.values,
                                                    mus_adj_np[pos]))
                        self.assertTrue(index.equals(mus_adj_df.index))
                        self.assertTrue(clusters.equals(mus_adj_df.columns))
                        # Run calc_f_obs_df.
                        f_obs_df = calc_f_obs_df(mus_adj_df, section, gap)
                        # Run calc_f_obs_numpy.
                        f_obs_np = calc_f_obs_numpy(mus_adj_np, gap)
                        # Compare the results.
                        self.assertIsInstance(f_obs_df, pd.Series)
                        self.assertTrue(np.allclose(f_obs_df.values,
                                                    f_obs_np))
                        self.assertTrue(clusters.equals(f_obs_df.index))


class TestCalcSeries(ut.TestCase):
    """ Test `mu.calc_mu_adj_series` and mu.calc_f_obs_series. """

    def test_equals_numpy(self):
        """ Check if the output of `calc_mu_adj_df` equals that of
        `calc_mu_adj_numpy` """
        max_mu = 0.1
        start = 1
        gaps = [0, 3]
        for length in range(1, 10):
            # Generate a random reference sequence.
            refseq = DNA.random(length)
            # Make a section for the sequence.
            section = Section("myref", refseq)
            for n_pos in range(length):
                # Choose a random set of positions, and sort them.
                pos = np.sort(rng.choice(length, n_pos, replace=False))
                # Make an index from those positions.
                index = seq_pos_to_index(refseq, pos + start, start)
                # Generate random mutation rates.
                mus_obs_values = max_mu * rng.random(n_pos)
                mus_obs_series = pd.Series(mus_obs_values, index=index)
                # To run calc_mu_adj_numpy, create an array of the
                # mutation rates where values in missing positions
                # are set to 0.
                mus_obs_np = np.zeros(length)
                for i_value, i_numpy in enumerate(pos):
                    mus_obs_np[i_numpy] = mus_obs_values[i_value]
                for gap in gaps:
                    # Run calc_mu_adj_series.
                    mus_adj_series = calc_mu_adj_series(mus_obs_series,
                                                        section, gap)
                    # Run calc_mu_adj_numpy.
                    mus_adj_np = calc_mu_adj_numpy(mus_obs_np, gap)
                    # Compare the results.
                    self.assertIsInstance(mus_adj_series, pd.Series)
                    self.assertTrue(np.array_equal(mus_adj_series.values,
                                                   mus_adj_np[pos]))
                    self.assertTrue(index.equals(mus_adj_series.index))
                    # Run calc_f_obs_series.
                    f_obs_series = calc_f_obs_series(mus_adj_series,
                                                     section, gap)
                    # Run calc_f_obs_numpy.
                    f_obs_np = calc_f_obs_numpy(mus_adj_np, gap)
                    # Compare the results.
                    self.assertIsInstance(f_obs_series, float)
                    self.assertIsInstance(f_obs_np, float)
                    self.assertEqual(f_obs_series, f_obs_np)


class TestGetMuQuantile(ut.TestCase):
    """ Test the function `mu.get_mu_quantile`. """

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
                values = np.array([get_mu_quantile(mus, quantile)
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
                    values = np.array([get_mu_quantile(mus, quantile)
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
            mu_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_mu_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                mu_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_empty(self):
        """ Test that an empty array always returns NaN values. """
        # Make mus an empty array.
        mus = np.array([], dtype=float)
        for n in [5, 11, 19]:
            quantiles = np.linspace(0., 1., n)
            # Temporarily suppress warnings about NaN values.
            mu_logger.addFilter(self.nan_filter)
            try:
                values = np.array([get_mu_quantile(mus, quantile)
                                   for quantile in quantiles])
            finally:
                # Re-enable warnings about NaN values.
                mu_logger.removeFilter(self.nan_filter)
            # Test that all quantile values are NaN.
            self.assertTrue(np.all(np.isnan(values)))

    def test_invalid_quantiles(self):
        """ Test that invalid quantiles raise errors. """
        n = 11
        mus = rng.random(n)
        errmsg = "Quantiles must be in the range [[]0, 1[]]"
        # Test that negative quantiles are invalid.
        for quantile in np.linspace(0., -1., n)[1:]:
            self.assertRaisesRegex(ValueError, errmsg,
                                   get_mu_quantile, mus, quantile)
        # Test that quantiles greater than 1 are invalid.
        for quantile in np.linspace(1., 2., n)[1:]:
            self.assertRaisesRegex(ValueError, errmsg,
                                   get_mu_quantile, mus, quantile)
        # Test that NaN is an invalid quantile.
        self.assertRaisesRegex(ValueError, errmsg,
                               get_mu_quantile, mus, np.nan)


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
