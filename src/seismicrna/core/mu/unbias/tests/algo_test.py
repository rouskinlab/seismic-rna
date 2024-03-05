import unittest as ut
from logging import Filter, LogRecord

import numpy as np

from seismicrna.core.mu.unbias.algo import (_calc_p_noclose_given_ends,
                                            _calc_p_mut_given_span_noclose,
                                            calc_mu_adj_numpy,
                                            calc_p_noclose_given_ends_numpy,
                                            _clip,
                                            logger as algo_logger)

rng = np.random.default_rng()


def has_close_muts(read: np.ndarray, min_gap: int):
    """ Return True if the read has two mutations separated by fewer
    than `min_gap` non-mutated positions, otherwise False. """
    if read.ndim != 1:
        raise ValueError(f"read must have 1 dimension, but got {read.ndim}")
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
    dists = np.diff(np.flatnonzero(read))
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
        """ Test that values outside [0, 1] are clipped. """
        n_pos = 64
        n_nan = 16
        min_scale = 1
        max_scale = 10
        nan_indexes = rng.choice(n_pos, n_nan, replace=False)

        class ClipFilter(Filter):
            """ Suppress warnings about invalid mutation rates. """

            def filter(self, rec: LogRecord):
                """ Suppress warnings about invalid mutation rates. """
                msg = f"Mutation rates outside [0, 1]"
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
            algo_logger.addFilter(clip_filter := ClipFilter())
            try:
                # Clip the mutation rates to the bounds.
                clipped = _clip(mus)
            finally:
                # Re-enable the warnings for mu.clip().
                algo_logger.removeFilter(clip_filter)
            # Test that all clipped mutation rates are in [0, 1].
            self.assertTrue(np.all(clipped >= 0.) and np.all(clipped <= 1.))
            # Test that NaN values in mus become 0 values in clipped.
            self.assertEqual(np.count_nonzero(clipped[nan_indexes]), 0)
            self.assertFalse(np.any(np.isnan(clipped)))
            self.assertTrue(np.all(np.isnan(mus[nan_indexes])))

    def test_without_clip(self):
        """ Test that values inside [0, 1] are not clipped. """
        mus = rng.random(64, dtype=float)
        self.assertTrue(np.allclose(mus, _clip(mus)))


class TestCalcPNoCloseGivenEnds(ut.TestCase):

    def test_negative_min_gap(self):
        for min_gap in range(-3, 0):
            self.assertRaisesRegex(ValueError,
                                   f"min_gap must be ≥ 0, but got {min_gap}",
                                   _calc_p_noclose_given_ends,
                                   rng.random((1, 1)),
                                   min_gap)

    def test_min_gap_0(self):
        for npos in range(10):
            with self.subTest(npos=npos):
                d, w = _calc_p_noclose_given_ends(rng.random((npos, 1)), 0)
                d_expect = np.ones((npos, npos, 1))
                self.assertTrue(np.array_equal(d, d_expect))
                w_expect = np.ones((1, npos + 1, 1))
                self.assertTrue(np.array_equal(w, w_expect))

    def test_length_0(self):
        mu = np.zeros((0, 1))
        for min_gap in range(4):
            with self.subTest(min_gap=min_gap):
                d, w = _calc_p_noclose_given_ends(mu, min_gap)
                d_expect = np.ones((0, 0, 1))
                self.assertTrue(np.array_equal(d, d_expect))
                w_expect = np.ones((1, 1, 1))
                self.assertTrue(np.array_equal(w, w_expect))

    def test_length_1(self):
        for x in [0., 0.01, 0.1, 1.]:
            mu = np.array([[x]])
            for min_gap in range(4):
                with self.subTest(x=x, min_gap=min_gap):
                    d, w = _calc_p_noclose_given_ends(mu, min_gap)
                    d_expect = np.ones((1, 1, 1))
                    self.assertTrue(np.array_equal(d, d_expect))
                    w_expect = np.ones((1, 2, 1))
                    self.assertTrue(np.array_equal(w, w_expect))

    def _check_dw_1_cluster(self, d, w, d_expect, w_expect):
        self.assertEqual(d.shape, d_expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(d, d_expect),
                                             np.isnan(d_expect))))
        self.assertEqual(w.shape, w_expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(w, w_expect),
                                             np.isnan(w_expect))))

    def test_length_2_min_gap_1(self):
        mu = np.array([[0.1, 0.2]]).reshape((2, 1))
        min_gap = 1
        d_expect = np.array([
            [1., 0.98],
            [np.nan, 1.],
        ]).reshape((2, 2, 1))
        w_expect = np.array([
            [1., 1., 1.],
            [np.nan, 0.9, 0.8],
        ]).reshape((2, 3, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_length_3_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 1
        d_expect = np.array([
            [1., 0.98, 0.926],
            [1., 1., 0.94],
            [np.nan, 1., 1.],
        ]).reshape((3, 3, 1))
        w_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
        ]).reshape((2, 4, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_length_3_min_gap_2(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 2
        d_expect = np.array([
            [1., 0.98, 0.902],
            [1., 1., 0.94],
            [1., 1., 1.],
        ]).reshape((3, 3, 1))
        w_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
            [np.nan, np.nan, 0.72, 0.56],
        ]).reshape((3, 4, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_length_4_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 1
        d_expect = np.array([
            [1., 0.98, 0.926, 0.83],
            [1., 1., 0.94, 0.844],
            [np.nan, 1., 1., 0.88],
            [np.nan, np.nan, 1., 1.],
        ]).reshape((4, 4, 1))
        w_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
        ]).reshape((2, 5, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_length_4_min_gap_2(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 2
        d_expect = np.array([
            [1., 0.98, 0.902, 0.7652],
            [1., 1., 0.94, 0.788],
            [1., 1., 1., 0.88],
            [np.nan, 1., 1., 1.],
        ]).reshape((4, 4, 1))
        w_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
        ]).reshape((3, 5, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_length_4_min_gap_3(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 3
        d_expect = np.array([
            [1., 0.98, 0.902, 0.7428],
            [1., 1., 0.94, 0.788],
            [1., 1., 1., 0.88],
            [1., 1., 1., 1.],
        ]).reshape((4, 4, 1))
        w_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
            [np.nan, np.nan, np.nan, 0.504, 0.336],
        ]).reshape((4, 5, 1))
        d, w = _calc_p_noclose_given_ends(mu, min_gap)
        self._check_dw_1_cluster(d, w, d_expect, w_expect)

    def test_clusters(self):
        for ncls in range(4):
            for npos in range(5):
                for min_gap in range(npos):
                    with self.subTest(npos=npos, ncls=ncls, min_gap=min_gap):
                        mu = rng.random((npos, ncls))
                        d, w = _calc_p_noclose_given_ends(mu, min_gap)
                        self.assertEqual(d.shape, (npos, npos, ncls))
                        self.assertEqual(w.shape, (min_gap + 1, npos + 1, ncls))
                        for k in range(ncls):
                            mu_k = mu[:, k].reshape((npos, 1))
                            d_k, w_k = _calc_p_noclose_given_ends(mu_k, min_gap)
                            self.assertEqual(d_k.shape,
                                             (npos, npos, 1))
                            self.assertEqual(w_k.shape,
                                             (min_gap + 1, npos + 1, 1))
                            self.assertTrue(np.allclose(np.triu(d[:, :, k]),
                                                        np.triu(d_k[:, :, 0])))
                            self.assertTrue(np.allclose(np.triu(w[:, :, k]),
                                                        np.triu(w_k[:, :, 0])))


class TestCalcPNoCloseGivenEndsNumPy(ut.TestCase):

    # @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_simulated(self):
        """ Test that this function accurately predicts the fraction of
        bit vectors without mutations that are too close. """
        n_pos = 16
        min_m, max_m = 0.01, 0.1
        # Choose the number of reads to simulate as follows:
        # The number of reads with mutations that are too close follows
        # a binomial distribution with std. dev. sqrt(p * (1 - p) * n).
        # The proportion of reads is the above divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of reads) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.isclose
        n_reads = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Choose random mutation rates.
        mus = min_m + rng.random(n_pos) * (max_m - min_m)
        # Generate random reads from the mutation rates.
        reads = np.less(rng.random((n_reads, n_pos)), mus)
        # Compute the mutation rates after simulating.
        mus_sim = np.mean(reads, axis=0)
        # Test each minimum gap between mutations (min_gap).
        for min_gap in [0, 3, n_pos]:
            with self.subTest(min_gap=min_gap):
                # Determine how many reads have mutations too close.
                n_close = sum(has_close_muts(reads[i], min_gap)
                              for i in range(n_reads))
                # Find the observed fraction of reads in the simulation.
                p_noclose_simulation = 1. - n_close / n_reads
                # Calculate the theoretical observed fraction.
                p_noclose_theory = calc_p_noclose_given_ends_numpy(mus_sim,
                                                                   min_gap)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.isclose(p_noclose_simulation,
                                           p_noclose_theory[0, -1],
                                           atol=tol,
                                           rtol=0.))

    def test_1_dim(self):
        """ Test that giving a 1D array returns a 2D array. """
        max_n = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for gap in range(max_g + 1):
                p_noclose_given_ends = calc_p_noclose_given_ends_numpy(
                    rng.random(n_pos), gap
                )
                self.assertIsInstance(p_noclose_given_ends, np.ndarray)
                self.assertEqual(p_noclose_given_ends.shape, (n_pos, n_pos))

    def test_2_dim(self):
        """ Test that giving a 2D array returns a 3D array. """
        max_n = 5
        max_c = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for n_clust in range(max_c + 1):
                for gap in range(max_g + 1):
                    p_noclose_given_ends = calc_p_noclose_given_ends_numpy(
                        rng.random((n_pos, n_clust)), gap
                    )
                    self.assertIsInstance(p_noclose_given_ends, np.ndarray)
                    self.assertEqual(p_noclose_given_ends.shape,
                                     (n_pos, n_pos, n_clust))

    def test_invalid_dim(self):
        """ Test that other dimensionalities raise an error. """
        for n_dim in range(5):
            if n_dim == 1 or n_dim == 2:
                # Skip the dimensions that are valid.
                continue
            err_msg = ("Expected p_mut_given_span to have 1 or 2 dimensions, "
                       f"but got {n_dim}")
            for size in range(5):
                dims = (size,) * n_dim
                for gap in range(5):
                    self.assertRaisesRegex(ValueError,
                                           err_msg,
                                           calc_p_noclose_given_ends_numpy,
                                           rng.random(dims),
                                           gap)


class TestCalcPMutGivenSpanNoClose(ut.TestCase):

    # @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_simulated(self):
        n_pos = 10
        min_m, max_m = 0.01, 0.1
        # Choose the number of reads to simulate as follows:
        # The number of reads with mutations that are too close follows
        # a binomial distribution with std. dev. sqrt(p * (1 - p) * n).
        # The proportion of reads is the above divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of reads) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.isclose
        n_reads = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Generate random real mutation rates.
        mus = min_m + rng.random(n_pos) * (max_m - min_m)
        # Generate random reads from the mutation rates.
        reads = np.less(rng.random((n_reads, n_pos)), mus)
        # Compute the mutation rates after simulating.
        mus_sim = np.mean(reads, axis=0).reshape((-1, 1))
        # Assume every read is full length.
        p_ends = np.zeros((n_pos, n_pos))
        p_ends[0, -1] = 1.
        # Test each minimum gap between mutations (min_gap).
        for min_gap in [0, 3, n_pos]:
            with self.subTest(min_gap=min_gap):
                # Drop reads with mutations too close.
                reads_gap = drop_close_muts(reads, min_gap)
                # Compute the mutation rates given no close mutations.
                mus_noclose = np.mean(reads_gap, axis=0)
                # Calculate the theoretical observed mutation rates.
                mus_theory = _calc_p_mut_given_span_noclose(
                    mus_sim,
                    p_ends,
                    *_calc_p_noclose_given_ends(mus_sim, min_gap)
                ).reshape(-1)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.allclose(mus_noclose,
                                            mus_theory,
                                            atol=tol,
                                            rtol=0.))

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


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
