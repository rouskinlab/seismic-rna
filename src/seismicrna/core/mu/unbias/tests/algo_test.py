import unittest as ut
from itertools import product

import numpy as np

from seismicrna.core.mu.unbias.algo import (_calc_p_noclose_given_ends,
                                            _calc_p_mut_given_span_noclose,
                                            calc_p_mut_p_ends_numpy,
                                            calc_p_noclose_given_ends_numpy,
                                            _calc_p_nomut_window,
                                            _clip,
                                            _adjust_min_gap,
                                            _triu_sum,
                                            _triu_norm,
                                            _triu_div,
                                            _triu_allclose)

rng = np.random.default_rng()


def no_close_muts(read: np.ndarray, min_gap: int):
    """ Return True if the read has no two mutations separated by fewer
    than `min_gap` non-mutated positions, otherwise False. """
    n_pos, = read.shape
    if min_gap <= 0 or n_pos <= 1:
        return True
    # Close mutations are separated from each other by less than min_gap
    # non-mutated bits. Equivalently, the difference in their positions
    # (their distance) is < min_gap + 1, or ≤ min_gap. These distances
    # are computed as the differences (using np.diff) in the positions
    # of consecutive mutations (using np.flatnonzero).
    diffs = np.diff(np.flatnonzero(read))
    return diffs.size == 0 or np.min(diffs) > min_gap


def label_no_close_muts(muts: np.ndarray, min_gap: int):
    """ Return a 1D vector that is True for every row in `muts` that
    has no two mutations that are too close, otherwise False. """
    n_reads, n_pos = muts.shape
    return np.fromiter((no_close_muts(muts[i], min_gap)
                        for i in range(n_reads)),
                       dtype=bool)


def simulate_reads(n_reads: int, p_mut: np.ndarray, p_ends: np.ndarray):
    """ Simulate `n_reads` reads based on the mutation rates (`p_mut`)
    and the distributions of end coordinates (`p_ends`). """
    n_pos, = p_mut.shape
    if p_ends.shape != (n_pos, n_pos):
        raise ValueError(f"p_ends must have dimensions {n_pos, n_pos}, "
                         f"but got {p_ends.shape}")
    # To sample end coordinates, first find the row and column indexes
    # of the upper triangle of the matrix.
    rows, cols = np.triu_indices(n_pos)
    # Extract the end coordinate probabilities from the upper triangle.
    p_ends_triu = p_ends[rows, cols]
    # Sample the end coordinates with replacement.
    end_indexes = rng.choice(p_ends_triu.size,
                             n_reads,
                             p=p_ends_triu,
                             replace=True)
    end5s = rows[end_indexes]
    end3s = cols[end_indexes]
    # Make a boolean array like `muts` where each element indicates
    # whether the read covers the position.
    info = np.empty((n_reads, n_pos))
    for pos in range(n_pos):
        info[:, pos] = np.logical_and(end5s <= pos, pos <= end3s)
    # Simulate reads as a boolean array where each row is a read, each
    # column is a position, and each element indicates whether the read
    # has a mutation at the position.
    muts = np.logical_and(np.less(rng.random((n_reads, n_pos)), p_mut), info)
    return muts, info, end5s, end3s


def simulate_params(n_pos: int, n_cls: int, p_mut_max: float = 1.):
    """ Return `p_mut` and `p_ends` parameters for `simulate_reads`. """
    p_mut = p_mut_max * rng.random((n_pos, n_cls))
    p_ends = np.triu(1. - rng.random((n_pos, n_pos)))
    p_ends /= np.sum(p_ends)
    return p_mut, p_ends


class TestNoCloseMuts(ut.TestCase):

    def test_no_muts(self):
        n_pos = 5
        bitvec = np.zeros(n_pos, dtype=bool)
        for g in range(bitvec.size):
            self.assertTrue(no_close_muts(bitvec, g))

    def test_one_mut(self):
        n_pos = 5
        for i in range(n_pos):
            bitvec = np.zeros(n_pos, dtype=bool)
            bitvec[i] = 1
            for g in range(bitvec.size):
                self.assertTrue(no_close_muts(bitvec, g))

    def test_more_muts(self):
        mut_gaps = {
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
        for read, mut_gap in mut_gaps.items():
            for min_gap in range(len(read)):
                self.assertEqual(no_close_muts(np.array(read), min_gap),
                                 mut_gap >= min_gap)


class TestClip(ut.TestCase):

    def test_with_clip(self):
        n_pos = 64
        n_nan = 16
        min_scale = 1
        max_scale = 10
        nan_indexes = rng.choice(n_pos, n_nan, replace=False)
        for scale in range(min_scale, max_scale + 1):
            # Generate random mutation rates, some of which may not be
            # in the range [0, 1].
            mus = scale * rng.random(n_pos, dtype=float) - (scale - 1) / 2.
            # Set several values to NaN to ensure that clip can also
            # replace missing values with 0.
            mus[nan_indexes] = np.nan
            # Confirm that the values are NaN.
            self.assertEqual(np.count_nonzero(np.isnan(mus)), n_nan)
            clipped = _clip(mus)
            # Test that all clipped mutation rates are in [0, 1].
            self.assertTrue(np.all(clipped >= 0.) and np.all(clipped <= 1.))
            # Test that NaN values in mus become 0 values in clipped.
            self.assertEqual(np.count_nonzero(clipped[nan_indexes]), 0)
            self.assertFalse(np.any(np.isnan(clipped)))
            self.assertTrue(np.all(np.isnan(mus[nan_indexes])))

    def test_without_clip(self):
        mus = rng.random(64, dtype=float)
        # All values are in [0, 1] and so should not be clipped.
        self.assertTrue(np.allclose(mus, _clip(mus)))


class TestTriuAllClose(ut.TestCase):

    def test_equal(self):
        for npos in range(4):
            for extra in range(3):
                for extras in product(range(4), repeat=extra):
                    dims = (npos, npos) + extras
                    a = rng.random(dims)
                    self.assertTrue(_triu_allclose(a, a.copy()))

    def test_triu(self):
        for npos in range(2, 4):
            a = rng.random((npos, npos))
            self.assertTrue(_triu_allclose(a, np.triu(a)))

    def test_tril(self):
        for npos in range(2, 4):
            a = rng.random((npos, npos))
            self.assertFalse(_triu_allclose(a, np.tril(a)))


class TestTriuDiv(ut.TestCase):

    def test_1x1(self):
        n = rng.random()
        d = 1. - rng.random()
        expect = np.array([[n / d]])
        self.assertTrue(np.array_equal(_triu_div(np.array([[n]]),
                                                 np.array([[d]])),
                                       expect))

    def test_1x1x1(self):
        n = rng.random()
        d = 1. - rng.random()
        expect = np.array([[[n / d]]])
        self.assertTrue(np.array_equal(_triu_div(np.array([[[n]]]),
                                                 np.array([[[d]]])),
                                       expect))

    def test_2x2(self):
        numer = np.array([[12., 3.],
                          [20., 56.]])
        denom = np.array([[2., 3.],
                          [5., 7.]])
        expect = np.array([[6., 1.],
                           [np.nan, 8.]])
        self.assertTrue(_triu_allclose(_triu_div(numer, denom), expect))

    def test_2x2x2(self):
        numer = np.array([[[2., 20.], [2., 36.]],
                          [[3., 7.], [12., 40.]]])
        denom = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[2., 4.], [1., 6.]],
                           [[np.nan, np.nan], [3., 5.]]])
        self.assertTrue(_triu_allclose(_triu_div(numer, denom), expect))


class TestTriuSum(ut.TestCase):

    def test_all_zero(self):
        for ndim in range(2, 6):
            array = rng.random((0,) * ndim)
            expect = np.zeros((0,) * (ndim - 2))
            self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_1x1(self):
        x = rng.random()
        array = np.array([[x]])
        expect = np.array(x)
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_1x1x1(self):
        x = rng.random()
        array = np.array([[[x]]])
        expect = np.array([x])
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_1x1x2(self):
        x = rng.random()
        y = rng.random()
        array = np.array([[[x, y]]])
        expect = np.array([x, y])
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array(7.)
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([7.])
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([7., 19.])
        self.assertTrue(np.array_equal(_triu_sum(array), expect))

    def test_2x2x2x2(self):
        array = np.array([[[[1., 5.],
                            [9., 13.]], [[2., 6.],
                                         [10., 14.]]],
                          [[[3., 7.],
                            [11., 15.]], [[4., 8.],
                                          [12., 16.]]]])
        expect = np.array([[7., 19.],
                           [31., 43.]])
        self.assertTrue(np.array_equal(_triu_sum(array), expect))


class TestTriuNorm(ut.TestCase):

    def test_1x1(self):
        array = rng.random((1, 1))
        expect = np.ones_like(array)
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_1x1x1(self):
        array = rng.random((1, 1, 1))
        expect = np.ones_like(array)
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_1x1x2(self):
        array = rng.random((1, 1, 2))
        expect = np.ones_like(array)
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array([[1. / 7., 2. / 7.],
                           [3. / 7., 4. / 7.]])
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([[[1. / 7.], [2. / 7.]],
                           [[3. / 7.], [4. / 7.]]])
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[1. / 7., 5. / 19.], [2. / 7., 6. / 19.]],
                           [[3. / 7., 7. / 19.], [4. / 7., 8. / 19.]]])
        self.assertTrue(np.array_equal(_triu_norm(array), expect))

    def test_2x2x2x2(self):
        array = np.array([[[[1., 5.],
                            [9., 13.]], [[2., 6.],
                                         [10., 14.]]],
                          [[[3., 7.],
                            [11., 15.]], [[4., 8.],
                                          [12., 16.]]]])
        expect = np.array([[[[1. / 7., 5. / 19.],
                             [9. / 31., 13. / 43]], [[2. / 7., 6. / 19.],
                                                     [10. / 31, 14. / 43]]],
                           [[[3. / 7., 7. / 19.],
                             [11. / 31, 15. / 43]], [[4. / 7., 8. / 19.],
                                                     [12. / 31, 16. / 43]]]])
        self.assertTrue(np.array_equal(_triu_norm(array), expect))


class TestAdjustMinGap(ut.TestCase):

    def test_zero_le_gap_lt_npos(self):
        for npos in range(10):
            for gap in range(npos):
                with self.subTest(npos=npos, gap=gap):
                    self.assertEqual(_adjust_min_gap(npos, gap), gap)

    def test_zero_lt_npos_le_gap(self):
        for gap in range(10):
            for npos in range(1, gap + 1):
                with self.subTest(npos=npos, gap=gap):
                    self.assertEqual(_adjust_min_gap(npos, gap), npos - 1)

    def test_npos_le_zero_le_gap(self):
        for npos in range(-4, 1):
            for gap in range(5):
                with self.subTest(npos=npos, gap=gap):
                    self.assertEqual(_adjust_min_gap(npos, gap), 0)

    def test_gap_le_zero_le_npos(self):
        for npos in range(5):
            for gap in range(-4, 1):
                with self.subTest(npos=npos, gap=gap):
                    self.assertEqual(_adjust_min_gap(npos, gap), 0)


class TestCalcPNoCloseGivenEnds(ut.TestCase):

    def test_min_gap_0(self):
        for npos in range(10):
            with self.subTest(npos=npos):
                mu = rng.random((npos, 1))
                nm = _calc_p_nomut_window(mu, 0)
                nm_expect = np.ones((1, npos + 1, 1))
                self.assertTrue(np.array_equal(nm, nm_expect))
                nc = _calc_p_noclose_given_ends(mu, nm)
                nc_expect = np.ones((npos, npos, 1))
                self.assertTrue(np.array_equal(nc, nc_expect))

    def test_length_0(self):
        mu = np.zeros((0, 1))
        for min_gap in range(4):
            with self.subTest(min_gap=min_gap):
                nm = _calc_p_nomut_window(mu, min_gap)
                nm_expect = np.ones((1, 1, 1))
                self.assertTrue(np.array_equal(nm, nm_expect))
                nc = _calc_p_noclose_given_ends(mu, nm)
                nc_expect = np.ones((0, 0, 1))
                self.assertTrue(np.array_equal(nc, nc_expect))

    def test_length_1(self):
        for x in [0., 0.01, 0.1, 1.]:
            mu = np.array([[x]])
            for min_gap in range(4):
                with self.subTest(x=x, min_gap=min_gap):
                    nm = _calc_p_nomut_window(mu, min_gap)
                    nm_expect = np.ones((1, 2, 1))
                    self.assertTrue(np.array_equal(nm, nm_expect))
                    nc = _calc_p_noclose_given_ends(mu, nm)
                    nc_expect = np.ones((1, 1, 1))
                    self.assertTrue(np.array_equal(nc, nc_expect))

    def _check_1_cluster(self, x, x_expect):
        self.assertEqual(x.shape, x_expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(x, x_expect),
                                             np.isnan(x_expect))))

    def test_length_2_min_gap_1(self):
        mu = np.array([[0.1, 0.2]]).reshape((2, 1))
        min_gap = 1
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1.],
            [np.nan, 0.9, 0.8],
        ]).reshape((2, 3, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98],
            [np.nan, 1.],
        ]).reshape((2, 2, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_3_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 1
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
        ]).reshape((2, 4, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.926],
            [1., 1., 0.94],
            [np.nan, 1., 1.],
        ]).reshape((3, 3, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_3_min_gap_2(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 2
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
            [np.nan, np.nan, 0.72, 0.56],
        ]).reshape((3, 4, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.902],
            [1., 1., 0.94],
            [1., 1., 1.],
        ]).reshape((3, 3, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_4_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 1
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
        ]).reshape((2, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.926, 0.83],
            [1., 1., 0.94, 0.844],
            [np.nan, 1., 1., 0.88],
            [np.nan, np.nan, 1., 1.],
        ]).reshape((4, 4, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_4_min_gap_2(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 2
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
        ]).reshape((3, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.902, 0.7652],
            [1., 1., 0.94, 0.788],
            [1., 1., 1., 0.88],
            [np.nan, 1., 1., 1.],
        ]).reshape((4, 4, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_4_min_gap_3(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 3
        nm = _calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
            [np.nan, np.nan, np.nan, 0.504, 0.336],
        ]).reshape((4, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = _calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.902, 0.7428],
            [1., 1., 0.94, 0.788],
            [1., 1., 1., 0.88],
            [1., 1., 1., 1.],
        ]).reshape((4, 4, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_clusters(self):
        for ncls in range(4):
            for npos in range(5):
                for min_gap in range(npos):
                    with self.subTest(npos=npos, ncls=ncls, min_gap=min_gap):
                        mu = rng.random((npos, ncls))
                        nm = _calc_p_nomut_window(mu, min_gap)
                        self.assertEqual(nm.shape, (min_gap + 1, npos + 1, ncls))
                        nc = _calc_p_noclose_given_ends(mu, nm)
                        self.assertEqual(nc.shape, (npos, npos, ncls))
                        for k in range(ncls):
                            mu_k = mu[:, k].reshape((npos, 1))
                            nm_k = _calc_p_nomut_window(mu_k, min_gap)
                            self.assertEqual(nm_k.shape,
                                             (min_gap + 1, npos + 1, 1))
                            self.assertTrue(np.allclose(np.triu(nm[:, :, k]),
                                                        np.triu(nm_k[:, :, 0])))
                            nc_k = _calc_p_noclose_given_ends(mu_k, nm_k)
                            self.assertEqual(nc_k.shape,
                                             (npos, npos, 1))
                            self.assertTrue(np.allclose(np.triu(nc[:, :, k]),
                                                        np.triu(nc_k[:, :, 0])))


class TestCalcPNoCloseGivenEndsNumPy(ut.TestCase):

    def test_1_dim(self):
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
        for n_dim in range(5):
            if n_dim == 1 or n_dim == 2:
                # Skip the dimensions that are valid.
                continue
            err_msg = ("p_mut_given_span must have 1 or 2 dimensions, "
                       f"but got {n_dim}")
            for size in range(5):
                dims = (size,) * n_dim
                for gap in range(5):
                    self.assertRaisesRegex(ValueError,
                                           err_msg,
                                           calc_p_noclose_given_ends_numpy,
                                           rng.random(dims),
                                           gap)


class TestBiasCalculation(ut.TestCase):

    def test_simulated(self):
        from scipy.stats import binom
        confidence = 0.9999
        # Simulate reads.
        n_pos = 6
        p_mut_max = 1.0
        n_reads = 1_000_000
        p_mut, p_ends = simulate_params(n_pos, 1, p_mut_max)
        muts, info, end5s, end3s = simulate_reads(n_reads, p_mut[:, 0], p_ends)
        # Test each minimum gap between mutations (min_gap).
        for min_gap in [0, 3, n_pos]:
            with self.subTest(min_gap=min_gap):
                # Determine which reads have no two mutations too close.
                noclose = label_no_close_muts(muts, min_gap)
                # Calculate the theoretical probability of no mutations
                # in each window.
                p_nomut_window_theory = _calc_p_nomut_window(p_mut, min_gap)
                # Calculate the theoretical probability of no mutations
                # being too close for each pair of end coordinates.
                p_noclose_given_ends_theory = _calc_p_noclose_given_ends(
                    p_mut,
                    p_nomut_window_theory
                ).reshape((n_pos, n_pos))
                # Compare the simulated probability of no mutations
                # being too close for each pair of end coordinates.
                for end5 in range(n_pos):
                    for end3 in range(end5, n_pos):
                        if end5 == end3:
                            # The read must have no close mutations.
                            self.assertEqual(
                                p_noclose_given_ends_theory[end5, end3], 1.
                            )
                        else:
                            # Count the reads with those coordinates.
                            end53s = np.logical_and(end5s == end5,
                                                    end3s == end3)
                            # Find the fraction of simualted reads with
                            # those coordinates without close mutations.
                            p_noclose_given_ends_sim = (
                                    np.sum(np.logical_and(end53s, noclose))
                                    / np.sum(end53s)
                            )
                            # Compute the confidence interval for the
                            # proportion of reads.
                            n_expect = round(n_reads * p_ends[end5, end3])
                            p_expect = p_noclose_given_ends_theory[end5, end3]
                            inter_lo, inter_up = binom.interval(confidence,
                                                                n_expect,
                                                                p_expect)
                            # Ensure the interval based on theoretical
                            # probability agrees with the simulation.
                            self.assertGreaterEqual(p_noclose_given_ends_sim,
                                                    inter_lo / n_expect)
                            self.assertLessEqual(p_noclose_given_ends_sim,
                                                 inter_up / n_expect)
                # Compute the expected coverage at each position.
                n_expect = np.zeros(n_pos, dtype=int)
                for end5 in range(n_pos):
                    for end3 in range(end5, n_pos):
                        n_expect[end5: end3 + 1] += round(n_reads
                                                          * p_ends[end5, end3])
                # Calculate the theoretical probability of a mutation at
                # each position given no mutations are too close.
                p_mut_given_noclose_theory = _calc_p_mut_given_span_noclose(
                    p_mut,
                    p_ends,
                    p_noclose_given_ends_theory.reshape((n_pos, n_pos, 1)),
                    p_nomut_window_theory
                ).reshape(n_pos)
                # Compare the simulated probabilities of mutations given
                # no mutations are too close.
                p_mut_given_noclose_sim = (np.sum(muts[noclose], axis=0)
                                           / np.sum(info[noclose], axis=0))
                inter_lo, inter_up = binom.interval(confidence,
                                                    n_expect,
                                                    p_mut_given_noclose_theory)
                self.assertTrue(np.all(p_mut_given_noclose_sim
                                       >= inter_lo / n_expect))
                self.assertTrue(np.all(p_mut_given_noclose_sim
                                       <= inter_up / n_expect))


class TestCalcPMutGivenSpanNoClose(ut.TestCase):

    def test_mu_multiplex(self):
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
                    mus_adj = calc_p_mut_p_ends_numpy(mus_obs, g)
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
                    mus_adj_sim = calc_p_mut_p_ends_numpy(mus_obs, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_adj_sep = np.empty_like(mus_obs)
                    for i in range(k):
                        obs_i = mus_obs[:, i].reshape((n_pos, 1))
                        adj_i = calc_p_mut_p_ends_numpy(obs_i, g).reshape(n_pos)
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
                    mus_adj = calc_p_mut_p_ends_numpy(mus_obs, g)
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
