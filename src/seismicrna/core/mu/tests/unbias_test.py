import unittest as ut
from itertools import product

import numpy as np

from seismicrna.core.mu.unbias import (_clip,
                                       _normalize,
                                       _adjust_min_gap,
                                       _triu_log,
                                       _triu_sum,
                                       _triu_cumsum,
                                       _triu_norm,
                                       _triu_dot,
                                       _triu_div,
                                       _triu_allclose,
                                       _calc_rectangular_sum,
                                       _calc_p_ends_observed,
                                       _calc_p_nomut_window,
                                       _calc_p_noclose_given_ends,
                                       _calc_p_mut_given_span_noclose,
                                       _calc_p_mut_given_span,
                                       _calc_p_ends,
                                       _slice_p_ends,
                                       _find_split_positions,
                                       calc_p_clust,
                                       calc_p_noclose,
                                       calc_p_noclose_given_ends,
                                       calc_p_ends_given_noclose,
                                       calc_p_clust_given_noclose,
                                       calc_params)

rng = np.random.default_rng(seed=0)


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
    """ Return `p_mut`, `p_ends`, and `p_cls` parameters. """
    p_mut = p_mut_max * rng.random((n_pos, n_cls))
    p_ends = np.triu(1. - rng.random((n_pos, n_pos)))
    p_ends /= np.sum(p_ends)
    p_cls = 1. - rng.random(n_cls)
    p_cls /= np.sum(p_cls)
    return p_mut, p_ends, p_cls


class TestCalcPEndsGivenNoClose(ut.TestCase):

    def test_npos0_ncls1(self):
        p_ends = np.ones((0, 0))
        p_noclose_given_ends = np.ones((0, 0, 1))
        self.assertRaisesRegex(ValueError,
                               "Size of dimension 'positions' must be ≥ 1, "
                               "but got 0",
                               calc_p_ends_given_noclose,
                               p_ends,
                               p_noclose_given_ends)

    def test_npos1_ncls0(self):
        p_ends = np.ones((1, 1))
        p_noclose_given_ends = np.ones((1, 1, 0))
        self.assertRaisesRegex(ValueError,
                               "Size of dimension 'clusters' must be ≥ 1, "
                               "but got 0",
                               calc_p_ends_given_noclose,
                               p_ends,
                               p_noclose_given_ends)

    def test_npos1_ncls1(self):
        p_ends = np.ones((1, 1))
        expect = np.ones((1, 1, 1))
        for p in [0.01, 0.1, 1.0]:
            p_noclose_given_ends = np.full((1, 1, 1), p)
            result = calc_p_ends_given_noclose(p_ends, p_noclose_given_ends)
            self.assertTrue(np.array_equal(result, expect))

    def test_npos2_ncls1(self):
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_noclose_given_ends = np.array([[[0.9], [0.6]],
                                         [[0.0], [0.8]]])
        expect = np.array([[[3. / 12.], [5. / 12.]],
                           [[np.nan], [4. / 12.]]])
        result = calc_p_ends_given_noclose(p_ends, p_noclose_given_ends)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))


class TestCalcPClustGivenNoClose(ut.TestCase):

    def test_ncls1(self):
        p_clust = np.ones((1,))
        expect = np.ones((1,))
        for npos in range(1, 5):
            with self.subTest(npos=npos):
                p_ends = np.ones((npos, npos))
                p_noclose_given_ends = np.ones((npos, npos, 1))
                p_noclose = calc_p_noclose(p_ends, p_noclose_given_ends)
                result = calc_p_clust_given_noclose(p_clust, p_noclose)
                self.assertEqual(result.shape, expect.shape)
                self.assertTrue(np.allclose(result, expect))

    def test_ncls2(self):
        p_clust = np.ones((1,))
        expect = np.ones((1,))
        for npos in range(1, 5):
            with self.subTest(npos=npos):
                p_ends = np.ones((npos, npos))
                p_noclose_given_ends = np.ones((npos, npos, 1))
                p_noclose = calc_p_noclose(p_ends, p_noclose_given_ends)
                result = calc_p_clust_given_noclose(p_clust, p_noclose)
                self.assertEqual(result.shape, expect.shape)
                self.assertTrue(np.allclose(result, expect))


class TestNoCloseMuts(ut.TestCase):

    def test_no_muts(self):
        n_pos = 5
        read = np.zeros(n_pos, dtype=bool)
        for g in range(read.size):
            self.assertTrue(no_close_muts(read, g))

    def test_one_mut(self):
        n_pos = 5
        for i in range(n_pos):
            read = np.zeros(n_pos, dtype=bool)
            read[i] = 1
            for g in range(read.size):
                self.assertTrue(no_close_muts(read, g))

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


class TestNormalize(ut.TestCase):

    def test_sum_positive(self):
        for ndim in range(1, 4):
            for dims in product(range(5), repeat=ndim):
                x = 1. - rng.random(dims)
                result = _normalize(x)
                self.assertEqual(result.shape, x.shape)
                if x.size > 0:
                    self.assertTrue(np.isclose(result.sum(), 1.))
                    ratio = x / result
                    self.assertTrue(np.allclose(ratio, ratio[(0,) * ndim]))

    def test_sum_zero(self):
        for ndim in range(1, 4):
            for dims in product(range(5), repeat=ndim):
                x = np.zeros(dims)
                result = _normalize(x)
                self.assertEqual(result.shape, x.shape)
                if x.size > 0:
                    self.assertTrue(np.isclose(result.sum(), 1.))
                    ratio = x / result
                    self.assertTrue(np.allclose(ratio, ratio[(0,) * ndim]))


class TestTriuLog(ut.TestCase):

    def compare(self, result: np.ndarray, expect: np.ndarray):
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))

    def test_2d(self):
        a = np.array([[1., 2.],
                      [0., 3.]])
        expect = np.array([[0., np.log(2.)],
                           [np.nan, np.log(3.)]])
        result = _triu_log(a)
        self.compare(result, expect)

    def test_3d(self):
        a = np.array([[[1., 4.], [2., 5.]],
                      [[0., 1.], [3., 6.]]])
        expect = np.array([[[0., np.log(4.)], [np.log(2.), np.log(5.)]],
                           [[np.nan, np.nan], [np.log(3.), np.log(6.)]]])
        result = _triu_log(a)
        self.compare(result, expect)


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


class TestTriuDot(ut.TestCase):

    def test_1x1(self):
        a = np.array([[2.]])
        b = np.array([[4.]])
        expect = np.array(8.)
        self.assertTrue(np.array_equal(_triu_dot(a, b), expect))

    def test_1x1x1(self):
        a = np.array([[[2.]]])
        b = np.array([[[4.]]])
        expect = np.array([8.])
        self.assertTrue(np.array_equal(_triu_dot(a, b), expect))

    def test_2x2(self):
        a = np.array([[2., 3.],
                      [5., 7.]])
        b = np.array([[4., 8.],
                      [16., 32.]])
        expect = np.array(256.)
        self.assertTrue(np.array_equal(_triu_dot(a, b), expect))

    def test_2x2x2(self):
        a = np.array([[[2., 20.], [3., 30.]],
                      [[5., 50.], [7., 70.]]])
        b = np.array([[[4., 40.], [8., 80.]],
                      [[16., 160.], [32., 320.]]])
        expect = np.array([256., 25600.])
        self.assertTrue(np.array_equal(_triu_dot(a, b), expect))


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


class TestTriuCumSum(ut.TestCase):

    def test_all_0(self):
        for ndim in range(2, 6):
            array = rng.random((0,) * ndim)
            result = _triu_cumsum(array)
            self.assertEqual(result.shape, array.shape)
            self.assertTrue(_triu_allclose(result, array))

    def test_all_1(self):
        for ndim in range(2, 6):
            array = rng.random((1,) * ndim)
            result = _triu_cumsum(array)
            self.assertEqual(result.shape, array.shape)
            self.assertTrue(_triu_allclose(result, array))

    def test_1x1x2(self):
        x = rng.random()
        y = rng.random()
        array = np.array([[[x, y]]])
        self.assertTrue(_triu_allclose(_triu_cumsum(array), array))

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array([[3., 2.],
                           [10., 6.]])
        self.assertTrue(_triu_allclose(_triu_cumsum(array), expect))

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([[[3.], [2.]],
                           [[10.], [6.]]])
        self.assertTrue(_triu_allclose(_triu_cumsum(array), expect))

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[3., 11.], [2., 6.]],
                           [[0., 0.], [6., 14.]]])
        self.assertTrue(_triu_allclose(_triu_cumsum(array), expect))

    def test_3x3(self):
        array = np.array([[3., 4., 6.],
                          [7., 8., 9.],
                          [1., 2., 5.]])
        expect = np.array([[13., 10., 6.],
                           [0., 27., 15.],
                           [0., 0., 20.]])
        self.assertTrue(_triu_allclose(_triu_cumsum(array), expect))

    def test_explicit_sum(self):
        for npos in range(8):
            array = rng.random((npos, npos))
            result = _triu_cumsum(array)
            for row in range(npos):
                for col in range(npos):
                    if row > col:
                        continue
                    self.assertTrue(np.isclose(result[row, col],
                                               array[:row + 1, col:].sum()))


class TestTriuNorm(ut.TestCase):

    def compare(self, result: np.ndarray, expect: np.ndarray):
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))

    def test_0x0(self):
        array = rng.random((0, 0))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_0x0x1(self):
        array = rng.random((0, 0, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1(self):
        array = rng.random((1, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1x1(self):
        array = rng.random((1, 1, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1x2(self):
        array = rng.random((1, 1, 2))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array([[1 / 7, 2 / 7],
                           [np.nan, 4 / 7]])
        self.compare(_triu_norm(array), expect)

    def test_2x2_zero(self):
        array = np.array([[0., 0.],
                          [3., 0.]])
        expect = np.array([[1 / 3, 1 / 3],
                           [np.nan, 1 / 3]])
        self.compare(_triu_norm(array), expect)

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([[[1 / 7], [2 / 7]],
                           [[np.nan], [4 / 7]]])
        self.compare(_triu_norm(array), expect)

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[1 / 7, 5 / 19], [2 / 7, 6 / 19]],
                           [[np.nan, np.nan], [4 / 7, 8 / 19]]])
        self.compare(_triu_norm(array), expect)

    def test_2x2x2_zero(self):
        array = np.array([[[1., 0.], [2., 0.]],
                          [[3., 5.], [4., 0.]]])
        expect = np.array([[[1 / 7, 1 / 3], [2 / 7, 1 / 3]],
                           [[np.nan, np.nan], [4 / 7, 1 / 3]]])
        self.compare(_triu_norm(array), expect)

    def test_2x2x2x2(self):
        array = np.array([[[[1., 5.],
                            [9., 13.]], [[2., 6.],
                                         [10., 14.]]],
                          [[[3., 7.],
                            [11., 15.]], [[4., 8.],
                                          [12., 16.]]]])
        expect = np.array([[[[1 / 7, 5 / 19],
                             [9 / 31, 13 / 43]], [[2 / 7, 6 / 19],
                                                  [10 / 31, 14 / 43]]],
                           [[[np.nan, np.nan],
                             [np.nan, np.nan]], [[4 / 7, 8 / 19],
                                                 [12 / 31, 16 / 43]]]])
        self.compare(_triu_norm(array), expect)

    def test_2x2x2x2_zero(self):
        array = np.array([[[[0., 5.],
                            [9., 0.]], [[0., 6.],
                                        [10., 0.]]],
                          [[[3., 7.],
                            [0., 15.]], [[0., 8.],
                                         [12., 0.]]]])
        expect = np.array([[[[1 / 3, 5 / 19],
                             [9 / 31, 1 / 3]], [[1 / 3, 6 / 19],
                                                [10 / 31, 1 / 3]]],
                           [[[np.nan, np.nan],
                             [np.nan, np.nan]], [[1 / 3, 8 / 19],
                                                 [12 / 31, 1 / 3]]]])
        self.compare(_triu_norm(array), expect)


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


class TestPrivateCalcPNoCloseGivenEnds(ut.TestCase):
    """ Test _calc_p_nomut_window and _calc_p_noclose_given_ends. """

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


class TestPublicCalcPNoCloseGivenEnds(ut.TestCase):

    def test_1_dim(self):
        max_n = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for gap in range(max_g + 1):
                p_noclose_given_ends = calc_p_noclose_given_ends(
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
                    p_noclose_given_ends = calc_p_noclose_given_ends(
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
                                           calc_p_noclose_given_ends,
                                           rng.random(dims),
                                           gap)


class TestCalcRectangularSum(ut.TestCase):

    @staticmethod
    def calc_spanning_sum_slow(array: np.ndarray):
        spanning_sum = np.empty(array.shape[1:])
        for j in range(array.shape[0]):
            spanning_sum[j] = np.sum(array[:(j + 1), j:], axis=(0, 1))
        return spanning_sum

    def test_2d(self):
        for n in range(5):
            with self.subTest(n=n):
                array = rng.random((n, n))
                fast_sum = _calc_rectangular_sum(array)
                slow_sum = self.calc_spanning_sum_slow(array)
                self.assertEqual(fast_sum.shape, (n,))
                self.assertEqual(slow_sum.shape, (n,))
                self.assertTrue(np.allclose(fast_sum, slow_sum))

    def test_3d(self):
        for n in range(5):
            for k in range(3):
                with self.subTest(n=n, k=k):
                    array = rng.random((n, n, k))
                    fast_sum = _calc_rectangular_sum(array)
                    slow_sum = self.calc_spanning_sum_slow(array)
                    self.assertEqual(fast_sum.shape, (n, k))
                    self.assertEqual(slow_sum.shape, (n, k))
                    self.assertTrue(np.allclose(fast_sum, slow_sum))


class TestCalcPMutGivenSpanNoClose(ut.TestCase):

    def test_simulated(self):
        from scipy.stats import binom
        confidence = 0.999
        # Simulate reads.
        n_pos = 6
        n_reads = 1_000_000
        p_mut_max = 0.5
        p_mut, p_ends, p_cls = simulate_params(n_pos, 1, p_mut_max)
        muts, info, end5s, end3s = simulate_reads(n_reads, p_mut[:, 0], p_ends)
        # Test each minimum gap between mutations (min_gap).
        for min_gap in [0, 3, n_pos]:
            with self.subTest(min_gap=min_gap):
                # Determine which reads have no two mutations too close.
                has_no_close = label_no_close_muts(muts, min_gap)
                # Calculate the theoretical probability of no mutations
                # in each window.
                p_nomut_window_theory = _calc_p_nomut_window(p_mut, min_gap)
                # Calculate the theoretical probability of no mutations
                # being too close for each pair of end coordinates.
                p_noclose_given_ends_theory = _calc_p_noclose_given_ends(
                    p_mut,
                    p_nomut_window_theory
                )
                # Compare the simulated probability of no mutations
                # being too close for each pair of end coordinates.
                for end5 in range(n_pos):
                    for end3 in range(end5, n_pos):
                        if end5 == end3:
                            # The read must have no close mutations.
                            self.assertEqual(
                                p_noclose_given_ends_theory[end5, end3, 0],
                                1.
                            )
                        else:
                            # Count the reads with those coordinates.
                            end53s = np.logical_and(end5s == end5,
                                                    end3s == end3)
                            # Find the fraction of simualted reads with
                            # those coordinates without close mutations.
                            p_noclose_given_ends_simulated = (
                                    np.sum(np.logical_and(end53s, has_no_close))
                                    / np.sum(end53s)
                            )
                            # Compute the confidence interval for the
                            # proportion of reads.
                            n_expect = round(n_reads * p_ends[end5, end3])
                            p_expect = p_noclose_given_ends_theory[
                                end5, end3, 0
                            ]
                            inter_lo, inter_up = binom.interval(confidence,
                                                                n_expect,
                                                                p_expect)
                            # Ensure the interval based on theoretical
                            # probability agrees with the simulation.
                            self.assertGreaterEqual(
                                p_noclose_given_ends_simulated,
                                inter_lo / n_expect
                            )
                            self.assertLessEqual(
                                p_noclose_given_ends_simulated,
                                inter_up / n_expect
                            )
                # Compute the theoretical proportion of reads aligning
                # to each pair of 5' and 3' coordinates.
                p_ends_given_noclose_theory = calc_p_ends_given_noclose(
                    p_ends,
                    p_noclose_given_ends_theory
                ).reshape((n_pos, n_pos))
                # Compare to the simulated proportion of reads aligning
                # to each pair of coordinates.
                p_noclose = calc_p_noclose(p_ends, p_noclose_given_ends_theory)
                n_expect = np.round(p_ends_given_noclose_theory
                                    * (p_noclose * n_reads))
                inter_lo, inter_up = binom.interval(confidence,
                                                    n_expect,
                                                    p_ends_given_noclose_theory)
                inter_lo = _triu_div(inter_lo, n_expect)
                inter_up = _triu_div(inter_up, n_expect)
                end5s_noclose = end5s[has_no_close]
                end3s_noclose = end3s[has_no_close]
                for end5 in range(n_pos):
                    for end3 in range(end5, n_pos):
                        p_ends_given_noclose_simulated = np.mean(np.logical_and(
                            end5s_noclose == end5,
                            end3s_noclose == end3
                        ))
                        self.assertGreaterEqual(p_ends_given_noclose_simulated,
                                                inter_lo[end5, end3])
                        self.assertLessEqual(p_ends_given_noclose_simulated,
                                             inter_up[end5, end3])
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
                    p_noclose_given_ends_theory,
                    p_nomut_window_theory
                ).reshape(n_pos)
                # Compare the simulated probabilities of mutations given
                # no mutations are too close.
                p_mut_given_noclose_simulated = (np.sum(muts[has_no_close], axis=0)
                                                 /
                                                 np.sum(info[has_no_close], axis=0))
                inter_lo, inter_up = binom.interval(confidence,
                                                    n_expect,
                                                    p_mut_given_noclose_theory)
                self.assertTrue(np.all(p_mut_given_noclose_simulated
                                       >= inter_lo / n_expect))
                self.assertTrue(np.all(p_mut_given_noclose_simulated
                                       <= inter_up / n_expect))

    def test_clusters(self):
        n_pos = 16
        max_clusters = 3
        max_gap = 3
        # Test each number of clusters.
        for n_cls in range(1, max_clusters + 1):
            p_mut, p_ends, p_cls = simulate_params(n_pos, n_cls)
            # Test each minimum gap between mutations.
            for min_gap in range(max_gap + 1):
                with self.subTest(n_cls=n_cls, min_gap=min_gap):
                    # Compute the mutation rates with no two mutations
                    # too close for all clusters simultaneously.
                    p_nomut_window = _calc_p_nomut_window(p_mut, min_gap)
                    p_noclose_given_ends = _calc_p_noclose_given_ends(
                        p_mut, p_nomut_window
                    )
                    p_mut_given_span_noclose = _calc_p_mut_given_span_noclose(
                        p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                    )
                    # Compare to computing each cluster individually.
                    for k in range(n_cls):
                        self.assertTrue(np.allclose(
                            p_mut_given_span_noclose[:, k],
                            _calc_p_mut_given_span_noclose(
                                p_mut[:, [k]],
                                p_ends,
                                p_noclose_given_ends[:, :, [k]],
                                p_nomut_window[:, :, [k]],
                            ).reshape(n_pos)
                        ))


class TestSlicePEnds(ut.TestCase):

    def test_0x0(self):
        p_ends = np.empty((0, 0))
        self.assertTrue(np.array_equal(_slice_p_ends(p_ends, p_ends, 0, 0),
                                       p_ends))

    def test_slice_3x3(self):
        p_ends = np.array([[1., 2., 3.],
                           [4., 5., 6.],
                           [7., 8., 9.]])
        p_ends_cumsum = _triu_cumsum(p_ends)
        result = _slice_p_ends(p_ends, p_ends_cumsum, 1, 2)
        expect = np.array([[16.]])
        self.assertTrue(_triu_allclose(result, expect))

    def test_slice_5x5(self):
        p_ends = np.arange(25.).reshape((5, 5))
        p_ends_cumsum = _triu_cumsum(p_ends)
        result = _slice_p_ends(p_ends, p_ends_cumsum, 1, 4)
        expect = np.array([[7., 9., 24.],
                           [0., 12., 27.],
                           [0., 0., 37.]])
        self.assertTrue(_triu_allclose(result, expect))


class TestFindSplitPositions(ut.TestCase):

    def test_0(self):
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(np.array([[]]), min_gap, 0.),
                np.array([], dtype=int)
            ))

    def test_thresh0(self):
        p_mut = 1. - rng.random((10, 2))
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(p_mut, min_gap, 0.),
                np.array([], dtype=int)
            ))

    def test_thresh1(self):
        p_mut = 1. - rng.random((10, 2))
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(p_mut, min_gap, 1.),
                np.array([], dtype=int)
            ))

    def test_gap0(self):
        p_mut = 1. - rng.random((10, 2))
        for thresh in np.linspace(0., 1., 5):
            self.assertTrue(np.array_equal(
                _find_split_positions(p_mut, 0, thresh),
                np.array([], dtype=int)
            ))

    def test_gap1_split1_single_mid(self):
        p_mut = np.array([[0.2, 0.1, 0.1, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([3, 4])
        ))

    def test_gap1_single_end5(self):
        p_mut = np.array([[0.0, 0.1, 0.1, 0.3, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([1], dtype=int)
        ))

    def test_gap2_single_end5(self):
        p_mut = np.array([[0.0, 0.0, 0.1, 0.3, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 2, 0.),
            np.array([2], dtype=int)
        ))

    def test_gap3_single_end5(self):
        p_mut = np.array([[0.0, 0.0, 0.0, 0.3, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 3, 0.),
            np.array([3], dtype=int)
        ))

    def test_gap1_single_end3(self):
        p_mut = np.array([[0.2, 0.1, 0.1, 0.3, 0.1, 0.0]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([5], dtype=int)
        ))

    def test_gap2_single_end3(self):
        p_mut = np.array([[0.2, 0.1, 0.1, 0.3, 0.0, 0.0]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 2, 0.),
            np.array([4], dtype=int)
        ))

    def test_gap3_single_end3(self):
        p_mut = np.array([[0.2, 0.1, 0.1, 0.0, 0.0, 0.0]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 3, 0.),
            np.array([3], dtype=int)
        ))

    def test_gap1_split1_double(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([2, 4])
        ))

    def test_gap1_split1_triple(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.0, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([2, 5])
        ))

    def test_gap1_split0_quadruple(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.0, 0.0]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([2], dtype=int)
        ))

    def test_gap1_split2(self):
        p_mut = np.array([[0.2, 0.0, 0.1, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 1, 0.),
            np.array([1, 2, 3, 4])
        ))

    def test_gap2_split0(self):
        p_mut = np.array([[0.2, 0.0, 0.1, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 2, 0.),
            np.array([], dtype=int)
        ))

    def test_gap2_split1(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 2, 0.),
            np.array([2, 4])
        ))

    def test_gap4_split0(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 4, 0.),
            np.array([], dtype=int)
        ))

    def test_gap4_split1(self):
        p_mut = np.array([[0.2, 0.1, 0.0, 0.0, 0.1, 0.2]]).T
        self.assertTrue(np.array_equal(
            _find_split_positions(p_mut, 4, 0.1),
            np.array([1, 5])
        ))

    def test_generic_split(self):

        def find_split_positions_slow(p_mut: np.ndarray,
                                      min_gap: int,
                                      threshold: int):
            p_mut = p_mut.max(axis=1)
            n_pos = p_mut.size
            splits = list()
            if n_pos > 0 and min_gap > 0:
                current_below_count = 0
                for pos in range(n_pos):
                    below_thresh = p_mut[pos] <= threshold
                    if below_thresh:
                        if current_below_count == 0 and pos > 0:
                            # Begin a below-threshold stretch.
                            splits.append(pos)
                        current_below_count += 1
                    else:
                        if current_below_count > 0:
                            if current_below_count >= min_gap:
                                # End a sufficiently long stretch.
                                splits.append(pos)
                            elif splits:
                                # End an insufficiently long stretch.
                                splits.pop()
                            current_below_count = 0
                if 0 < current_below_count < min_gap and splits:
                    # End an insufficiently long stretch.
                    splits.pop()
            return np.array(splits, dtype=int)

        p = rng.random((100, 2))
        for gap in range(4):
            for thresh in np.linspace(0., 1., 6):
                self.assertTrue(np.array_equal(
                    _find_split_positions(p, gap, thresh),
                    find_split_positions_slow(p, gap, thresh)
                ))


class TestQuickUnbias(ut.TestCase):
    """ Test that the quick unbiasing algorithm produces results very
    similar to the exact algorithm. """

    ABS_TOL = 5.e-4
    REL_TOL = 3.e-2

    @staticmethod
    def random_params(npos: int, ncls: int):
        p_mut = rng.beta(0.5, 9.5, (npos, ncls))
        p_ends = 1. - rng.random((npos, npos, ncls))
        p_cls = 1. - rng.random(ncls)
        return p_mut, p_ends, p_cls

    def test_threshold_0(self):
        t = 0.
        n_pos = 500
        for n_cls in [1, 2]:
            for min_gap in [1, 3]:
                for f_below in [0.2, 0.5, 0.8]:
                    p_mut, p_ends, p_cls = self.random_params(n_pos, n_cls)
                    # Set some mutation rates to below the threshold.
                    n_below = round(n_pos * f_below)
                    below_t = rng.choice(n_pos, n_below, replace=False)
                    is_below_t = np.zeros(n_pos, dtype=bool)
                    is_below_t[below_t] = True
                    p_mut[is_below_t] = t
                    # Compare quick vs. exact unbias.
                    (p_mut_quick,
                     p_ends_quick,
                     p_clust_quick) = calc_params(p_mut, p_ends, p_cls, min_gap,
                                                  quick_unbias=True,
                                                  quick_unbias_thresh=t)
                    (p_mut_exact,
                     p_ends_exact,
                     p_clust_exact) = calc_params(p_mut, p_ends, p_cls, min_gap,
                                                  quick_unbias=False,
                                                  quick_unbias_thresh=t)
                    with self.subTest(n_cls=n_cls,
                                      min_gap=min_gap,
                                      f_below=f_below):
                        self.assertTrue(np.allclose(p_mut_quick[is_below_t],
                                                    0.))
                        self.assertTrue(np.allclose(p_mut_quick[~is_below_t],
                                                    p_mut_exact[~is_below_t],
                                                    atol=self.ABS_TOL,
                                                    rtol=self.REL_TOL))
                        self.assertTrue(np.allclose(np.triu(p_ends_quick),
                                                    np.triu(p_ends_exact),
                                                    atol=self.ABS_TOL,
                                                    rtol=self.REL_TOL))
                        self.assertTrue(np.allclose(p_clust_quick,
                                                    p_clust_exact,
                                                    atol=self.ABS_TOL,
                                                    rtol=self.REL_TOL))

    def test_threshold_0p001(self):
        t = 0.001
        n_pos = 500
        for n_cls in [1, 2]:
            for min_gap in [1, 3]:
                p_mut, p_ends, p_cls = self.random_params(n_pos, n_cls)
                # Compare quick vs. exact unbias.
                (p_mut_quick,
                 p_ends_quick,
                 p_clust_quick) = calc_params(p_mut, p_ends, p_cls, min_gap,
                                              quick_unbias=True,
                                              quick_unbias_thresh=t)
                (p_mut_exact,
                 p_ends_exact,
                 p_clust_exact) = calc_params(p_mut, p_ends, p_cls, min_gap,
                                              quick_unbias=False,
                                              quick_unbias_thresh=t)
                with self.subTest(n_cls=n_cls, min_gap=min_gap):
                    self.assertTrue(np.allclose(p_mut_quick,
                                                p_mut_exact,
                                                atol=self.ABS_TOL,
                                                rtol=self.REL_TOL))
                    self.assertTrue(np.allclose(np.triu(p_ends_quick),
                                                np.triu(p_ends_exact),
                                                atol=self.ABS_TOL,
                                                rtol=self.REL_TOL))
                    self.assertTrue(np.allclose(p_clust_quick,
                                                p_clust_exact,
                                                atol=self.ABS_TOL,
                                                rtol=self.REL_TOL))


class TestCalcPMutGivenSpan(ut.TestCase):

    def test_invert(self):
        """ Test the inverse of `_calc_p_mut_given_span_noclose`. """
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, p_cls = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the mutation rates without mutations too close.
                p_nomut_window = _calc_p_nomut_window(p_mut, min_gap)
                p_noclose_given_ends = _calc_p_noclose_given_ends(
                    p_mut, p_nomut_window
                )
                p_mut_given_span_noclose = _calc_p_mut_given_span_noclose(
                    p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                )
                p_mut_given_span = _calc_p_mut_given_span(
                    p_mut_given_span_noclose,
                    min_gap,
                    p_ends,
                    p_mut_given_span_noclose,
                )
                self.assertEqual(p_mut_given_span.shape, p_mut.shape)
                self.assertTrue(np.allclose(p_mut_given_span,
                                            p_mut,
                                            atol=1.e-4,
                                            rtol=1.e-3))


class TestCalcPEnds(ut.TestCase):

    def test_invert(self):
        """ Test the inverse of `calc_p_ends_given_noclose`. """
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, p_cls = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the coordinate distributions without mutations
                # too close.
                p_noclose_given_ends = _calc_p_noclose_given_ends(
                    p_mut, _calc_p_nomut_window(p_mut, min_gap)
                )
                p_ends_given_noclose = calc_p_ends_given_noclose(
                    p_ends, p_noclose_given_ends
                )
                # Infer the original distribution of all reads.
                p_ends_inferred = _calc_p_ends(p_ends_given_noclose,
                                               p_noclose_given_ends,
                                               p_mut,
                                               p_cls)
                self.assertEqual(p_ends_inferred.shape, p_ends.shape)
                self.assertTrue(_triu_allclose(p_ends_inferred, p_ends))


class TestCalcPNoClose(ut.TestCase):

    def test_p_noclose(self):
        p_ends = np.array([[0.2, 0.5],
                           [0.4, 0.3]])
        p_noclose_given_ends = np.array([[[0.2, 0.4], [0.1, 0.2]],
                                         [[0.3, 0.6], [0.4, 0.8]]])
        expect = np.array([0.21, 0.42])
        result = calc_p_noclose(p_ends, p_noclose_given_ends)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcPClust(ut.TestCase):

    def test_p_clust(self):
        p_clust_given_observed = np.array([0.2, 0.5, 0.3])
        p_noclose = np.array([0.8, 0.4, 0.5])
        expect = np.array([0.25, 1.25, 0.60]) / 2.1
        result = calc_p_clust(p_clust_given_observed, p_noclose)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcPEndsObserved(ut.TestCase):

    def test_calc_p_ends_given_observed(self):
        npos = 3
        end5s = np.array([2, 0, 0, 1, 1, 2, 0, 0])
        end3s = np.array([2, 1, 2, 2, 1, 2, 0, 1])
        weights = np.array([[0.1, 0.2],
                            [0.2, 0.3],
                            [0.3, 0.4],
                            [0.4, 0.5],
                            [0.5, 0.6],
                            [0.6, 0.7],
                            [0.7, 0.8],
                            [0.8, 0.9]])
        expect = np.array([[[0.7, 0.8], [1.0, 1.2], [0.3, 0.4]],
                           [[0.0, 0.0], [0.5, 0.6], [0.4, 0.5]],
                           [[0.0, 0.0], [0.0, 0.0], [0.7, 0.9]]])
        result = _calc_p_ends_observed(npos, end5s, end3s, weights)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcParams(ut.TestCase):

    def test_infer(self):
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, p_cls = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the mutation rates and distributions of end
                # coordinates without mutations too close.
                p_nomut_window = _calc_p_nomut_window(p_mut, min_gap)
                p_noclose_given_ends = _calc_p_noclose_given_ends(
                    p_mut, p_nomut_window
                )
                p_noclose = calc_p_noclose(p_ends, p_noclose_given_ends)
                p_mut_given_span_noclose = _calc_p_mut_given_span_noclose(
                    p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                )
                p_ends_given_noclose = calc_p_ends_given_noclose(
                    p_ends, p_noclose_given_ends
                )
                p_clust_given_noclose = calc_p_clust_given_noclose(
                    p_cls, p_noclose
                )
                # Infer the original parameters using those of the reads
                # without mutations too close.
                (p_mut_inferred,
                 p_ends_inferred,
                 p_cls_inferred) = calc_params(
                    p_mut_given_span_noclose,
                    p_ends_given_noclose,
                    p_clust_given_noclose,
                    min_gap,
                )
                self.assertEqual(p_mut_inferred.shape, p_mut.shape)
                self.assertTrue(np.allclose(p_mut_inferred,
                                            p_mut,
                                            atol=1.e-4,
                                            rtol=1.e-2))
                self.assertEqual(p_ends_inferred.shape, p_ends.shape)
                self.assertTrue(_triu_allclose(p_ends_inferred, p_ends))
                self.assertEqual(p_cls_inferred.shape, p_cls.shape)
                self.assertTrue(np.allclose(p_cls_inferred,
                                            p_cls,
                                            atol=1.e-3,
                                            rtol=1.e-2))


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
