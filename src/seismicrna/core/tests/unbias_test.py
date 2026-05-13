import unittest as ut
from itertools import product

import numpy as np
from numba import jit

from seismicrna.core.unbias import (_clip,
                                    _normalize,
                                    _adjust_min_gap,
                                    _triu_cumsum,
                                    triu_dot,
                                    _triu_norm,
                                    _triu_mul,
                                    _triu_sum,
                                    _triu_div,
                                    calc_p_ends_observed,
                                    calc_p_nomut_window,
                                    calc_p_noclose_given_ends,
                                    calc_p_noclose_given_ends_auto,
                                    calc_p_mut_given_span_dropped,
                                    calc_p_mut_given_span,
                                    calc_p_ends,
                                    _slice_p_ends,
                                    _find_split_positions,
                                    calc_p_clust,
                                    calc_p_noclose_given_clust,
                                    calc_p_noclose,
                                    calc_p_ends_given_clust_noclose,
                                    calc_p_ends_given_noclose,
                                    calc_p_clust_given_ends_noclose,
                                    calc_p_clust_given_noclose,
                                    calc_params,
                                    calc_rectangular_sum,
                                    _calc_rectangular_sum_weighted,
                                    require_square_atleast2d,
                                    require_same_square_atleast2d,
                                    _calc_p_mut_given_span_merged)
from seismicrna.core.validate import require_equal



def triu_sum(a: np.ndarray):
    """ Calculate the sum over the upper triangle(s) of array `a`.

    Parameters
    ----------
    a: np.ndarray
        Array whose upper triangle to sum.

    Returns
    -------
    np.ndarray
        Sum of the upper triangle(s), with the same shape as the third
        and subsequent dimensions of `a`.
    """
    require_square_atleast2d("a", a)
    return _triu_sum(a)


@jit()
def _triu_allclose(a: np.ndarray,
                   b: np.ndarray,
                   rtol: float = 1.e-3,
                   atol: float = 1.e-6):
    """ Whether the upper triangles of `a` and `b` are all close.

    This function is meant to be called by another function that has
    validated the arguments; hence, this function makes assumptions:

    -   `a` and `b` both have at least 2 dimensions.
    -   The first and second dimensions of `a` are equal.
    -   The first and second dimensions of `b` are equal.

    Parameters
    ----------
    a: np.ndarray
        Array 1.
    b: np.ndarray
        Array 2.
    rtol: float = 1.0e-3
        Relative tolerance.
    atol: float = 1.0e-6
        Absolute tolerance.

    Returns
    -------
    bool
        Whether all elements of the upper triangles of `a` and `b` are
        close using the function `np.allclose`.
    """
    for j in range(a.shape[0]):
        if not np.allclose(a[j, j:], b[j, j:], rtol=rtol, atol=atol):
            return False
    return True


def triu_allclose(a: np.ndarray | float,
                  b: np.ndarray | float,
                  rtol: float = 1.e-3,
                  atol: float = 1.e-6):
    """ Whether the upper triangles of `a` and `b` are all close.

    Parameters
    ----------
    a: np.ndarray | float
        Array 1.
    b: np.ndarray | float
        Array 2.
    rtol: float = 1.0e-3
        Relative tolerance.
    atol: float = 1.0e-6
        Absolute tolerance.

    Returns
    -------
    bool
        Whether all elements of the upper triangles of `a` and `b` are
        close using the function `np.allclose`.
    """
    # Ensure a and b are arrays with the same dimensions.
    a, b = np.broadcast_arrays(np.asarray(a), np.asarray(b))
    # Explicitly set the writeable flag to False in order to suppress
    # "FutureWarning: future versions will not create a writeable array
    # from broadcast_array. Set the writable flag explicitly to avoid
    # this warning."
    # These lines may not be necessary with a future version of NumPy.
    a.flags.writeable = False
    b.flags.writeable = False
    # Ensure a and b are both square in their first 2 dimensions.
    require_same_square_atleast2d(a, b)
    return _triu_allclose(a, b, rtol, atol)


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


def merge_mutations_right_to_left(muts: np.ndarray, min_gap: int):
    """ Merge close mutations row by row from right to left.

    For each row, scan from the last column to the first. When a True
    value is encountered, set the next `min_gap` positions to the left
    to False, then continue scanning.
    """
    if not isinstance(muts, np.ndarray) or muts.ndim != 2:
        raise ValueError("muts must be a 2D NumPy array")
    if muts.dtype != np.bool_:
        raise ValueError("muts must be a boolean NumPy array")
    if not isinstance(min_gap, int) or min_gap < 0:
        raise ValueError("min_gap must be a non-negative integer")
    merged = muts.copy()
    n_reads, n_pos = merged.shape
    # If min_gap is 0 or n_pos ≤ 1, then there is nothing to do.
    # A performance optimized function would return here, but this
    # function is meant for testing and so we want to test that the
    # algorithm produces the correct result in this case as well.
    for i in range(n_reads):
        j = n_pos - 1
        while j >= 0:
            if merged[i, j]:
                # Set the next min_gap positions to the left to False.
                left = max(0, j - min_gap)
                merged[i, left: j] = False
                j = left - 1
            else:
                j -= 1
    return merged


def simulate_reads(n_reads: int, p_mut: np.ndarray, p_ends: np.ndarray):
    """ Simulate `n_reads` reads based on the mutation rates (`p_mut`)
    and the distributions of end coordinates (`p_ends`). """
    rng = np.random.default_rng(seed=0)
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
    """ Return `p_mut`, `p_ends`, and `p_clust` parameters. """
    rng = np.random.default_rng(seed=0)
    p_mut = p_mut_max * rng.random((n_pos, n_cls))
    p_ends = np.triu(1. - rng.random((n_pos, n_pos)))
    p_ends /= np.sum(p_ends)
    p_clust = 1. - rng.random(n_cls)
    p_clust /= np.sum(p_clust)
    return p_mut, p_ends, p_clust


class TestCalcPEndsGivenClustNoClose(ut.TestCase):

    def test_npos0_ncls1(self):
        p_ends = np.ones((0, 0))
        p_noclose_given_ends = np.ones((0, 0, 1))
        self.assertRaisesRegex(ValueError,
                               r"Must have size\(positions\) ≥ 1, but got 0",
                               calc_p_ends_given_clust_noclose,
                               p_ends,
                               p_noclose_given_ends)

    def test_npos1_ncls0(self):
        p_ends = np.ones((1, 1))
        p_noclose_given_ends = np.ones((1, 1, 0))
        self.assertRaisesRegex(ValueError,
                               r"Must have size\(clusters\) ≥ 1, but got 0",
                               calc_p_ends_given_clust_noclose,
                               p_ends,
                               p_noclose_given_ends)

    def test_npos1_ncls1(self):
        p_ends = np.ones((1, 1))
        expect = np.ones((1, 1, 1))
        for p in [0.01, 0.1, 1.0]:
            p_noclose_given_ends = np.full((1, 1, 1), p)
            result = calc_p_ends_given_clust_noclose(p_ends,
                                                     p_noclose_given_ends)
            self.assertTrue(np.array_equal(result, expect))

    def test_npos2_ncls1(self):
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_noclose_given_ends = np.array([[[0.9], [0.6]],
                                         [[0.0], [0.8]]])
        expect = np.array([[[3. / 12.], [5. / 12.]],
                           [[np.nan], [4. / 12.]]])
        result = calc_p_ends_given_clust_noclose(p_ends, p_noclose_given_ends)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))


class TestCalcPEndsGivenNoClose(ut.TestCase):

    def test_calc_p_ends_given_noclose(self):
        p_ends_given_clust_noclose = np.array(
            [[[0.2, 0.1, 0.5], [0.5, 0.7, 0.3]],
             [[0.0, 0.0, 0.0], [0.3, 0.2, 0.2]]]
        )
        p_clust_given_noclose = np.array([0.5, 0.3, 0.2])
        expect = np.array([[0.23, 0.52],
                           [np.nan, 0.25]])
        result = calc_p_ends_given_noclose(p_ends_given_clust_noclose,
                                           p_clust_given_noclose)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))


class TestCalcPClustGivenEndsNoClose(ut.TestCase):

    def test_calc_p_clust_given_ends_noclose(self):
        p_ends_given_clust_noclose = np.array(
            [[[0.2, 0.1, 0.5], [0.5, 0.7, 0.3]],
             [[0.0, 0.0, 0.0], [0.3, 0.2, 0.2]]]
        )
        p_clust_given_noclose = np.array([0.5, 0.3, 0.2])
        expect = np.array([[[0.43478261, 0.13043478, 0.43478261],
                            [0.48076923, 0.40384615, 0.11538462]],
                           [[np.nan, np.nan, np.nan],
                            [0.6, 0.24, 0.16]]])
        result = calc_p_clust_given_ends_noclose(p_ends_given_clust_noclose,
                                                 p_clust_given_noclose)
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
                p_noclose_given_clust = calc_p_noclose_given_clust(
                    p_ends, p_noclose_given_ends
                )
                result = calc_p_clust_given_noclose(p_clust,
                                                    p_noclose_given_clust)
                self.assertEqual(result.shape, expect.shape)
                self.assertTrue(np.allclose(result, expect))

    def test_ncls2(self):
        p_clust = np.ones((1,))
        expect = np.ones((1,))
        for npos in range(1, 5):
            with self.subTest(npos=npos):
                p_ends = np.ones((npos, npos))
                p_noclose_given_ends = np.ones((npos, npos, 1))
                p_noclose_given_clust = calc_p_noclose_given_clust(
                    p_ends, p_noclose_given_ends
                )
                result = calc_p_clust_given_noclose(p_clust,
                                                    p_noclose_given_clust)
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


class TestMergeMutationsRightToLeft(ut.TestCase):

    @staticmethod
    def bits_to_matrix(bits: str):
        return np.array(
            [[bit == "1" for bit in row] for row in bits.split()],
            dtype=bool
        )
    
    def test_bits_to_matrix(self):
        bits = "1001011011 0001110101"
        expect = np.array(
            [[True, False, False, True, False, True, True, False, True, True],
             [False, False, False, True, True, True, False, True, False, True]],
            dtype=bool
        )
        result = self.bits_to_matrix(bits)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.array_equal(result, expect))

    def test_zero_columns(self):
        for n_reads in range(5):
            with self.subTest(n_reads=n_reads):
                muts = np.zeros((n_reads, 0), dtype=bool)
                result = merge_mutations_right_to_left(muts, 3)
                self.assertEqual(result.shape, muts.shape)
                self.assertTrue(np.array_equal(result, muts))

    def test_one_column(self):
        muts = np.array([[False],
                         [True],
                         [False],
                         [True]],
                        dtype=bool)
        for min_gap in range(5):
            with self.subTest(min_gap=min_gap):
                result = merge_mutations_right_to_left(muts, min_gap)
                self.assertEqual(result.shape, muts.shape)
                self.assertTrue(np.array_equal(result, muts))
    
    def test_min_gap_0(self):
        rng = np.random.default_rng(seed=0)
        for n_reads in range(4):
            for n_pos in range(6):
                with self.subTest(n_reads=n_reads, n_pos=n_pos):
                    muts = rng.integers(0, 2, size=(n_reads, n_pos), dtype=bool)
                    result = merge_mutations_right_to_left(muts, 0)
                    self.assertEqual(result.shape, muts.shape)
                    self.assertTrue(np.array_equal(result, muts))
                    self.assertFalse(np.shares_memory(result, muts))
    
    def test_min_gaps_0_to_6(self):
        muts = self.bits_to_matrix("10010110110101")
        expect = {0: "10010110110101",
                  1: "10010010010101",
                  2: "10010010010001",
                  3: "10000100010001",
                  4: "00010000100001",
                  5: "10000010000001",
                  6: "00000010000001"}
        for min_gap, bits in expect.items():
            with self.subTest(min_gap=min_gap):
                self.assertTrue(np.array_equal(
                    merge_mutations_right_to_left(muts, min_gap),
                    self.bits_to_matrix(bits)
                ))


class TestClip(ut.TestCase):

    def test_with_clip(self):
        rng = np.random.default_rng(seed=0)
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
        rng = np.random.default_rng(seed=0)
        mus = rng.random(64, dtype=float)
        # All values are in [0, 1] and so should not be clipped.
        self.assertTrue(np.allclose(mus, _clip(mus)))


class TestNormalize(ut.TestCase):

    def test_sum_positive(self):
        rng = np.random.default_rng(seed=0)
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


class TestTriuAllClose(ut.TestCase):

    def test_equal(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(4):
            for extra in range(3):
                for extras in product(range(4), repeat=extra):
                    dims = (npos, npos) + extras
                    a = rng.random(dims)
                    self.assertTrue(triu_allclose(a, a.copy()))

    def test_triu(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(2, 4):
            a = rng.random((npos, npos))
            self.assertTrue(triu_allclose(a, np.triu(a)))

    def test_tril(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(2, 4):
            a = rng.random((npos, npos))
            self.assertFalse(triu_allclose(a, np.tril(a)))

    def test_float_equal(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(4):
            for extra in range(3):
                for extras in product(range(4), repeat=extra):
                    dims = (npos, npos) + extras
                    b = rng.random()
                    a = np.full(dims, b)
                    self.assertTrue(triu_allclose(a, b))
                    self.assertTrue(triu_allclose(b, a))

    def test_float_unequal(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(1, 4):
            for extra in range(3):
                for extras in product(range(1, 4), repeat=extra):
                    dims = (npos, npos) + extras
                    b = rng.random()
                    a = np.full(dims, b + 1.)
                    self.assertFalse(triu_allclose(a, b))
                    self.assertFalse(triu_allclose(b, a))


class TestTriuDot(ut.TestCase):

    def test_1x1(self):
        a = np.array([[2.]])
        b = np.array([[4.]])
        expect = np.array(8.)
        self.assertTrue(np.array_equal(triu_dot(a, b), expect))

    def test_1x1x1(self):
        a = np.array([[[2.]]])
        b = np.array([[[4.]]])
        expect = np.array([8.])
        self.assertTrue(np.array_equal(triu_dot(a, b), expect))

    def test_2x2(self):
        a = np.array([[2., 3.],
                      [5., 7.]])
        b = np.array([[4., 8.],
                      [16., 32.]])
        expect = np.array(256.)
        self.assertTrue(np.array_equal(triu_dot(a, b), expect))

    def test_2x2x2(self):
        a = np.array([[[2., 20.], [3., 30.]],
                      [[5., 50.], [7., 70.]]])
        b = np.array([[[4., 40.], [8., 80.]],
                      [[16., 160.], [32., 320.]]])
        expect = np.array([256., 25600.])
        self.assertTrue(np.array_equal(triu_dot(a, b), expect))


class TestTriuMul(ut.TestCase):

    def test_1x1(self):
        rng = np.random.default_rng(seed=0)
        a = rng.random()
        b = 1. - rng.random()
        expect = np.array([[a * b]])
        self.assertTrue(np.array_equal(_triu_mul(np.array([[a]]),
                                                 np.array([[b]])),
                                       expect))

    def test_1x1x1(self):
        rng = np.random.default_rng(seed=0)
        a = rng.random()
        b = 1. - rng.random()
        expect = np.array([[[a * b]]])
        self.assertTrue(np.array_equal(_triu_mul(np.array([[[a]]]),
                                                 np.array([[[b]]])),
                                       expect))

    def test_2x2(self):
        a = np.array([[12., 3.],
                      [20., 56.]])
        b = np.array([[2., 3.],
                      [5., 7.]])
        expect = np.array([[24., 9.],
                           [np.nan, 392.]])
        self.assertTrue(triu_allclose(_triu_mul(a, b), expect))

    def test_2x2x2(self):
        a = np.array([[[2., 20.], [2., 36.]],
                      [[3., 7.], [12., 40.]]])
        b = np.array([[[1., 5.], [2., 6.]],
                      [[3., 7.], [4., 8.]]])
        expect = np.array([[[2., 100.], [4., 216.]],
                           [[np.nan, np.nan], [48., 320.]]])
        self.assertTrue(triu_allclose(_triu_mul(a, b), expect))


class TestTriuDiv(ut.TestCase):

    def test_1x1(self):
        rng = np.random.default_rng(seed=0)
        n = rng.random()
        d = 1. - rng.random()
        expect = np.array([[n / d]])
        self.assertTrue(np.array_equal(_triu_div(np.array([[n]]),
                                                 np.array([[d]])),
                                       expect))

    def test_1x1x1(self):
        rng = np.random.default_rng(seed=0)
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
        self.assertTrue(triu_allclose(_triu_div(numer, denom), expect))

    def test_2x2x2(self):
        numer = np.array([[[2., 20.], [2., 36.]],
                          [[3., 7.], [12., 40.]]])
        denom = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[2., 4.], [1., 6.]],
                           [[np.nan, np.nan], [3., 5.]]])
        self.assertTrue(triu_allclose(_triu_div(numer, denom), expect))


class TestTriuSum(ut.TestCase):

    def test_all_zero(self):
        rng = np.random.default_rng(seed=0)
        for ndim in range(2, 6):
            array = rng.random((0,) * ndim)
            expect = np.zeros((0,) * (ndim - 2))
            self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_1x1(self):
        rng = np.random.default_rng(seed=0)
        x = rng.random()
        array = np.array([[x]])
        expect = np.array(x)
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_1x1x1(self):
        rng = np.random.default_rng(seed=0)
        x = rng.random()
        array = np.array([[[x]]])
        expect = np.array([x])
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_1x1x2(self):
        rng = np.random.default_rng(seed=0)
        x = rng.random()
        y = rng.random()
        array = np.array([[[x, y]]])
        expect = np.array([x, y])
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array(7.)
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([7.])
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([7., 19.])
        self.assertTrue(np.array_equal(triu_sum(array), expect))

    def test_2x2x2x2(self):
        array = np.array([[[[1., 5.],
                            [9., 13.]], [[2., 6.],
                                         [10., 14.]]],
                          [[[3., 7.],
                            [11., 15.]], [[4., 8.],
                                          [12., 16.]]]])
        expect = np.array([[7., 19.],
                           [31., 43.]])
        self.assertTrue(np.array_equal(triu_sum(array), expect))


class TestTriuCumSum(ut.TestCase):

    def test_all_0(self):
        rng = np.random.default_rng(seed=0)
        for ndim in range(2, 6):
            array = rng.random((0,) * ndim)
            result = _triu_cumsum(array)
            self.assertEqual(result.shape, array.shape)
            self.assertTrue(triu_allclose(result, array))

    def test_all_1(self):
        rng = np.random.default_rng(seed=0)
        for ndim in range(2, 6):
            array = rng.random((1,) * ndim)
            result = _triu_cumsum(array)
            self.assertEqual(result.shape, array.shape)
            self.assertTrue(triu_allclose(result, array))

    def test_1x1x2(self):
        rng = np.random.default_rng(seed=0)
        x = rng.random()
        y = rng.random()
        array = np.array([[[x, y]]])
        self.assertTrue(triu_allclose(_triu_cumsum(array), array))

    def test_2x2(self):
        array = np.array([[1., 2.],
                          [3., 4.]])
        expect = np.array([[3., 2.],
                           [10., 6.]])
        self.assertTrue(triu_allclose(_triu_cumsum(array), expect))

    def test_2x2x1(self):
        array = np.array([[[1.], [2.]],
                          [[3.], [4.]]])
        expect = np.array([[[3.], [2.]],
                           [[10.], [6.]]])
        self.assertTrue(triu_allclose(_triu_cumsum(array), expect))

    def test_2x2x2(self):
        array = np.array([[[1., 5.], [2., 6.]],
                          [[3., 7.], [4., 8.]]])
        expect = np.array([[[3., 11.], [2., 6.]],
                           [[0., 0.], [6., 14.]]])
        self.assertTrue(triu_allclose(_triu_cumsum(array), expect))

    def test_3x3(self):
        array = np.array([[3., 4., 6.],
                          [7., 8., 9.],
                          [1., 2., 5.]])
        expect = np.array([[13., 10., 6.],
                           [0., 27., 15.],
                           [0., 0., 20.]])
        self.assertTrue(triu_allclose(_triu_cumsum(array), expect))

    def test_explicit_sum(self):
        rng = np.random.default_rng(seed=0)
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
        rng = np.random.default_rng(seed=0)
        array = rng.random((0, 0))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_0x0x1(self):
        rng = np.random.default_rng(seed=0)
        array = rng.random((0, 0, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1(self):
        rng = np.random.default_rng(seed=0)
        array = rng.random((1, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1x1(self):
        rng = np.random.default_rng(seed=0)
        array = rng.random((1, 1, 1))
        expect = np.ones_like(array)
        self.compare(_triu_norm(array), expect)

    def test_1x1x2(self):
        rng = np.random.default_rng(seed=0)
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


class TestCalcPNoCloseGivenEnds(ut.TestCase):
    """ Test calc_p_nomut_window and calc_p_noclose_given_ends. """

    def test_min_gap_0(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(10):
            with self.subTest(npos=npos):
                mu = rng.random((npos, 1))
                nm = calc_p_nomut_window(mu, 0)
                nm_expect = np.ones((1, npos + 1, 1))
                self.assertTrue(np.array_equal(nm, nm_expect))
                nc = calc_p_noclose_given_ends(mu, nm)
                nc_expect = np.ones((npos, npos, 1))
                self.assertTrue(np.array_equal(nc, nc_expect))

    def test_length_0(self):
        mu = np.zeros((0, 1))
        for min_gap in range(4):
            with self.subTest(min_gap=min_gap):
                nm = calc_p_nomut_window(mu, min_gap)
                nm_expect = np.ones((1, 1, 1))
                self.assertTrue(np.array_equal(nm, nm_expect))
                nc = calc_p_noclose_given_ends(mu, nm)
                nc_expect = np.ones((0, 0, 1))
                self.assertTrue(np.array_equal(nc, nc_expect))

    def test_length_1(self):
        for x in [0., 0.01, 0.1, 1.]:
            mu = np.array([[x]])
            for min_gap in range(4):
                with self.subTest(x=x, min_gap=min_gap):
                    nm = calc_p_nomut_window(mu, min_gap)
                    nm_expect = np.ones((1, 2, 1))
                    self.assertTrue(np.array_equal(nm, nm_expect))
                    nc = calc_p_noclose_given_ends(mu, nm)
                    nc_expect = np.ones((1, 1, 1))
                    self.assertTrue(np.array_equal(nc, nc_expect))

    def _check_1_cluster(self, x, x_expect):
        self.assertEqual(x.shape, x_expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(x, x_expect),
                                             np.isnan(x_expect))))

    def test_length_2_min_gap_1(self):
        mu = np.array([[0.1, 0.2]]).reshape((2, 1))
        min_gap = 1
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1.],
            [np.nan, 0.9, 0.8],
        ]).reshape((2, 3, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98],
            [np.nan, 1.],
        ]).reshape((2, 2, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_3_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 1
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
        ]).reshape((2, 4, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.926],
            [1., 1., 0.94],
            [np.nan, 1., 1.],
        ]).reshape((3, 3, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_3_min_gap_2(self):
        mu = np.array([[0.1, 0.2, 0.3]]).reshape((3, 1))
        min_gap = 2
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7],
            [np.nan, np.nan, 0.72, 0.56],
        ]).reshape((3, 4, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.902],
            [1., 1., 0.94],
            [1., 1., 1.],
        ]).reshape((3, 3, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_length_4_min_gap_1(self):
        mu = np.array([[0.1, 0.2, 0.3, 0.4]]).reshape((4, 1))
        min_gap = 1
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
        ]).reshape((2, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
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
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
        ]).reshape((3, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
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
        nm = calc_p_nomut_window(mu, min_gap)
        nm_expect = np.array([
            [1., 1., 1., 1., 1.],
            [np.nan, 0.9, 0.8, 0.7, 0.6],
            [np.nan, np.nan, 0.72, 0.56, 0.42],
            [np.nan, np.nan, np.nan, 0.504, 0.336],
        ]).reshape((4, 5, 1))
        self._check_1_cluster(nm, nm_expect)
        nc = calc_p_noclose_given_ends(mu, nm)
        nc_expect = np.array([
            [1., 0.98, 0.902, 0.7428],
            [1., 1., 0.94, 0.788],
            [1., 1., 1., 0.88],
            [1., 1., 1., 1.],
        ]).reshape((4, 4, 1))
        self._check_1_cluster(nc, nc_expect)

    def test_clusters(self):
        rng = np.random.default_rng(seed=0)
        for ncls in range(4):
            for npos in range(5):
                for min_gap in range(npos):
                    with self.subTest(npos=npos, ncls=ncls, min_gap=min_gap):
                        mu = rng.random((npos, ncls))
                        nm = calc_p_nomut_window(mu, min_gap)
                        self.assertEqual(nm.shape, (min_gap + 1, npos + 1, ncls))
                        nc = calc_p_noclose_given_ends(mu, nm)
                        self.assertEqual(nc.shape, (npos, npos, ncls))
                        for k in range(ncls):
                            mu_k = mu[:, k].reshape((npos, 1))
                            nm_k = calc_p_nomut_window(mu_k, min_gap)
                            self.assertEqual(nm_k.shape,
                                             (min_gap + 1, npos + 1, 1))
                            self.assertTrue(np.allclose(np.triu(nm[:, :, k]),
                                                        np.triu(nm_k[:, :, 0])))
                            nc_k = calc_p_noclose_given_ends(mu_k, nm_k)
                            self.assertEqual(nc_k.shape,
                                             (npos, npos, 1))
                            self.assertTrue(np.allclose(np.triu(nc[:, :, k]),
                                                        np.triu(nc_k[:, :, 0])))


class TestCalcPNoCloseGivenEndsAuto(ut.TestCase):

    def test_1_dim(self):
        rng = np.random.default_rng(seed=0)
        max_n = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for gap in range(max_g + 1):
                mu = rng.random(n_pos)
                p_noclose_given_ends = calc_p_noclose_given_ends_auto(mu,
                                                                      gap)
                self.assertIsInstance(p_noclose_given_ends, np.ndarray)
                self.assertEqual(p_noclose_given_ends.shape, (n_pos, n_pos))

    def test_2_dim(self):
        rng = np.random.default_rng(seed=0)
        max_n = 5
        max_c = 5
        max_g = 4
        for n_pos in range(max_n + 1):
            for n_clust in range(max_c + 1):
                for gap in range(max_g + 1):
                    mu = rng.random((n_pos, n_clust))
                    p_noclose_given_ends = calc_p_noclose_given_ends_auto(mu,
                                                                          gap)
                    self.assertIsInstance(p_noclose_given_ends, np.ndarray)
                    self.assertEqual(p_noclose_given_ends.shape,
                                     (n_pos, n_pos, n_clust))

    def test_invalid_dim(self):
        rng = np.random.default_rng(seed=0)
        for n_dim in range(5):
            if n_dim == 1 or n_dim == 2:
                # Skip the dimensions that are valid.
                continue
            err_msg = ("p_mut_given_span must have 2 dimensions, "
                       f"but got {n_dim}")
            for size in range(5):
                dims = (size,) * n_dim
                for gap in range(5):
                    self.assertRaisesRegex(ValueError,
                                           err_msg,
                                           calc_p_noclose_given_ends,
                                           rng.random(dims),
                                           gap)


class TestCalcPMutGivenSpanDropped(ut.TestCase):

    def test_simulated(self):
        from scipy.stats import binom
        confidence = 0.999
        # Simulate reads.
        n_pos = 6
        n_reads = 1_000_000
        p_mut_max = 0.5
        p_mut, p_ends, _ = simulate_params(n_pos, 1, p_mut_max)
        muts, info, end5s, end3s = simulate_reads(n_reads, p_mut[:, 0], p_ends)
        # Test each minimum gap between mutations (min_gap).
        for min_gap in [0, 3, n_pos]:
            with self.subTest(min_gap=min_gap):
                # Determine which reads have no two mutations too close.
                has_no_close = label_no_close_muts(muts, min_gap)
                # Calculate the theoretical probability of no mutations
                # in each window.
                p_nomut_window_theory = calc_p_nomut_window(p_mut, min_gap)
                # Calculate the theoretical probability of no mutations
                # being too close for each pair of end coordinates.
                p_noclose_given_ends_theory = calc_p_noclose_given_ends(
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
                p_ends_given_clust_noclose_theory = calc_p_ends_given_clust_noclose(
                    p_ends,
                    p_noclose_given_ends_theory
                ).reshape((n_pos, n_pos))
                # Compare to the simulated proportion of reads aligning
                # to each pair of coordinates.
                p_noclose_given_clust = calc_p_noclose_given_clust(
                    p_ends, p_noclose_given_ends_theory
                )
                n_expect = np.round(p_ends_given_clust_noclose_theory
                                    * (p_noclose_given_clust * n_reads))
                inter_lo, inter_up = binom.interval(confidence,
                                                    n_expect,
                                                    p_ends_given_clust_noclose_theory)
                inter_lo = _triu_div(inter_lo, n_expect)
                inter_up = _triu_div(inter_up, n_expect)
                end5s_noclose = end5s[has_no_close]
                end3s_noclose = end3s[has_no_close]
                for end5 in range(n_pos):
                    for end3 in range(end5, n_pos):
                        p_ends_given_clust_noclose_simulated = np.mean(np.logical_and(
                            end5s_noclose == end5,
                            end3s_noclose == end3
                        ))
                        self.assertGreaterEqual(p_ends_given_clust_noclose_simulated,
                                                inter_lo[end5, end3])
                        self.assertLessEqual(p_ends_given_clust_noclose_simulated,
                                             inter_up[end5, end3])
                # Calculate the theoretical probability of a mutation at
                # each position given no mutations are too close.
                p_mut_given_span_dropped_theory = calc_p_mut_given_span_dropped(
                    p_mut,
                    p_ends,
                    p_noclose_given_ends_theory,
                    p_nomut_window_theory
                ).reshape(n_pos)
                # Compare the simulated probabilities of mutations given
                # no mutations are too close.
                info_sum = np.sum(info[has_no_close], axis=0)
                p_mut_given_span_dropped_simulated = (
                    np.sum(muts[has_no_close], axis=0) / info_sum
                )
                inter_lo, inter_up = binom.interval(confidence,
                                                    info_sum,
                                                    p_mut_given_span_dropped_theory)
                self.assertTrue(np.all(p_mut_given_span_dropped_simulated
                                       >= inter_lo / info_sum))
                self.assertTrue(np.all(p_mut_given_span_dropped_simulated
                                       <= inter_up / info_sum))

    def test_clusters(self):
        n_pos = 16
        max_clusters = 3
        max_gap = 3
        # Test each number of clusters.
        for n_cls in range(1, max_clusters + 1):
            p_mut, p_ends, _ = simulate_params(n_pos, n_cls)
            # Test each minimum gap between mutations.
            for min_gap in range(max_gap + 1):
                with self.subTest(n_cls=n_cls, min_gap=min_gap):
                    # Compute the mutation rates with no two mutations
                    # too close for all clusters simultaneously.
                    p_nomut_window = calc_p_nomut_window(p_mut, min_gap)
                    p_noclose_given_ends = calc_p_noclose_given_ends(
                        p_mut, p_nomut_window
                    )
                    p_mut_given_span_noclose = calc_p_mut_given_span_dropped(
                        p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                    )
                    # Compare to computing each cluster individually.
                    for k in range(n_cls):
                        self.assertTrue(np.allclose(
                            p_mut_given_span_noclose[:, k],
                            calc_p_mut_given_span_dropped(
                                p_mut[:, [k]],
                                p_ends,
                                p_noclose_given_ends[:, :, [k]],
                                p_nomut_window[:, :, [k]],
                            ).reshape(n_pos)
                        ))


class TestCalcPMutGivenSpanMerged(ut.TestCase):

    @staticmethod
    def _simulate_after_merge(p_before: np.ndarray,
                              p_ends: np.ndarray,
                              min_gap: int,
                              n_reads: int):
        # p_before shape: (n_pos, n_cls); simulate each column independently.
        _, n_cls = p_before.shape
        # Simulate mutations and merge them.
        result = np.empty_like(p_before)
        for k in range(n_cls):
            muts, info, _, _ = simulate_reads(n_reads, 
                                              p_before[:, k],
                                              p_ends)
            muts_merged = merge_mutations_right_to_left(muts, min_gap)
            result[:, k] = (np.count_nonzero(muts_merged, axis=0)
                            / np.count_nonzero(info, axis=0))
        return result
    
    def assert_not_statistically_significant(self,
                                             p_expect: np.ndarray,
                                             p_result: np.ndarray,
                                             coverage: np.ndarray,
                                             z_max: float = 6.0):
        self.assertEqual(p_result.shape, p_expect.shape)
        variance = p_expect * (1.0 - p_expect)
        se = np.sqrt(np.maximum(variance, 1.e-12) / coverage)
        z = np.abs(p_result - p_expect) / se
        self.assertTrue(
            np.all(z <= z_max),
            msg=(f"Max z-score {np.max(z):.3f} exceeded {z_max}"
                 f"; max abs diff = {np.max(np.abs(p_result - p_expect)):.3e}")
        )

    def test_matches_simulated(self):
        cases = [
            (1, 1, 0, 201),
            (5, 1, 1, 202),
            (7, 2, 2, 203),
            (8, 3, 4, 204),
        ]
        n_reads = 100_000
        for n_pos, n_cls, min_gap, seed in cases:
            with self.subTest(n_pos=n_pos,
                              n_cls=n_cls,
                              min_gap=min_gap,
                              seed=seed):
                case_rng = np.random.default_rng(seed)
                p_mut_given_span = case_rng.random((n_pos, n_cls))
                p_ends = np.triu(case_rng.random((n_pos, n_pos)))
                p_ends /= p_ends.sum()
                p_expect = _calc_p_mut_given_span_merged(p_mut_given_span,
                                                         p_ends,
                                                         min_gap)
                p_result = self._simulate_after_merge(
                    p_mut_given_span,
                    p_ends,
                    min_gap,
                    n_reads=n_reads
                )
                coverage = n_reads * np.array([p_ends[:j + 1, j:].sum()
                                               for j in range(n_pos)])
                self.assert_not_statistically_significant(p_expect,
                                                          p_result,
                                                          coverage[:, np.newaxis])


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
        self.assertTrue(triu_allclose(result, expect))

    def test_slice_5x5(self):
        p_ends = np.arange(25.).reshape((5, 5))
        p_ends_cumsum = _triu_cumsum(p_ends)
        result = _slice_p_ends(p_ends, p_ends_cumsum, 1, 4)
        expect = np.array([[7., 9., 24.],
                           [0., 12., 27.],
                           [0., 0., 37.]])
        self.assertTrue(triu_allclose(result, expect))


class TestFindSplitPositions(ut.TestCase):

    def test_0(self):
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(np.array([[]]), min_gap, 0.),
                np.array([], dtype=int)
            ))

    def test_thresh0(self):
        rng = np.random.default_rng(seed=0)
        p_mut = 1. - rng.random((10, 2))
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(p_mut, min_gap, 0.),
                np.array([], dtype=int)
            ))

    def test_thresh1(self):
        rng = np.random.default_rng(seed=0)
        p_mut = 1. - rng.random((10, 2))
        for min_gap in range(4):
            self.assertTrue(np.array_equal(
                _find_split_positions(p_mut, min_gap, 1.),
                np.array([], dtype=int)
            ))

    def test_gap0(self):
        rng = np.random.default_rng(seed=0)
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

        rng = np.random.default_rng(seed=0)
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

    REL_TOL = 0.01

    @staticmethod
    def random_params(npos: int, ncls: int, rng: np.random.Generator):
        p_mut = rng.beta(0.5, 9.5, (npos, ncls))
        p_ends = 1. - rng.random((npos, npos, ncls))
        p_clust = 1. - rng.random(ncls)
        return p_mut, p_ends, p_clust

    def test_quick_unbias(self):
        rng = np.random.default_rng(seed=0)
        n_pos = 256
        for t in [0., 0.001]:
            for n_cls in [1, 2]:
                for min_gap in [1, 3]:
                    for f_below in [0.1, 0.5, 0.9]:
                        p_mut, p_ends, p_clust = self.random_params(n_pos, n_cls, rng)
                        # Set some mutation rates to the threshold.
                        n_below = round(n_pos * f_below)
                        t_pos = rng.choice(n_pos, n_below, replace=False)
                        equals_t = np.zeros(n_pos, dtype=bool)
                        equals_t[t_pos] = True
                        p_mut[equals_t] = t
                        for mut_collisions in ["drop", "merge"]:
                            # Compare quick vs. exact unbias.
                            (p_mut_quick,
                            p_ends_quick,
                            p_clust_quick) = calc_params(
                                p_mut, p_ends, p_clust, min_gap, mut_collisions,
                                quick_unbias=True,
                                quick_unbias_thresh=t
                            )
                            (p_mut_exact,
                            p_ends_exact,
                            p_clust_exact) = calc_params(
                                p_mut, p_ends, p_clust, min_gap, mut_collisions,
                                quick_unbias=False,
                                quick_unbias_thresh=t
                            )
                            with self.subTest(
                                threshold=t,
                                n_cls=n_cls,
                                min_gap=min_gap,
                                f_below=f_below,
                                mut_collisions=mut_collisions
                            ):
                                self.assertTrue(
                                    np.allclose(p_mut_quick,
                                                p_mut_exact,
                                                atol=(2. * t),
                                                rtol=self.REL_TOL)
                                )
                                self.assertTrue(
                                    np.allclose(np.triu(p_ends_quick),
                                                np.triu(p_ends_exact),
                                                atol=(2. * t),
                                                rtol=self.REL_TOL)
                                )
                                self.assertTrue(
                                    np.allclose(p_clust_quick,
                                                p_clust_exact,
                                                atol=(2. * t),
                                                rtol=self.REL_TOL)
                                )


class TestCalcPMutGivenSpan(ut.TestCase):

    def test_invert(self):
        """ Test the inverse of `calc_p_mut_given_span_noclose`. """
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, _ = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the mutation rates without mutations too close.
                p_nomut_window = calc_p_nomut_window(p_mut, min_gap)
                p_noclose_given_ends = calc_p_noclose_given_ends(
                    p_mut, p_nomut_window
                )
                p_mut_given_span_dropped = calc_p_mut_given_span_dropped(
                    p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                )
                p_mut_given_span = calc_p_mut_given_span(
                    p_mut_given_span_dropped,
                    p_ends,
                    min_gap,
                    "drop",
                    p_mut_given_span_dropped,
                )
                self.assertEqual(p_mut_given_span.shape, p_mut.shape)
                self.assertTrue(np.allclose(p_mut_given_span,
                                            p_mut,
                                            atol=1.e-4,
                                            rtol=1.e-3))


class TestCalcPEnds(ut.TestCase):

    def test_invert(self):
        """ Test the inverse of `calc_p_ends_given_clust_noclose`. """
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, p_clust = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the coordinate distributions without mutations
                # too close.
                p_noclose_given_ends = calc_p_noclose_given_ends(
                    p_mut, calc_p_nomut_window(p_mut, min_gap)
                )
                p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
                    p_ends, p_noclose_given_ends
                )
                # Infer the original distribution of all reads.
                p_ends_inferred = calc_p_ends(p_ends_given_clust_noclose,
                                              p_noclose_given_ends,
                                              p_mut,
                                              p_clust)
                self.assertEqual(p_ends_inferred.shape, p_ends.shape)
                self.assertTrue(triu_allclose(p_ends_inferred, p_ends))


class TestCalcPNoCloseGivenClust(ut.TestCase):

    def test_p_noclose_given_clust(self):
        p_ends = np.array([[0.2, 0.5],
                           [0.4, 0.3]])
        p_noclose_given_ends = np.array([[[0.2, 0.4], [0.1, 0.2]],
                                         [[0.3, 0.6], [0.4, 0.8]]])
        expect = np.array([0.21, 0.42])
        result = calc_p_noclose_given_clust(p_ends, p_noclose_given_ends)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcPNoClose(ut.TestCase):

    def test_p_noclose(self):
        p_clust = np.array([0.6, 0.3, 0.1])
        p_noclose_given_clust = np.array([0.9, 0.8, 0.7])
        p_noclose = calc_p_noclose(p_clust, p_noclose_given_clust)
        self.assertIsInstance(p_noclose, float)
        self.assertTrue(np.isclose(p_noclose, 0.85))


class TestCalcPClust(ut.TestCase):

    def test_p_clust(self):
        p_clust_given_observed = np.array([0.2, 0.5, 0.3])
        p_noclose_given_clust = np.array([0.8, 0.4, 0.5])
        expect = np.array([0.25, 1.25, 0.60]) / 2.1
        result = calc_p_clust(p_clust_given_observed, p_noclose_given_clust)
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
        numer = np.array([[[0.7, 0.8], [1.0, 1.2], [0.3, 0.4]],
                          [[np.nan, np.nan], [0.5, 0.6], [0.4, 0.5]],
                          [[np.nan, np.nan], [np.nan, np.nan], [0.7, 0.9]]])
        denom = np.array([3.6, 4.4])
        expect = numer / denom
        result = calc_p_ends_observed(npos, end5s, end3s, weights)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.all(np.logical_or(np.isclose(result, expect),
                                             np.isnan(expect))))


class TestCalcParams(ut.TestCase):

    def test_infer(self):
        n_pos = 16
        max_p_mut = 0.5
        max_clusters = 3
        max_gap = 3
        for n_cls in range(1, max_clusters):
            p_mut, p_ends, p_clust = simulate_params(n_pos, n_cls, max_p_mut)
            for min_gap in range(max_gap + 1):
                # Compute the mutation rates and distributions of end
                # coordinates without mutations too close.
                p_nomut_window = calc_p_nomut_window(p_mut, min_gap)
                p_noclose_given_ends = calc_p_noclose_given_ends(
                    p_mut, p_nomut_window
                )
                p_noclose_given_clust = calc_p_noclose_given_clust(
                    p_ends, p_noclose_given_ends
                )
                p_mut_given_span_noclose = calc_p_mut_given_span_dropped(
                    p_mut, p_ends, p_noclose_given_ends, p_nomut_window
                )
                p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
                    p_ends, p_noclose_given_ends
                )
                p_clust_given_noclose = calc_p_clust_given_noclose(
                    p_clust, p_noclose_given_clust
                )
                # Infer the original parameters using those of the reads
                # without mutations too close.
                (p_mut_inferred,
                 p_ends_inferred,
                 p_clust_inferred) = calc_params(
                    p_mut_given_span_noclose,
                    p_ends_given_clust_noclose,
                    p_clust_given_noclose,
                    min_gap,
                    "drop"
                )
                self.assertEqual(p_mut_inferred.shape, p_mut.shape)
                self.assertTrue(np.allclose(p_mut_inferred,
                                            p_mut,
                                            atol=1.e-4,
                                            rtol=1.e-2))
                self.assertEqual(p_ends_inferred.shape, p_ends.shape)
                self.assertTrue(triu_allclose(p_ends_inferred, p_ends))
                self.assertEqual(p_clust_inferred.shape, p_clust.shape)
                self.assertTrue(np.allclose(p_clust_inferred,
                                            p_clust,
                                            atol=1.e-3,
                                            rtol=1.e-2))


class TestCalcRectangularSum(ut.TestCase):

    @staticmethod
    def calc_spanning_sum_slow(array: np.ndarray):
        spanning_sum = np.empty(array.shape[1:])
        for j in range(array.shape[0]):
            spanning_sum[j] = np.sum(array[:(j + 1), j:], axis=(0, 1))
        return spanning_sum

    def test_2d(self):
        rng = np.random.default_rng(seed=0)
        for n in range(5):
            with self.subTest(n=n):
                array = rng.random((n, n))
                fast_sum = calc_rectangular_sum(array)
                slow_sum = self.calc_spanning_sum_slow(array)
                self.assertEqual(fast_sum.shape, (n,))
                self.assertEqual(slow_sum.shape, (n,))
                self.assertTrue(np.allclose(fast_sum, slow_sum))

    def test_3d(self):
        rng = np.random.default_rng(seed=0)
        for n in range(5):
            for k in range(3):
                with self.subTest(n=n, k=k):
                    array = rng.random((n, n, k))
                    fast_sum = calc_rectangular_sum(array)
                    slow_sum = self.calc_spanning_sum_slow(array)
                    self.assertEqual(fast_sum.shape, (n, k))
                    self.assertEqual(slow_sum.shape, (n, k))
                    self.assertTrue(np.allclose(fast_sum, slow_sum))


class TestCalcRectangularSumWeighted(ut.TestCase):

    @staticmethod
    def calc_weighted_slow(array: np.ndarray, weights: np.ndarray):
        return calc_rectangular_sum(weights[:, :, np.newaxis] * array)

    def test_3d(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(5):
            for ncls in range(3):
                with self.subTest(npos=npos, ncls=ncls):
                    array = rng.random((npos, npos, ncls))
                    weights = rng.random((npos, npos))
                    fast = _calc_rectangular_sum_weighted(array, weights)
                    slow = self.calc_weighted_slow(array, weights)
                    self.assertEqual(fast.shape, (npos, ncls))
                    self.assertTrue(np.allclose(fast, slow))

    def test_unit_weights(self):
        rng = np.random.default_rng(seed=0)
        for npos in range(5):
            for ncls in range(3):
                with self.subTest(npos=npos, ncls=ncls):
                    array = rng.random((npos, npos, ncls))
                    weights = np.ones((npos, npos))
                    fast = _calc_rectangular_sum_weighted(array, weights)
                    unweighted = calc_rectangular_sum(array)
                    self.assertEqual(fast.shape, (npos, ncls))
                    self.assertTrue(np.allclose(fast, unweighted))


if __name__ == "__main__":
    ut.main(verbosity=2)
