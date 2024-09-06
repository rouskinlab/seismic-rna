import unittest as ut
from itertools import product

import numpy as np
from scipy.stats import binom

from seismicrna.core.random import (combinations_array,
                                    combinations_bools,
                                    stochastic_round)

rng = np.random.default_rng()


class TestCombinationsArray(ut.TestCase):

    def test_neg_n(self):
        for n in [-2, -1]:
            for m in [-2, -1, 1, 2]:
                self.assertRaisesRegex(ValueError,
                                       f"n must be ≥ 0, but got {n}",
                                       combinations_array,
                                       n, m)

    def test_0_m(self):
        for m in [-2, -1, 1, 2]:
            self.assertRaisesRegex(ValueError,
                                   f"m must be ≥ 0 and ≤ 0, but got {m}",
                                   combinations_array,
                                   0, m)

    def test_n_0(self):
        for n in range(5):
            result = combinations_array(n, 0)
            expect = np.zeros((1, 0), dtype=bool)
            self.assertTrue(np.array_equal(result, expect))

    def test_n_n(self):
        for n in range(5):
            result = combinations_array(n, n)
            expect = np.arange(n).reshape((1, n))
            self.assertTrue(np.array_equal(result, expect))

    def test_n_1(self):
        for n in range(1, 5):
            result = combinations_array(n, 1)
            expect = np.arange(n).reshape((n, 1))
            self.assertTrue(np.array_equal(result, expect))

    def test_3_2(self):
        result = combinations_array(3, 2)
        expect = np.array([[0, 1],
                           [0, 2],
                           [1, 2]])
        self.assertTrue(np.array_equal(result, expect))

    def test_4_2(self):
        result = combinations_array(4, 2)
        expect = np.array([[0, 1],
                           [0, 2],
                           [0, 3],
                           [1, 2],
                           [1, 3],
                           [2, 3]])
        self.assertTrue(np.array_equal(result, expect))

    def test_4_3(self):
        result = combinations_array(4, 3)
        expect = np.array([[0, 1, 2],
                           [0, 1, 3],
                           [0, 2, 3],
                           [1, 2, 3]])
        self.assertTrue(np.array_equal(result, expect))


class TestCombinationsBools(ut.TestCase):

    def test_1_0(self):
        result = combinations_bools(1, 0)
        expect = np.array([[False]])
        self.assertTrue(np.array_equal(result, expect))

    def test_1_1(self):
        result = combinations_bools(1, 1)
        expect = np.array([[True]])
        self.assertTrue(np.array_equal(result, expect))

    def test_2_0(self):
        result = combinations_bools(2, 0)
        expect = np.array([[False, False]])
        self.assertTrue(np.array_equal(result, expect))

    def test_2_1(self):
        result = combinations_bools(2, 1)
        expect = np.array([[True, False],
                           [False, True]])
        self.assertTrue(np.array_equal(result, expect))

    def test_2_2(self):
        result = combinations_bools(2, 2)
        expect = np.array([[True, True]])
        self.assertTrue(np.array_equal(result, expect))

    def test_3_0(self):
        result = combinations_bools(3, 0)
        expect = np.array([[False, False, False]])
        self.assertTrue(np.array_equal(result, expect))

    def test_3_1(self):
        result = combinations_bools(3, 1)
        expect = np.array([[True, False, False],
                           [False, True, False],
                           [False, False, True]])
        self.assertTrue(np.array_equal(result, expect))

    def test_3_2(self):
        result = combinations_bools(3, 2)
        expect = np.array([[True, True, False],
                           [True, False, True],
                           [False, True, True]])
        self.assertTrue(np.array_equal(result, expect))

    def test_3_3(self):
        result = combinations_bools(3, 3)
        expect = np.array([[True, True, True]])
        self.assertTrue(np.array_equal(result, expect))


class TestStochasticRound(ut.TestCase):

    def test_int_arrays(self):
        n_trials = 10
        max_ndim = 3
        max_dim = 3
        low = 0
        high = 10
        for preserve_sum in [False, True]:
            for ndim in range(max_ndim):
                for dims in product(range(max_dim), repeat=ndim):
                    for _ in range(n_trials):
                        values = rng.integers(low, high, dims)
                        self.assertTrue(np.array_equal(
                            stochastic_round(values, preserve_sum),
                            values
                        ))

    def test_float_arrays(self):
        confidence = 0.9995
        n_trials = 10000
        floor = 8
        for dims in [(), (0,), (1,), (2,), (100,), (10, 10)]:
            mantissas = rng.random(dims)
            values = floor + mantissas
            trials = [stochastic_round(values)
                      for _ in range(n_trials)]
            # Check that every value was rounded to one of the nearest
            # integers.
            for trial in trials:
                self.assertTrue(np.all(np.logical_or(trial == floor,
                                                     trial == floor + 1)))
            # Check the fraction of trials that rounded each value down.
            floor_ci_lo, floor_ci_up = binom.interval(
                confidence,
                n=n_trials,
                p=(1. - mantissas)
            )
            n_floor = np.sum(np.array(trials) == floor, axis=0)
            self.assertTrue(np.all(np.greater_equal(n_floor, floor_ci_lo)))
            self.assertTrue(np.all(np.less_equal(n_floor, floor_ci_up)))

    def test_float_arrays_preserve_sum(self):
        confidence = 0.9995
        n_trials = 5000
        floor = 8
        for dims in [(), (0,), (1,), (2,), (4,), (2, 2)]:
            mantissas = rng.random(dims)
            values = floor + mantissas
            trials = [stochastic_round(values, preserve_sum=True)
                      for _ in range(n_trials)]
            # Check that every value was rounded to one of the nearest
            # integers.
            for trial in trials:
                self.assertTrue(np.all(np.logical_or(trial == floor,
                                                     trial == floor + 1)))
            # Check the fraction of trials that rounded each value down.
            floor_ci_lo, floor_ci_up = binom.interval(
                confidence,
                n=n_trials,
                p=(1. - mantissas)
            )
            n_floor = np.sum(np.array(trials) == floor, axis=0)
            self.assertTrue(np.all(np.greater_equal(n_floor, floor_ci_lo)))
            self.assertTrue(np.all(np.less_equal(n_floor, floor_ci_up)))
            # Check that every trial sums to either values_sum_floor or
            # (values_sum_floor + 1), due to preserve_sum=True.
            values_sum = values.sum()
            values_sum_floor = int(values_sum)
            sums = np.array([trial.sum() for trial in trials])
            self.assertTrue(np.all(np.logical_or(sums == values_sum_floor,
                                                 sums == values_sum_floor + 1)))
            # Check the fraction of trials that sum to values_sum_floor
            # is within the expected range, due to preserve_sum=True.
            floor_ci_lo, floor_ci_up = binom.interval(
                confidence,
                n=n_trials,
                p=(1. - (values_sum - values_sum_floor))
            )
            n_floor = np.count_nonzero(sums == values_sum_floor)
            self.assertGreaterEqual(n_floor, floor_ci_lo)
            self.assertLessEqual(n_floor, floor_ci_up)
