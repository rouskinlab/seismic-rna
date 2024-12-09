import unittest as ut
from itertools import product

import numpy as np
from scipy.stats import binom

from seismicrna.core.random import stochastic_round

rng = np.random.default_rng()


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
        n_trials = 10000
        floor = 8
        for dims in [(), (0,), (1,), (2,), (5,), (5, 5)]:
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
