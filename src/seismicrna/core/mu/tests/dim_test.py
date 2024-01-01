import unittest as ut
from itertools import product

import numpy as np
import pandas as pd

from seismicrna.core.mu.dim import (count_pos,
                                    counts_pos,
                                    counts_pos_consensus)

rng = np.random.default_rng()


class TestCountPos(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "A 0-D array has no positional axis",
                               count_pos,
                               np.array(0.))

    def test_array1d(self):
        for n in range(5):
            self.assertEqual(count_pos(rng.random(n)), n)

    def test_array2d(self):
        for n in range(5):
            for m in range(3):
                self.assertEqual(count_pos(rng.random((n, m))), n)

    def test_series(self):
        for n in range(5):
            self.assertEqual(count_pos(pd.Series(rng.random(n))), n)

    def test_dataframe(self):
        for n in range(5):
            for m in range(3):
                self.assertEqual(count_pos(pd.DataFrame(rng.random((n, m)))), n)


class TestCountsPos(ut.TestCase):

    def test_array1d(self):
        for a in range(3):
            for ns in product(range(5), repeat=a):
                arrays = (rng.random(n) for n in ns)
                self.assertEqual(counts_pos(*arrays), ns)

    def test_series(self):
        for a in range(3):
            for ns in product(range(5), repeat=a):
                arrays = (pd.Series(rng.random(n)) for n in ns)
                self.assertEqual(counts_pos(*arrays), ns)

    def test_dataframe(self):
        for a in range(3):
            for ns in product(range(5), repeat=a):
                for m in range(3):
                    arrays = (pd.DataFrame(rng.random((n, m))) for n in ns)
                    self.assertEqual(counts_pos(*arrays), ns)


class TestCountsPosConsensus(ut.TestCase):

    def test_empty(self):
        self.assertRaisesRegex(ValueError,
                               "requires at least 1 array, but got 0 arrays",
                               counts_pos_consensus)

    def test_equal(self):
        for a in range(1, 5):
            for n in range(5):
                arrays = (rng.random(n) for _ in range(a))
                self.assertEqual(counts_pos_consensus(*arrays), n)

    def test_unequal(self):
        for n in range(5):
            array1 = rng.random(n)
            array2 = rng.random(n + 1)
            self.assertRaisesRegex(ValueError,
                                   "requires all arrays to have the same",
                                   counts_pos_consensus,
                                   array1,
                                   array2)


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
