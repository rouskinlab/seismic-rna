import unittest as ut
import warnings
from itertools import permutations

import numpy as np

from seismicrna.cluster.emk import assign_clusterings
from seismicrna.core.array import calc_inverse

rng = np.random.default_rng()


class TestAssignClusterings(ut.TestCase):

    def test_0_clusters(self):
        for n in range(5):
            x = np.empty((n, 0))
            y = np.empty((n, 0))
            self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                           np.empty(0, dtype=int)))

    def test_0_positions(self):
        for n in range(5):
            x = np.empty((0, n))
            y = np.empty((0, n))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                with np.errstate(invalid="ignore"):
                    self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                                   np.arange(n)))

    def test_1_cluster(self):
        for n in range(1, 10):
            x = rng.random((n, 1))
            y = rng.random((n, 1))
            self.assertTrue(np.array_equal(assign_clusterings(x, y),
                                           np.array([0])))

    def test_more_clusters(self):
        for ncls in range(2, 5):
            for npos in range(2, 8):
                with self.subTest(ncls=ncls, npos=npos):
                    x = rng.random((npos, ncls))
                    for permutation in permutations(range(ncls)):
                        cols = np.array(permutation)
                        y = x[:, cols]
                        self.assertTrue(np.array_equal(
                            assign_clusterings(x, y),
                            calc_inverse(cols)
                        ))
                        self.assertTrue(np.array_equal(
                            assign_clusterings(y, x),
                            cols
                        ))


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
