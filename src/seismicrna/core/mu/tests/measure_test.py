import unittest as ut

import numpy as np

from seismicrna.core.mu.measure import calc_gini, calc_signal_noise

rng = np.random.default_rng()


class TestCalcGini(ut.TestCase):

    def test_all_zero_1d(self):
        for npos in range(1, 5):
            self.assertTrue(np.isnan(calc_gini(np.zeros(npos))))

    def test_all_equal_1d(self):
        for npos in range(1, 5):
            for x in np.linspace(0.001, 1., 5):
                self.assertEqual(float(calc_gini(np.full(npos, x))), 0.)

    def test_all_zero_2d(self):
        for npos in range(1, 5):
            for ncls in range(1, 5):
                self.assertTrue(np.all(np.isnan(
                    calc_gini(np.zeros((npos, ncls))))
                ))

    def test_all_equal_2d(self):
        for npos in range(1, 5):
            for ncls in range(1, 5):
                self.assertTrue(np.array_equal(
                    calc_gini(np.broadcast_to(
                        np.linspace(0.001, 1., ncls)[np.newaxis, :],
                        (npos, ncls)
                    )),
                    np.zeros(ncls)
                ))

    def test_monopoly_1d(self):
        for npos in range(1, 11):
            for pos in range(npos):
                values = np.zeros(npos)
                for x in np.linspace(0.001, 1., 5):
                    values[pos] = x
                    self.assertTrue(np.isclose(float(calc_gini(values)),
                                               1. - 1. / npos))

    def test_distributed_1d(self):
        self.assertTrue(np.isclose(float(calc_gini(np.linspace(0., 1., 11))),
                                   4. / 11.))

    def test_empty_1d(self):
        self.assertTrue(np.isnan(calc_gini(np.empty(0))))


class TestCalcSignalNoise(ut.TestCase):

    def test_1d(self):
        for ns in range(1, 5):
            signal = 1. - rng.random(ns)
            smean = signal.mean()
            for nn in range(1, 5):
                npos = ns + nn
                noise = 1. - rng.random(nn)
                nmean = noise.mean()
                is_signal = rng.permutation(npos) < ns
                values = np.empty(npos)
                values[is_signal] = signal
                values[~is_signal] = noise
                self.assertTrue(np.isclose(calc_signal_noise(values,
                                                             is_signal),
                                           smean / nmean))

    def test_2d(self):
        for ncls in range(1, 5):
            for ns in range(1, 5):
                signal = 1. - rng.random((ns, ncls))
                smean = signal.mean(axis=0)
                for nn in range(1, 5):
                    npos = ns + nn
                    noise = 1. - rng.random((nn, ncls))
                    nmean = noise.mean(axis=0)
                    is_signal = rng.permutation(npos) < ns
                    values = np.empty((npos, ncls))
                    values[is_signal] = signal
                    values[~is_signal] = noise
                    self.assertTrue(np.allclose(calc_signal_noise(values,
                                                                  is_signal),
                                                smean / nmean))

    def test_empty_1d(self):
        self.assertTrue(np.isnan(calc_signal_noise(np.empty(0),
                                                   np.empty(0, dtype=bool))))


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
