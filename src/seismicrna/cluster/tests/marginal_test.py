import unittest as ut

import numpy as np

from seismicrna.cluster.marginal import calc_marginal_resps


class TestMarginalResps(ut.TestCase):

    def compare(self,
                p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                end5s: np.ndarray,
                end3s: np.ndarray,
                unmasked: np.ndarray,
                muts_per_pos: list[np.ndarray],
                min_mut_gap: int,
                expect_log_marginals: np.ndarray,
                expect_resps: np.ndarray):
        result_log_marginals, result_resps = calc_marginal_resps(p_mut,
                                                                 p_ends,
                                                                 p_clust,
                                                                 end5s,
                                                                 end3s,
                                                                 unmasked,
                                                                 muts_per_pos,
                                                                 min_mut_gap)
        self.assertEqual(result_log_marginals.shape, expect_log_marginals.shape)
        self.assertTrue(np.allclose(result_log_marginals, expect_log_marginals))
        self.assertEqual(expect_resps.shape, result_resps.shape)
        self.assertTrue(np.allclose(expect_resps, result_resps))

    def test_1pos(self):
        p_mut = np.array([[0.1]])
        p_ends = np.array([[1.]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0])
        end3s = np.array([0, 0])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1])]
        min_mut_gap = 0
        expect_log_marginals = np.log([0.9, 0.1])
        expect_resps = np.ones((2, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_resps)

    def test_2pos_gap0(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3, 5]), np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.36, 0.04, 0.09, 0.01, 0.24, 0.06]
        )
        expect_resps = np.ones((8, 1))
        min_mut_gap = 0
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_resps)

    def test_2pos_gap1(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3]), np.array([4, 6])]
        expect_log_marginals = np.log([0.18 / 0.99,
                                       0.02 / 0.99,
                                       0.36 / 0.99,
                                       0.04 / 0.99,
                                       0.09 / 0.99,
                                       0.24 / 0.99,
                                       0.06 / 0.99])
        expect_resps = np.ones((7, 1))
        min_mut_gap = 1
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_resps)

    def test_2pos_masked0(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([1])
        muts_per_pos = [np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.20, 0.20, 0.40, 0.40, 0.10, 0.10, 0.24, 0.06]
        )
        expect_resps = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         expect_log_marginals,
                         expect_resps)

    def test_2pos_masked1(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1, 3, 5])]
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.45, 0.05, 0.45, 0.05, 0.30, 0.30]
        )
        expect_resps = np.ones((8, 1))
        for min_mut_gap in [0, 1]:
            self.compare(p_mut,
                         p_ends,
                         p_clust,
                         end5s,
                         end3s,
                         unmasked,
                         muts_per_pos,
                         min_mut_gap,
                         expect_log_marginals,
                         expect_resps)

    def test_1pos_2clusters(self):
        p_mut = np.array([[0.1, 0.2]])
        p_ends = np.array([[1.]])
        p_clust = np.array([0.4, 0.6])
        end5s = np.array([0, 0])
        end3s = np.array([0, 0])
        unmasked = np.array([0])
        muts_per_pos = [np.array([1])]
        min_mut_gap = 0
        expect_log_marginals = np.log([0.84, 0.16])
        expect_resps = np.array([[3. / 7., 4. / 7.],
                                 [1. / 4., 3. / 4.]])
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_resps)

    def test_2pos_gap0_2clusters(self):
        p_mut = np.array([[0.1, 0.4],
                          [0.2, 0.3]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([0.4, 0.6])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3, 5]), np.array([4, 5, 7])]
        expect_log_marginals = np.log(
            [0.144, 0.056, 0.270, 0.100, 0.090, 0.040, 0.222, 0.078]
        )
        expect_resps = np.array([[1. / 2., 1. / 2.],
                                 [1. / 7., 6. / 7.],
                                 [8. / 15., 7. / 15.],
                                 [4. / 25., 21. / 25.],
                                 [2. / 5., 3. / 5.],
                                 [1. / 10., 9. / 10.],
                                 [16. / 37., 21. / 37.],
                                 [4. / 13., 9. / 13.]])
        min_mut_gap = 0
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_resps)


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
