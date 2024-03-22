import unittest as ut

import numpy as np

from seismicrna.cluster.em import _expectation


class TestExpectation(ut.TestCase):

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
                expect_memberships: np.ndarray):
        result_log_marginals, result_memberships = _expectation(p_mut,
                                                                p_ends,
                                                                p_clust,
                                                                end5s,
                                                                end3s,
                                                                unmasked,
                                                                muts_per_pos,
                                                                min_mut_gap)
        # print("EXPECT MARGINALS")
        # print(np.exp(expect_log_marginals))
        # print("RESULT MARGINALS")
        # print(np.exp(result_log_marginals))
        self.assertEqual(result_log_marginals.shape, expect_log_marginals.shape)
        self.assertTrue(np.allclose(result_log_marginals, expect_log_marginals))
        self.assertEqual(expect_memberships.shape, result_memberships.shape)
        self.assertTrue(np.allclose(expect_memberships, result_memberships))

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
        expect_memberships = np.ones((2, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

    def test_2pos(self):
        p_mut = np.array([[0.1],
                          [0.2]])
        p_ends = np.array([[0.2, 0.5],
                           [0.0, 0.3]])
        p_clust = np.array([1.])
        end5s = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        end3s = np.array([0, 0, 1, 1, 1, 1, 1, 1])
        unmasked = np.array([0, 1])
        muts_per_pos = [np.array([1, 3, 5]), np.array([4, 5, 7])]
        min_mut_gap = 0
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.36, 0.04, 0.09, 0.01, 0.24, 0.06]
        )
        expect_memberships = np.ones((8, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

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
        min_mut_gap = 0
        expect_log_marginals = np.log(
            [0.20, 0.20, 0.40, 0.40, 0.10, 0.10, 0.24, 0.06]
        )
        expect_memberships = np.ones((8, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)

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
        min_mut_gap = 0
        expect_log_marginals = np.log(
            [0.18, 0.02, 0.45, 0.05, 0.45, 0.05, 0.30, 0.30]
        )
        expect_memberships = np.ones((8, 1))
        self.compare(p_mut,
                     p_ends,
                     p_clust,
                     end5s,
                     end3s,
                     unmasked,
                     muts_per_pos,
                     min_mut_gap,
                     expect_log_marginals,
                     expect_memberships)


if __name__ == "__main__":
    ut.main()
