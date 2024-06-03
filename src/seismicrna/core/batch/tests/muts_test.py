import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.ends import mask_segment_ends
from seismicrna.core.batch.muts import calc_muts_matrix
from seismicrna.core.rel.code import (DELET,
                                      MATCH,
                                      NOCOV,
                                      SUB_A,
                                      SUB_C,
                                      SUB_G,
                                      SUB_T)
from seismicrna.core.seq import DNA, Section


class TestCalcMutsMatrix(ut.TestCase):

    def test_full_reads_no_muts(self):
        for length in range(10):
            section = Section("myref", DNA.random(length))
            muts = dict()
            for num_reads in range(10):
                read_nums = np.arange(num_reads)
                seg_end5s = np.full((num_reads, 1), section.end5)
                seg_end3s = np.full((num_reads, 1), section.end3)
                with self.subTest(length=length, num_reads=num_reads):
                    result = calc_muts_matrix(section,
                                              read_nums,
                                              seg_end5s,
                                              seg_end3s,
                                              muts)
                    expect = pd.DataFrame(MATCH, read_nums, section.unmasked)
                    self.assertTrue(expect.equals(result))

    def test_full_reads_no_muts_some_masked(self):
        section = Section("myref", DNA("GTACTCAG"))
        section.mask_gu()
        muts = dict()
        for num_reads in range(10):
            read_nums = np.arange(num_reads)
            seg_end5s = np.full((num_reads, 1), section.end5)
            seg_end3s = np.full((num_reads, 1), section.end3)
            with self.subTest(num_reads=num_reads):
                result = calc_muts_matrix(section,
                                          read_nums,
                                          seg_end5s,
                                          seg_end3s,
                                          muts)
                expect = pd.DataFrame(MATCH, read_nums, section.unmasked)
                self.assertTrue(expect.equals(result))

    def test_partial_reads_no_muts(self):
        section = Section("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, NOCOV]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))

    def test_partial_reads_no_muts_some_masked(self):
        section = Section("myref", DNA("TAGCT"))
        section.mask_gu()
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV],
                               [NOCOV, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, NOCOV],
                               [NOCOV, NOCOV]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))

    def test_paired_reads_no_muts(self):
        section = Section("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        seg_end3s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, MATCH, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, MATCH, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, MATCH, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, MATCH, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))

    def test_paired_reads_masked_segments(self):
        section = Section("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        seg_end3s = np.array([[1, 0, 1, 0, 1, 0, 3, 0, 5],
                              [0, 2, 0, 4, 0, 5, 0, 5, 0]]).T
        seg_end5s, seg_end3s = mask_segment_ends(seg_end5s, seg_end3s)
        self.assertTrue(np.ma.is_masked(seg_end5s))
        self.assertTrue(np.ma.is_masked(seg_end3s))
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, MATCH, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, MATCH, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))

    def test_partial_reads_muts(self):
        section = Section("myref", DNA.random(5))
        muts = {1: {DELET: np.array([3]),
                    SUB_A: np.array([5]),
                    SUB_C: np.array([7]),
                    SUB_G: np.array([8]),
                    SUB_T: np.array([9])},
                2: {DELET: np.array([5]),
                    SUB_A: np.array([7]),
                    SUB_C: np.array([8]),
                    SUB_G: np.array([9]),
                    SUB_T: np.array([12])},
                3: {DELET: np.array([7]),
                    SUB_A: np.array([8]),
                    SUB_C: np.array([9]),
                    SUB_G: np.array([12]),
                    SUB_T: np.array([13])},
                4: {DELET: np.array([8]),
                    SUB_A: np.array([9]),
                    SUB_C: np.array([12]),
                    SUB_G: np.array([13]),
                    SUB_T: np.array([16])},
                5: {DELET: np.array([9]),
                    SUB_A: np.array([12]),
                    SUB_C: np.array([13]),
                    SUB_G: np.array([16]),
                    SUB_T: np.array([19])}}
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV, NOCOV, NOCOV, NOCOV],
                               [DELET, NOCOV, NOCOV, NOCOV, NOCOV],
                               [SUB_A, DELET, NOCOV, NOCOV, NOCOV],
                               [SUB_C, SUB_A, DELET, NOCOV, NOCOV],
                               [SUB_G, SUB_C, SUB_A, DELET, NOCOV],
                               [SUB_T, SUB_G, SUB_C, SUB_A, DELET],
                               [NOCOV, SUB_T, SUB_G, SUB_C, SUB_A],
                               [NOCOV, NOCOV, SUB_T, SUB_G, SUB_C],
                               [NOCOV, NOCOV, NOCOV, SUB_T, SUB_G],
                               [NOCOV, NOCOV, NOCOV, NOCOV, SUB_T],
                               [NOCOV, NOCOV, NOCOV, NOCOV, NOCOV]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))


if __name__ == "__main__":
    ut.main()
