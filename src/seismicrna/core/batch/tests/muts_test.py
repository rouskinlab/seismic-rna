import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.muts import calc_muts_matrix
from seismicrna.core.rel.code import MATCH, NOCOV
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
        read_nums = np.arange(9)
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5]]).T
        seg_end3s = np.array([[1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))

    def test_partial_reads_no_muts_some_masked(self):
        section = Section("myref", DNA("TAGCT"))
        section.mask_gu()
        muts = dict()
        read_nums = np.arange(9)
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5]]).T
        seg_end3s = np.array([[1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        result = calc_muts_matrix(section,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, NOCOV]],
                              read_nums,
                              section.unmasked)
        self.assertTrue(expect.equals(result))


if __name__ == "__main__":
    ut.main()
