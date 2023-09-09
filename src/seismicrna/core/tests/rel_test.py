"""

Tests for the Relate Core Module

========================================================================

"""

import unittest as ut

import numpy as np
import pandas as pd

from ..rel import (IRREC, INDEL, NOCOV, MATCH, DELET, INS_5, INS_3,
                   INS_8, MINS5, MINS3, ANY_8, SUB_A, SUB_C, SUB_G,
                   SUB_T, SUB_N, ANY_B, ANY_D, ANY_H, ANY_V, ANY_N,
                   CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                   CIG_DELET, CIG_INSRT, CIG_SCLIP,
                   MIN_QUAL, MAX_QUAL, encode_relate, encode_match,
                   parse_cigar, blank_relvec)
from ..sect import seq_pos_to_index
from ..seq import DNA


class TestConstants(ut.TestCase):
    """ Test constants of `rel` module. """

    def test_primary_codes(self):
        """ Test the primary relation codes. """
        for exp, code in enumerate([MATCH, DELET, INS_5, INS_3,
                                    SUB_A, SUB_C, SUB_G, SUB_T]):
            self.assertEqual(code, 2 ** exp)

    def test_derived_codes(self):
        """ Test the derived relation codes. """
        self.assertEqual(IRREC, 0)
        self.assertEqual(MINS5, 5)
        self.assertEqual(MINS3, 9)
        self.assertEqual(INS_8, 12)
        self.assertEqual(ANY_8, 13)
        self.assertEqual(INDEL, 14)
        self.assertEqual(SUB_N, 240)
        self.assertEqual(ANY_N, 241)
        self.assertEqual(NOCOV, 255)
        self.assertEqual(ANY_B, ANY_N - SUB_A)
        self.assertEqual(ANY_D, ANY_N - SUB_C)
        self.assertEqual(ANY_H, ANY_N - SUB_G)
        self.assertEqual(ANY_V, ANY_N - SUB_T)


class TestEncodeRelate(ut.TestCase):
    """ Test function `encode_relate`. """

    def test_encode_relate_hi_qual(self):
        """ Test when the quality is at least the minimum. """
        for ref, ref_sub in zip("ACGTN",
                                [ANY_B, ANY_D, ANY_H, ANY_V, ANY_N],
                                strict=True):
            for read, read_sub in zip("ACGTN",
                                      [SUB_A, SUB_C, SUB_G, SUB_T, SUB_N],
                                      strict=True):
                code = encode_relate(ref, read, MAX_QUAL, MAX_QUAL)
                if ref == 'N':
                    self.assertEqual(code, ANY_N)
                elif read == 'N':
                    self.assertEqual(code, ref_sub)
                else:
                    self.assertEqual(code, MATCH if read == ref else read_sub)

    def test_encode_relate_lo_qual(self):
        """ Test when the quality is less than the minimum. """
        for ref, sub in zip("ACGTN", [ANY_B, ANY_D, ANY_H, ANY_V, ANY_N],
                            strict=True):
            for read in "ACGTN":
                self.assertEqual(encode_relate(ref, read, MIN_QUAL, MAX_QUAL),
                                 sub)


class TestEncodeMatch(ut.TestCase):
    """ Test function `encode_match`. """

    def test_encode_match_hi_qual(self):
        """ Test when the quality is at least the minimum. """
        for read in "ACGT":
            self.assertEqual(encode_match(read, MAX_QUAL, MAX_QUAL), MATCH)
        self.assertEqual(encode_match('N', MAX_QUAL, MAX_QUAL), ANY_N)

    def test_encode_match_lo_qual(self):
        """ Test when the quality is less than the minimum. """
        for read, sub in zip("ACGTN", [ANY_B, ANY_D, ANY_H, ANY_V, ANY_N]):
            self.assertEqual(encode_match(read, MIN_QUAL, MAX_QUAL), sub)





class TestParseCigar(ut.TestCase):
    """ Test function `parse_cigar`. """

    def test_cigar_match_subst_valid(self):
        """ Parse a valid CIGAR string with match and subst codes. """
        cigar = "9S23=1X13=1D9=2I56=3S"
        expect = [(CIG_SCLIP, 9), (CIG_MATCH, 23), (CIG_SUBST, 1),
                  (CIG_MATCH, 13), (CIG_DELET, 1), (CIG_MATCH, 9),
                  (CIG_INSRT, 2), (CIG_MATCH, 56), (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)

    def test_cigar_align_valid(self):
        """ Parse a valid CIGAR string with align codes. """
        cigar = "9S37M1D9M2I56M3S"
        expect = [(CIG_SCLIP, 9), (CIG_ALIGN, 37), (CIG_DELET, 1),
                  (CIG_ALIGN, 9), (CIG_INSRT, 2), (CIG_ALIGN, 56),
                  (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)




class TestBlankRelvec(ut.TestCase):
    """ Test function `blank_relvec`. """

    def compare_numpy(self, result, expect: np.ndarray):
        self.assertIsInstance(result, np.ndarray)
        self.assertIs(type(result), type(expect))
        self.assertIs(result.dtype, expect.dtype)
        self.assertTrue(np.array_equal(result, expect))

    def compare_pandas(self, result, expect: pd.Series | pd.DataFrame):
        self.assertIsInstance(result, (pd.Series, pd.DataFrame))
        self.assertIs(type(result), type(expect))
        self.assertTrue(expect.equals(result))

    def test_numpy_1d(self):
        """ Test returning a 1D NumPy array. """
        for length in range(10):
            self.compare_numpy(blank_relvec(length),
                               np.full(shape=(length,),
                                       fill_value=NOCOV,
                                       dtype=np.uint8))

    def test_numpy_2d_int(self):
        """ Test returning a 2D NumPy array with integer reads. """
        for length in range(10):
            for n_reads in range(10):
                self.compare_numpy(blank_relvec(length, n_reads),
                                   np.full(shape=(n_reads, length),
                                           fill_value=NOCOV,
                                           dtype=np.uint8))

    def test_numpy_2d_list(self):
        """ Test returning a 2D NumPy array with a list of reads. """
        for length in range(10):
            for n_reads in range(10):
                for list_func in [list, np.array, pd.Index]:
                    reads = list_func(list(map(str, range(n_reads))))
                    self.compare_numpy(blank_relvec(length, reads),
                                       np.full(shape=(n_reads, length),
                                               fill_value=NOCOV,
                                               dtype=np.uint8))

    def test_pandas_series(self):
        """ Test returning a Pandas Series. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            self.compare_pandas(blank_relvec(seq),
                                pd.Series(NOCOV, index=index, dtype=np.uint8))

    def test_pandas_dataframe_int(self):
        """ Test returning a Pandas DataFrame with integer reads. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            for n_reads in range(10):
                reads = [f"Read_{i + 1}" for i in range(n_reads)]
                self.compare_pandas(blank_relvec(seq, n_reads),
                                    pd.DataFrame(NOCOV,
                                                 index=reads,
                                                 columns=index,
                                                 dtype=np.uint8))

    def test_pandas_dataframe_list(self):
        """ Test returning a Pandas DataFrame with integer reads. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            for n_reads in range(10):
                for list_func in [list, np.array, pd.Index]:
                    reads = list_func(list(map(str, range(n_reads))))
                    self.compare_pandas(blank_relvec(seq, reads),
                                        pd.DataFrame(NOCOV,
                                                     index=reads,
                                                     columns=index,
                                                     dtype=np.uint8))
