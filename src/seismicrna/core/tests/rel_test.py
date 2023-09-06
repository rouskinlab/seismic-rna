"""

Tests for the Relate Core Module

========================================================================

"""

import string
import unittest as ut
from itertools import chain
from typing import Generator, Sequence

import numpy as np
import pandas as pd

from ..rel import (IRREC, MATCH, DELET,
                   INS_5, INS_3, INS_8, MINS5, MINS3, ANY_8,
                   SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                   ANY_B, ANY_D, ANY_H, ANY_V, ANY_N,
                   INDEL, NOCOV,
                   MIN_QUAL, MAX_QUAL,
                   CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                   CIG_DELET, CIG_INSRT, CIG_SCLIP,
                   CigarOp, parse_cigar, count_cigar_muts, find_cigar_op_pos,
                   translate_relvec, blank_relvec, random_relvecs,
                   encode_relate, encode_match,
                   validate_relvec, iter_relvecs_q53, iter_relvecs_all,
                   relvec_to_read, ref_to_alignments, iter_alignments, as_sam)
from ..sect import seq_pos_to_index
from ..seq import DNA, expand_degenerate_seq, BASEA, BASEC, BASEG, BASET, BASEN


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


class TestCigarOp(ut.TestCase):
    """ Test class `CigarOp`. """

    ops = CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_DELET, CIG_INSRT, CIG_SCLIP

    def test_cigar_init_valid(self):
        """ Test that CigarOp instances can be created. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            self.assertIsInstance(cigar_op, CigarOp)
            self.assertEqual(cigar_op.op, op)
            self.assertEqual(cigar_op._len, 1)

    def test_cigar_init_invalid(self):
        """ Test that CigarOp instances cannot be created with invalid
        operations. """
        for op in string.printable:
            if op not in self.ops:
                self.assertRaisesRegex(ValueError, f"Invalid CIGAR operation: ",
                                       CigarOp, op)

    def test_cigar_lengthen(self):
        """ Test lengthening CIGAR operations. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            for length in range(1, 10):
                self.assertEqual(cigar_op._len, length)
                cigar_op.lengthen()

    def test_cigar_str(self):
        """ Test string representations of CIGAR operations. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            for length in range(1, 10):
                self.assertEqual(str(cigar_op), f"{length}{op}")
                cigar_op.lengthen()


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


class TestCountCigarMuts(ut.TestCase):
    """ Test function `count_cigar_muts`. """

    def test_cigar_match_subst_valid(self):
        """ Count mutations in a valid CIGAR string. """
        self.assertEqual(count_cigar_muts("9S23=1X13=1D9=2I56=3S"), 4)


class TestFindCigarOpPos(ut.TestCase):
    """ Test function `find_cigar_op_pos`. """

    def test_cigar_xeq_aln_valid(self):
        """ Find aligned positions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_ALIGN)),
                         [])

    def test_cigar_m_aln_valid(self):
        """ Find aligned positions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_ALIGN)),
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                          15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 28, 29,
                          30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                          42, 43, 44, 45, 46, 47, 48, 49, 52, 53, 54, 55,
                          56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                          68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                          80, 81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92,
                          93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                          104, 105, 106, 107, 108])

    def test_cigar_xeq_mat_valid(self):
        """ Find matches in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_MATCH)),
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                          15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30,
                          31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
                          43, 44, 45, 46, 47, 48, 49, 52, 53, 54, 55, 56,
                          57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                          69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
                          81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93,
                          94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
                          105, 106, 107, 108])

    def test_cigar_m_mat_valid(self):
        """ Find matches in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_MATCH)),
                         [])

    def test_cigar_xeq_sub_valid(self):
        """ Find substitutions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_SUBST)),
                         [24])

    def test_cigar_m_sub_valid(self):
        """ Find substitutions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_SUBST)),
                         [])

    def test_cigar_xeq_del_valid(self):
        """ Find deletions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_DELET)),
                         [])

    def test_cigar_m_del_valid(self):
        """ Find deletions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_DELET)),
                         [])

    def test_cigar_xeq_ins_valid(self):
        """ Find insertions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_INSRT)),
                         [25, 26, 27, 50, 51, 83])

    def test_cigar_m_ins_valid(self):
        """ Find insertions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_INSRT)),
                         [25, 26, 27, 50, 51, 83])

    def test_cigar_xeq_scl_valid(self):
        """ Find soft clippings in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_SCLIP)),
                         [])

    def test_cigar_m_scl_valid(self):
        """ Find soft clippings in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_SCLIP)),
                         [])


class TestValidateRelVec(ut.TestCase):
    """ Test function `validate_relvecs`. """

    def test_matches(self):
        """ Test that a relation vector of all matches is valid. """
        self.assertIsNotNone(validate_relvec(np.ones(8, np.uint8)))

    def test_nocov_margin(self):
        """ Test that a relation vector with non-covered positions on
        its margins is valid. """
        relvec = np.full(8, NOCOV, np.uint8)
        relvec[4] = 1
        self.assertIsNotNone(validate_relvec(relvec))

    def test_nocov_middle(self):
        """ Test that a relation vector with non-covered positions in
        the middle is invalid. """
        relvec = np.full(8, MATCH, np.uint8)
        relvec[4] = NOCOV
        self.assertRaises(ValueError, validate_relvec, relvec)

    def test_blank(self):
        """ Test that an entirely blank relation vector is invalid. """
        self.assertRaises(ValueError, validate_relvec,
                          np.full(8, NOCOV, np.uint8))

    def test_not_array(self):
        """ Test that a non-array relation vector is invalid. """
        self.assertRaises(TypeError, validate_relvec,
                          np.ones(8, np.uint8).tolist())

    def test_not_uint8(self):
        """ Test that a non-uint8 relation vector is invalid. """
        self.assertRaises(TypeError, validate_relvec, np.ones(8, np.int8))


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


class TestIterRelvecsQ53(ut.TestCase):
    """ Test function `iter_relvecs_q53`. """

    @staticmethod
    def list_rels(seq: str, low_qual: Sequence[int] = (),
                  end5: int | None = None, end3: int | None = None):
        """ Convenience function to run `rel.iter_relvecs_q53` from a
        sequence of str and return a list of lists of ints. """
        return list(map(np.ndarray.tolist,
                        iter_relvecs_q53(DNA(seq), low_qual, end5, end3)))

    def test_type(self):
        """ Test that the result is a Generator of NumPy arrays. """
        self.assertTrue(isinstance(iter_relvecs_q53(DNA('A')), Generator))
        self.assertTrue(all(isinstance(relvec, np.ndarray)
                            for relvec in iter_relvecs_q53(DNA('A'))))
        self.assertIs(list(iter_relvecs_q53(DNA('A')))[0].dtype.type, np.uint8)

    def test_a(self):
        """ Test with sequence 'A'. """
        self.assertEqual(self.list_rels('A'),
                         [[MATCH], [SUB_C], [SUB_G], [SUB_T]])

    def test_c(self):
        """ Test with sequence 'C'. """
        self.assertEqual(self.list_rels('C'),
                         [[SUB_A], [MATCH], [SUB_G], [SUB_T]])

    def test_g(self):
        """ Test with sequence 'G'. """
        self.assertEqual(self.list_rels('G'),
                         [[SUB_A], [SUB_C], [MATCH], [SUB_T]])

    def test_t(self):
        """ Test with sequence 'T'. """
        self.assertEqual(self.list_rels('T'),
                         [[SUB_A], [SUB_C], [SUB_G], [MATCH]])

    def test_n(self):
        """ Test with sequence 'N'. """
        self.assertEqual(self.list_rels('N'),
                         [[ANY_N], [ANY_N], [ANY_N], [ANY_N]])

    def test_aa(self):
        """ Test with sequence 'AA'. """
        self.assertEqual(self.list_rels("AA"),
                         [[MATCH, MATCH],
                          [MATCH + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_G],
                          [MATCH + INS_5, INS_3 + SUB_G],
                          [MATCH, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_C, MATCH],
                          [SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_G],
                          [SUB_C + INS_5, INS_3 + SUB_G],
                          [SUB_C, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_G, MATCH],
                          [SUB_G + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_G],
                          [SUB_G + INS_5, INS_3 + SUB_G],
                          [SUB_G, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_T],
                          [SUB_T, MATCH],
                          [SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_G],
                          [SUB_T + INS_5, INS_3 + SUB_G],
                          [SUB_T, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_T]])

    def test_low_qual(self):
        """ Test with each low-quality base. """
        for base, low_qual in zip("ACGTN", [ANY_B, ANY_D, ANY_H, ANY_V, ANY_N],
                                  strict=True):
            self.assertEqual(self.list_rels(base, [1]),
                             [[low_qual]])

    def test_low_qual_invalid(self):
        """ Test that invalid low-qual positions raise ValueError. """
        seq = "ACGTN"
        for n in range(1, len(seq) + 1):
            self.assertTrue(isinstance(self.list_rels(seq[: n], [1]), list))
            self.assertTrue(isinstance(self.list_rels(seq[: n], [n]), list))
            self.assertRaises(ValueError, self.list_rels, seq[: n], [0])
            self.assertRaises(ValueError, self.list_rels, seq[: n], [n + 1])

    def test_xaax(self):
        """ Test that bases with no coverage are marked. """
        self.assertEqual(self.list_rels("TAAG", end5=2, end3=3),
                         [[NOCOV, MATCH, MATCH, NOCOV],
                          [NOCOV, MATCH + INS_5, INS_3 + MATCH, NOCOV],
                          [NOCOV, MATCH, SUB_C, NOCOV],
                          [NOCOV, MATCH + INS_5, INS_3 + SUB_C, NOCOV],
                          [NOCOV, MATCH, SUB_G, NOCOV],
                          [NOCOV, MATCH + INS_5, INS_3 + SUB_G, NOCOV],
                          [NOCOV, MATCH, SUB_T, NOCOV],
                          [NOCOV, MATCH + INS_5, INS_3 + SUB_T, NOCOV],
                          [NOCOV, SUB_C, MATCH, NOCOV],
                          [NOCOV, SUB_C + INS_5, INS_3 + MATCH, NOCOV],
                          [NOCOV, SUB_C, SUB_C, NOCOV],
                          [NOCOV, SUB_C + INS_5, INS_3 + SUB_C, NOCOV],
                          [NOCOV, SUB_C, SUB_G, NOCOV],
                          [NOCOV, SUB_C + INS_5, INS_3 + SUB_G, NOCOV],
                          [NOCOV, SUB_C, SUB_T, NOCOV],
                          [NOCOV, SUB_C + INS_5, INS_3 + SUB_T, NOCOV],
                          [NOCOV, SUB_G, MATCH, NOCOV],
                          [NOCOV, SUB_G + INS_5, INS_3 + MATCH, NOCOV],
                          [NOCOV, SUB_G, SUB_C, NOCOV],
                          [NOCOV, SUB_G + INS_5, INS_3 + SUB_C, NOCOV],
                          [NOCOV, SUB_G, SUB_G, NOCOV],
                          [NOCOV, SUB_G + INS_5, INS_3 + SUB_G, NOCOV],
                          [NOCOV, SUB_G, SUB_T, NOCOV],
                          [NOCOV, SUB_G + INS_5, INS_3 + SUB_T, NOCOV],
                          [NOCOV, SUB_T, MATCH, NOCOV],
                          [NOCOV, SUB_T + INS_5, INS_3 + MATCH, NOCOV],
                          [NOCOV, SUB_T, SUB_C, NOCOV],
                          [NOCOV, SUB_T + INS_5, INS_3 + SUB_C, NOCOV],
                          [NOCOV, SUB_T, SUB_G, NOCOV],
                          [NOCOV, SUB_T + INS_5, INS_3 + SUB_G, NOCOV],
                          [NOCOV, SUB_T, SUB_T, NOCOV],
                          [NOCOV, SUB_T + INS_5, INS_3 + SUB_T, NOCOV]])

    def test_agg(self):
        """ Test with sequence 'AGG'. """
        rels = self.list_rels("AGG")
        self.assertEqual(rels,
                         [[MATCH, SUB_A, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_A, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_A + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_A, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_A, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_A + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_A, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_A, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_A + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_A + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_A, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_A, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_T],
                          [MATCH, SUB_A + INS_5, INS_3 + SUB_T],
                          [MATCH, SUB_C, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_C, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_C + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_C, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_C, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_C + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_C, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_C, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_C + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_C + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_C, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_C, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_T],
                          [MATCH, SUB_C + INS_5, INS_3 + SUB_T],
                          [MATCH, MATCH, SUB_A],
                          [MATCH + INS_5, INS_3 + MATCH, SUB_A],
                          [MATCH + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_A],
                          [MATCH, MATCH + INS_5, INS_3 + SUB_A],
                          [MATCH, MATCH, SUB_C],
                          [MATCH + INS_5, INS_3 + MATCH, SUB_C],
                          [MATCH + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_C],
                          [MATCH, MATCH + INS_5, INS_3 + SUB_C],
                          [MATCH, MATCH, MATCH],
                          [MATCH + INS_5, INS_3 + MATCH, MATCH],
                          [MATCH + INS_5, INS_3 + MATCH + INS_5, INS_3 + MATCH],
                          [MATCH, MATCH + INS_5, INS_3 + MATCH],
                          [MATCH, MATCH, SUB_T],
                          [MATCH + INS_5, INS_3 + MATCH, SUB_T],
                          [MATCH + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_T],
                          [MATCH, MATCH + INS_5, INS_3 + SUB_T],
                          [MATCH, SUB_T, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_T, SUB_A],
                          [MATCH + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_T + INS_5, INS_3 + SUB_A],
                          [MATCH, SUB_T, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_T, SUB_C],
                          [MATCH + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_T + INS_5, INS_3 + SUB_C],
                          [MATCH, SUB_T, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_T, MATCH],
                          [MATCH + INS_5, INS_3 + SUB_T + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_T + INS_5, INS_3 + MATCH],
                          [MATCH, SUB_T, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_T, SUB_T],
                          [MATCH + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_T],
                          [MATCH, SUB_T + INS_5, INS_3 + SUB_T],
                          [MATCH, DELET, SUB_A],
                          [MATCH, DELET, SUB_C],
                          [MATCH, DELET, MATCH],
                          [MATCH, DELET, SUB_T],
                          [SUB_C, SUB_A, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_A, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_A, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_A, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_A, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_A, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_A, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_A, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_C, SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_C, SUB_C, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_C, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_C, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_C, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_C, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_C, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_C, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_C, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_C, SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_C, MATCH, SUB_A],
                          [SUB_C + INS_5, INS_3 + MATCH, SUB_A],
                          [SUB_C + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_C, MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_C, MATCH, SUB_C],
                          [SUB_C + INS_5, INS_3 + MATCH, SUB_C],
                          [SUB_C + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_C, MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_C, MATCH, MATCH],
                          [SUB_C + INS_5, INS_3 + MATCH, MATCH],
                          [SUB_C + INS_5, INS_3 + MATCH + INS_5, INS_3 + MATCH],
                          [SUB_C, MATCH + INS_5, INS_3 + MATCH],
                          [SUB_C, MATCH, SUB_T],
                          [SUB_C + INS_5, INS_3 + MATCH, SUB_T],
                          [SUB_C + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_C, MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_C, SUB_T, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_T, SUB_A],
                          [SUB_C + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_C, SUB_T, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_T, SUB_C],
                          [SUB_C + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_C, SUB_T, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_T, MATCH],
                          [SUB_C + INS_5, INS_3 + SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_C, SUB_T, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_T, SUB_T],
                          [SUB_C + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_C, SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_C, DELET, SUB_A],
                          [SUB_C, DELET, SUB_C],
                          [SUB_C, DELET, MATCH],
                          [SUB_C, DELET, SUB_T],
                          [SUB_G, SUB_A, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_A, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_A, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_A, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_A, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_A, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_A, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_A, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_G, SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_G, SUB_C, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_C, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_C, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_C, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_C, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_C, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_C, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_C, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_G, SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_G, MATCH, SUB_A],
                          [SUB_G + INS_5, INS_3 + MATCH, SUB_A],
                          [SUB_G + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_G, MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_G, MATCH, SUB_C],
                          [SUB_G + INS_5, INS_3 + MATCH, SUB_C],
                          [SUB_G + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_G, MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_G, MATCH, MATCH],
                          [SUB_G + INS_5, INS_3 + MATCH, MATCH],
                          [SUB_G + INS_5, INS_3 + MATCH + INS_5, INS_3 + MATCH],
                          [SUB_G, MATCH + INS_5, INS_3 + MATCH],
                          [SUB_G, MATCH, SUB_T],
                          [SUB_G + INS_5, INS_3 + MATCH, SUB_T],
                          [SUB_G + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_G, MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_G, SUB_T, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_T, SUB_A],
                          [SUB_G + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_G, SUB_T, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_T, SUB_C],
                          [SUB_G + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_G, SUB_T, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_T, MATCH],
                          [SUB_G + INS_5, INS_3 + SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_G, SUB_T, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_T, SUB_T],
                          [SUB_G + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_G, SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_G, DELET, SUB_A],
                          [SUB_G, DELET, SUB_C],
                          [SUB_G, DELET, MATCH],
                          [SUB_G, DELET, SUB_T],
                          [SUB_T, SUB_A, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_A, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_A + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_A, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_A, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_A + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_A, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_A, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_A + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_A, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_A, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_T, SUB_A + INS_5, INS_3 + SUB_T],
                          [SUB_T, SUB_C, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_C, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_C + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_C, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_C, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_C + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_C, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_C, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_C + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_C, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_C, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_T, SUB_C + INS_5, INS_3 + SUB_T],
                          [SUB_T, MATCH, SUB_A],
                          [SUB_T + INS_5, INS_3 + MATCH, SUB_A],
                          [SUB_T + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_T, MATCH + INS_5, INS_3 + SUB_A],
                          [SUB_T, MATCH, SUB_C],
                          [SUB_T + INS_5, INS_3 + MATCH, SUB_C],
                          [SUB_T + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_T, MATCH + INS_5, INS_3 + SUB_C],
                          [SUB_T, MATCH, MATCH],
                          [SUB_T + INS_5, INS_3 + MATCH, MATCH],
                          [SUB_T + INS_5, INS_3 + MATCH + INS_5, INS_3 + MATCH],
                          [SUB_T, MATCH + INS_5, INS_3 + MATCH],
                          [SUB_T, MATCH, SUB_T],
                          [SUB_T + INS_5, INS_3 + MATCH, SUB_T],
                          [SUB_T + INS_5, INS_3 + MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_T, MATCH + INS_5, INS_3 + SUB_T],
                          [SUB_T, SUB_T, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_T, SUB_A],
                          [SUB_T + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_T + INS_5, INS_3 + SUB_A],
                          [SUB_T, SUB_T, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_T, SUB_C],
                          [SUB_T + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_T + INS_5, INS_3 + SUB_C],
                          [SUB_T, SUB_T, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_T, MATCH],
                          [SUB_T + INS_5, INS_3 + SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_T + INS_5, INS_3 + MATCH],
                          [SUB_T, SUB_T, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_T, SUB_T],
                          [SUB_T + INS_5, INS_3 + SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_T, SUB_T + INS_5, INS_3 + SUB_T],
                          [SUB_T, DELET, SUB_A],
                          [SUB_T, DELET, SUB_C],
                          [SUB_T, DELET, MATCH],
                          [SUB_T, DELET, SUB_T]])


class TestIterRelvecsAll(ut.TestCase):
    """ Test function `iter_relvecs_all`. """

    def assert_equal(self, ref: DNA, expects: list):
        """ Check that the expected and actual results match. """
        for exp, res in zip(chain(*expects),
                            iter_relvecs_all(ref),
                            strict=True):
            with self.subTest(exp=exp, res=res):
                self.assertTrue(np.all(exp == res))

    def test_length_1(self):
        """ Test with all length-1 DNA sequences. """
        for ref in expand_degenerate_seq(DNA("N")):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
            ]
            self.assert_equal(ref, expects)

    def test_length_2(self):
        """ Test with all length-2 DNA sequences. """
        for ref in expand_degenerate_seq(DNA("NN")):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
                iter_relvecs_q53(ref, [], 1, 2),
                iter_relvecs_q53(ref, [1], 1, 2),
                iter_relvecs_q53(ref, [2], 1, 2),
                iter_relvecs_q53(ref, [1, 2], 1, 2),
                iter_relvecs_q53(ref, [], 2, 2),
                iter_relvecs_q53(ref, [2], 2, 2),
            ]
            self.assert_equal(ref, expects)

    def test_length_3(self):
        """ Test with all length-3 DNA sequences. """
        for ref in expand_degenerate_seq(DNA("NNN")):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
                iter_relvecs_q53(ref, [], 1, 2),
                iter_relvecs_q53(ref, [1], 1, 2),
                iter_relvecs_q53(ref, [2], 1, 2),
                iter_relvecs_q53(ref, [1, 2], 1, 2),
                iter_relvecs_q53(ref, [], 1, 3),
                iter_relvecs_q53(ref, [1], 1, 3),
                iter_relvecs_q53(ref, [2], 1, 3),
                iter_relvecs_q53(ref, [3], 1, 3),
                iter_relvecs_q53(ref, [1, 2], 1, 3),
                iter_relvecs_q53(ref, [1, 3], 1, 3),
                iter_relvecs_q53(ref, [2, 3], 1, 3),
                iter_relvecs_q53(ref, [1, 2, 3], 1, 3),
                iter_relvecs_q53(ref, [], 2, 2),
                iter_relvecs_q53(ref, [2], 2, 2),
                iter_relvecs_q53(ref, [], 2, 3),
                iter_relvecs_q53(ref, [2], 2, 3),
                iter_relvecs_q53(ref, [3], 2, 3),
                iter_relvecs_q53(ref, [2, 3], 2, 3),
                iter_relvecs_q53(ref, [], 3, 3),
                iter_relvecs_q53(ref, [3], 3, 3),
            ]
            self.assert_equal(ref, expects)


class TestRelvecToRead(ut.TestCase):
    """ Test function `relvec_to_read`. """

    def assert_equal(self, ref: DNA,
                     relvecs: list[list[int]],
                     expects: list[tuple[str, str, str, int, int]]):
        """ Assert that the actual and expected outputs match. """
        for relvec, expect in zip(relvecs, expects, strict=True):
            with self.subTest(relvec=relvec, expect=expect):
                self.assertEqual(relvec_to_read(ref, np.array(relvec,
                                                              dtype=np.uint8),
                                                MAX_QUAL, MIN_QUAL),
                                 (DNA(expect[0]),) + expect[1:])

    def assert_raise(self, ref: DNA,
                     relvecs: list[list[int]],
                     error: type[Exception],
                     regex: str):
        """ Assert that the relation vectors raise an exception. """
        for relvec in relvecs:
            with self.subTest(relvec=relvec):
                self.assertRaisesRegex(error, regex, relvec_to_read,
                                       ref, np.array(relvec,
                                                     dtype=np.uint8),
                                       MAX_QUAL, MIN_QUAL)

    def test_all_match(self):
        """ Test when the read has four matching bases. """
        ref = DNA("ACGT")
        relvecs = [[MATCH, MATCH, MATCH, MATCH]]
        expects = [("ACGT", "IIII", "4=", 1, 4)]
        self.assert_equal(ref, relvecs, expects)

    def test_all_match_n(self):
        """ Test when the read has four matching bases and an ambiguous
        base. """
        ref = DNA("ACNGT")
        relvecs = [[MATCH, MATCH, ANY_N, MATCH, MATCH]]
        expects = [("ACNGT", "II!II", "2=1M2=", 1, 5)]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_valid(self):
        """ Test when the read does not cover one or both ends of the
        reference. """
        ref = DNA("ACGT")
        relvecs = [
            [NOCOV, MATCH, MATCH, MATCH],
            [MATCH, MATCH, MATCH, NOCOV],
            [NOCOV, MATCH, MATCH, NOCOV],
            [NOCOV, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, NOCOV],
        ]
        expects = [
            ("CGT", "III", "3=", 2, 4),
            ("ACG", "III", "3=", 1, 3),
            ("CG", "II", "2=", 2, 3),
            ("GT", "II", "2=", 3, 4),
            ("AC", "II", "2=", 1, 2),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_middle_invalid(self):
        """ Test when the read does not cover a middle position. """
        ref = DNA("ACGT")
        relvecs = [
            [MATCH, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, MATCH],
            [MATCH, NOCOV, NOCOV, MATCH],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Expected [0-9]+ base calls")

    def test_nocov_all_invalid(self):
        """ Test when the read does not cover any positions. """
        ref = DNA("ACGT")
        relvecs = [[NOCOV, NOCOV, NOCOV, NOCOV]]
        self.assert_raise(ref, relvecs, ValueError,
                          "Relation vector is blank")

    def test_low_qual_valid(self):
        """ Test when the read has a low-quality base. """
        ref = DNA("ACGT")
        relvecs = [
            [ANY_N - SUB_A, MATCH, MATCH, MATCH],
            [MATCH, ANY_N - SUB_C, MATCH, MATCH],
            [MATCH, MATCH, ANY_N - SUB_G, MATCH],
            [MATCH, MATCH, MATCH, ANY_N - SUB_T],
        ]
        expects = [
            ("NCGT", "!III", "1M3=", 1, 4),
            ("ANGT", "I!II", "1=1M2=", 1, 4),
            ("ACNT", "II!I", "2=1M1=", 1, 4),
            ("ACGN", "III!", "3=1M", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_valid(self):
        """ Test when the read has a substitution. """
        ref = DNA("ACGT")
        relvecs = [
            [SUB_C, MATCH, MATCH, MATCH],
            [SUB_G, MATCH, MATCH, MATCH],
            [SUB_T, MATCH, MATCH, MATCH],
            [MATCH, SUB_A, MATCH, MATCH],
            [MATCH, SUB_G, MATCH, MATCH],
            [MATCH, SUB_T, MATCH, MATCH],
            [MATCH, MATCH, SUB_A, MATCH],
            [MATCH, MATCH, SUB_C, MATCH],
            [MATCH, MATCH, SUB_T, MATCH],
            [MATCH, MATCH, MATCH, SUB_A],
            [MATCH, MATCH, MATCH, SUB_C],
            [MATCH, MATCH, MATCH, SUB_G],
        ]
        expects = [
            ("CCGT", "IIII", "1X3=", 1, 4),
            ("GCGT", "IIII", "1X3=", 1, 4),
            ("TCGT", "IIII", "1X3=", 1, 4),
            ("AAGT", "IIII", "1=1X2=", 1, 4),
            ("AGGT", "IIII", "1=1X2=", 1, 4),
            ("ATGT", "IIII", "1=1X2=", 1, 4),
            ("ACAT", "IIII", "2=1X1=", 1, 4),
            ("ACCT", "IIII", "2=1X1=", 1, 4),
            ("ACTT", "IIII", "2=1X1=", 1, 4),
            ("ACGA", "IIII", "3=1X", 1, 4),
            ("ACGC", "IIII", "3=1X", 1, 4),
            ("ACGG", "IIII", "3=1X", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_invalid(self):
        """ Test when the read has an invalid substitution. """
        ref = DNA("ACGT")
        relvecs = [
            [SUB_A, MATCH, MATCH, MATCH],
            [MATCH, SUB_C, MATCH, MATCH],
            [MATCH, MATCH, SUB_G, MATCH],
            [MATCH, MATCH, MATCH, SUB_T],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Cannot substitute [ACGT] to itself")

    def test_delete_valid(self):
        """ Test when the read has deletions. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 deletion
            [MATCH, DELET, MATCH, MATCH],
            [MATCH, MATCH, DELET, MATCH],
            # 2 deletions
            [MATCH, DELET, DELET, MATCH],
        ]
        expects = [
            # 1 deletion
            ("AGT", "III", "1=1D2=", 1, 4),
            ("ACT", "III", "2=1D1=", 1, 4),
            # 2 deletions
            ("AT", "II", "1=2D1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_delete_invalid(self):
        """ Test when the read has a deletion at either end. """
        ref = DNA("ACGT")
        relvecs = [
            [DELET, MATCH, MATCH, MATCH],
            [NOCOV, DELET, MATCH, MATCH],
            [NOCOV, NOCOV, DELET, MATCH],
            [NOCOV, NOCOV, NOCOV, DELET],
            [DELET, DELET, DELET, DELET],
            [DELET, MATCH, MATCH, DELET],
            [MATCH, MATCH, MATCH, DELET],
            [MATCH, MATCH, DELET, NOCOV],
            [MATCH, DELET, NOCOV, NOCOV],
            [DELET, NOCOV, NOCOV, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Deletion cannot be at position [0-9]+ in .+")

    def test_insert_valid(self):
        """ Test when the read has insertions. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion
            [MINS5, MINS3, MATCH, MATCH],
            [MATCH, MINS5, MINS3, MATCH],
            [MATCH, MATCH, MINS5, MINS3],
            [MINS5, MINS3, MATCH, NOCOV],
            [MATCH, MINS5, MINS3, NOCOV],
            [NOCOV, MINS5, MINS3, MATCH],
            [NOCOV, MATCH, MINS5, MINS3],
            [NOCOV, MINS5, MINS3, NOCOV],
            # 2 insertions, 1 base apart
            [MINS5, ANY_8, MINS3, MATCH],
            [MATCH, MINS5, ANY_8, MINS3],
            # 2 insertions, 2 bases apart
            [MINS5, MINS3, MINS5, MINS3],
            # 3 insertions, 1 base apart
            [MINS5, ANY_8, ANY_8, MINS3],
        ]
        expects = [
            # 1 insertion
            ("ANCGT", "IIIII", "1=1I3=", 1, 4),
            ("ACNGT", "IIIII", "2=1I2=", 1, 4),
            ("ACGNT", "IIIII", "3=1I1=", 1, 4),
            ("ANCG", "IIII", "1=1I2=", 1, 3),
            ("ACNG", "IIII", "2=1I1=", 1, 3),
            ("CNGT", "IIII", "1=1I2=", 2, 4),
            ("CGNT", "IIII", "2=1I1=", 2, 4),
            ("CNG", "III", "1=1I1=", 2, 3),
            # 2 insertions, 1 base apart
            ("ANCNGT", "IIIIII", "1=1I1=1I2=", 1, 4),
            ("ACNGNT", "IIIIII", "2=1I1=1I1=", 1, 4),
            # 2 insertions, 2 bases apart
            ("ANCGNT", "IIIIII", "1=1I2=1I1=", 1, 4),
            # 3 insertions, 1 base apart
            ("ANCNGNT", "IIIIIII", "1=1I1=1I1=1I1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert_end5_invalid(self):
        """ Test when the read has an insertion at the 5' end. """
        ref = DNA("ACGT")
        relvecs = [
            [MATCH, MATCH, MATCH, MINS5],
            [MATCH, MATCH, MINS5, ANY_8],
            [MATCH, MINS5, ANY_8, ANY_8],
            [MINS5, ANY_8, ANY_8, ANY_8],
            [NOCOV, MATCH, MINS5, NOCOV],
            [NOCOV, MINS5, ANY_8, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Position [0-9]+ in .+ cannot be 5' of an insertion")

    def test_insert_end3_invalid(self):
        """ Test when the read has an insertion at the 3' end. """
        ref = DNA("ACGT")
        relvecs = [
            [MINS3, MATCH, MATCH, MATCH],
            [ANY_8, MINS3, MATCH, MATCH],
            [ANY_8, ANY_8, MINS3, MATCH],
            [ANY_8, ANY_8, ANY_8, MINS3],
            [NOCOV, MINS3, MATCH, NOCOV],
            [NOCOV, ANY_8, MINS3, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Position [0-9]+ in .+ cannot be 3' of an insertion")

    def test_insert_dangling_5_invalid(self):
        """ Test when the read has an unmatched 5' insertion. """
        ref = DNA("ACGT")
        relvecs = [
            [MINS5, MATCH, MATCH, MATCH],
            [MATCH, MINS5, MATCH, MATCH],
            [MATCH, MATCH, MINS5, MATCH],
            [MINS5, MATCH, MINS3, MATCH],
            [MATCH, MINS5, MATCH, MINS3],
            [NOCOV, MINS5, MATCH, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Missing 3' ins at [0-9]+ in .+")

    def test_insert_dangling_3_invalid(self):
        """ Test when the read has an unmatched 3' insertion. """
        ref = DNA("ACGT")
        relvecs = [
            [MATCH, MINS3, MATCH, MATCH],
            [MATCH, MATCH, MINS3, MATCH],
            [MATCH, MATCH, MATCH, MINS3],
            [NOCOV, MATCH, MINS3, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Unexpected 3' ins at [0-9+] in .+")

    def test_insert_bare_invalid(self):
        """ Test when the read has bare insertions (with no underlying
        relationship). """
        ref = DNA("ACGT")
        relvecs = [
            # 1 bare insertion
            [MINS5, INS_3, MATCH, MATCH],
            [MATCH, MINS5, INS_3, MATCH],
            [MATCH, MATCH, MINS5, INS_3],
            [INS_5, MINS3, MATCH, MATCH],
            [MATCH, INS_5, MINS3, MATCH],
            [MATCH, MATCH, INS_5, MINS3],
            [NOCOV, MINS5, INS_3, NOCOV],
            [NOCOV, INS_5, MINS3, NOCOV],
            # 2 bare insertions
            [MINS5, INS_8, MINS3, MATCH],
            [MATCH, MINS5, INS_8, MINS3],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Invalid relation 0")

    def test_insert_deletion_invalid(self):
        """ Test when the read has an insertion next to a deletion. """
        ref = DNA("ACGT")
        relvecs = [
            [MINS5, INS_3 + DELET, MATCH, MATCH],
            [MATCH, MINS5, INS_3 + DELET, MATCH],
            [MATCH, MATCH, MINS5, INS_3 + DELET],
            [DELET + INS_5, MINS3, MATCH, MATCH],
            [MATCH, DELET + INS_5, MINS3, MATCH],
            [MATCH, MATCH, DELET + INS_5, MINS3],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Relation .+ is del and ins")

    def test_insert_non_match_valid(self):
        """ Test when the read has insertions next to substitutions or
        low-quality base calls. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion next to 1 substitution
            [SUB_C + INS_5, MINS3, MATCH, MATCH],
            [SUB_C, MINS5, MINS3, MATCH],
            [MINS5, INS_3 + SUB_T, MATCH, MATCH],
            [MATCH, SUB_T + INS_5, MINS3, MATCH],
            [MATCH, SUB_T, MINS5, MINS3],
            # 1 insertion next to 1 low-quality base call
            [ANY_N - SUB_A + INS_5, MINS3, MATCH, MATCH],
            [ANY_N - SUB_A, MINS5, MINS3, MATCH],
            [MINS5, INS_3 + ANY_N - SUB_C, MATCH, MATCH],
            [MATCH, ANY_N - SUB_C + INS_5, MINS3, MATCH],
            [MATCH, ANY_N - SUB_C, MINS5, MINS3],
        ]
        expects = [
            # 1 insertion next to 1 substitution
            ("CNCGT", "IIIII", "1X1I3=", 1, 4),
            ("CCNGT", "IIIII", "1X1=1I2=", 1, 4),
            ("ANTGT", "IIIII", "1=1I1X2=", 1, 4),
            ("ATNGT", "IIIII", "1=1X1I2=", 1, 4),
            ("ATGNT", "IIIII", "1=1X1=1I1=", 1, 4),
            # 1 insertion next to 1 low-quality base call
            ("NNCGT", "!IIII", "1M1I3=", 1, 4),
            ("NCNGT", "!IIII", "1M1=1I2=", 1, 4),
            ("ANNGT", "II!II", "1=1I1M2=", 1, 4),
            ("ANNGT", "I!III", "1=1M1I2=", 1, 4),
            ("ANGNT", "I!III", "1=1M1=1I1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)


class TestAsSam(ut.TestCase):
    """ Test function `as_sam`. """

    def test_line_in_sam_format(self):
        line = as_sam("FS10000136:77:BPG61616-2310:1:1101:1000:1300", 99,
                      "SARS2_FSE", 1, 42, "151=", "=", 133, 283,
                      DNA("CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTG"
                          "GAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTC"
                          "AGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC"),
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:")
        expect = ("FS10000136:77:BPG61616-2310:1:1101:1000:1300\t99\tSARS2_FSE"
                  "\t1\t42\t151=\t=\t133\t283\t"
                  "CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTT"
                  "ATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCG"
                  "TTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC\t"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:\n")
        self.assertEqual(line, expect)
