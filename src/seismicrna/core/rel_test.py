"""
Core -- Relation Testing Module
========================================================================
Auth: Matty

Unit tests for `core.rel`.
"""

import unittest as ut
from itertools import chain
from typing import Generator, Sequence

import numpy as np

from .rel import (IRREC, MATCH, DELET,
                  INS_5, INS_3, INS_8, MINS5, MINS3, ANY_8,
                  SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                  ANY_B, ANY_D, ANY_H, ANY_V, ANY_N,
                  INDEL, NOCOV,
                  MIN_QUAL, MAX_QUAL,
                  CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_DELET, CIG_INSRT, CIG_SCLIP,
                  parse_cigar, count_cigar_muts, find_cigar_op_pos,
                  validate_relvec, iter_relvecs_q53, iter_relvecs_all,
                  relvec_to_read, as_sam)
from .seq import DNA, expand_degenerate_seq


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


class TestParseCigar(ut.TestCase):
    """ Test function `parse_cigar`. """

    def test_cigar_match_subst_valid(self):
        """ Parse a valid CIGAR string with match and subst codes. """
        cigar = b"9S23=1X13=1D9=2I56=3S"
        expect = [(CIG_SCLIP, 9), (CIG_MATCH, 23), (CIG_SUBST, 1),
                  (CIG_MATCH, 13), (CIG_DELET, 1), (CIG_MATCH, 9),
                  (CIG_INSRT, 2), (CIG_MATCH, 56), (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)

    def test_cigar_align_valid(self):
        """ Parse a valid CIGAR string with align codes. """
        cigar = b"9S37M1D9M2I56M3S"
        expect = [(CIG_SCLIP, 9), (CIG_ALIGN, 37), (CIG_DELET, 1),
                  (CIG_ALIGN, 9), (CIG_INSRT, 2), (CIG_ALIGN, 56),
                  (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)


class TestCountCigarMuts(ut.TestCase):
    """ Test function `count_cigar_muts`. """

    def test_cigar_match_subst_valid(self):
        """ Count mutations in a valid CIGAR string. """
        self.assertEqual(count_cigar_muts(b"9S23=1X13=1D9=2I56=3S"), 4)


class TestFindCigarOpPos(ut.TestCase):
    """ Test function `find_cigar_op_pos`. """

    def test_cigar_xeq_ins_valid(self):
        """ Find insertions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos(b"9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_INSRT)),
                         [25, 26, 27, 50, 51, 83])


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


class TestIterRelvecsQ53(ut.TestCase):
    """ Test function `iter_relvecs_q53`. """

    @staticmethod
    def list_rels(seq: str, low_qual: Sequence[int] = (),
                  end5: int | None = None, end3: int | None = None):
        """ Convenience function to run `rel.iter_relvecs_q53` from a
        sequence of str and return a list of lists of ints. """
        return list(map(np.ndarray.tolist,
                        iter_relvecs_q53(DNA(seq.encode()),
                                         low_qual, end5, end3)))

    def test_type(self):
        """ Test that the result is a Generator of NumPy arrays. """
        self.assertTrue(isinstance(iter_relvecs_q53(DNA(b"A")), Generator))
        self.assertTrue(all(isinstance(relvec, np.ndarray)
                            for relvec in iter_relvecs_q53(DNA(b"A"))))
        self.assertIs(list(iter_relvecs_q53(DNA(b"A")))[0].dtype.type, np.uint8)

    def test_a(self):
        """ Test with sequence 'A'. """
        self.assertEqual(self.list_rels("A"),
                         [[MATCH], [SUB_C], [SUB_G], [SUB_T]])

    def test_c(self):
        """ Test with sequence 'C'. """
        self.assertEqual(self.list_rels("C"),
                         [[SUB_A], [MATCH], [SUB_G], [SUB_T]])

    def test_g(self):
        """ Test with sequence 'G'. """
        self.assertEqual(self.list_rels("G"),
                         [[SUB_A], [SUB_C], [MATCH], [SUB_T]])

    def test_t(self):
        """ Test with sequence 'T'. """
        self.assertEqual(self.list_rels("T"),
                         [[SUB_A], [SUB_C], [SUB_G], [MATCH]])

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
        for base, low_qual in zip("ACGT", [ANY_B, ANY_D, ANY_H, ANY_V]):
            self.assertEqual(self.list_rels(base, [1]),
                             [[low_qual]])

    def test_low_qual_invalid(self):
        """ Test that invalid low-qual positions raise ValueError. """
        seq = "ACGT"
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
        for ref in expand_degenerate_seq(b"N"):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
            ]
            self.assert_equal(ref, expects)

    def test_length_2(self):
        """ Test with all length-2 DNA sequences. """
        for ref in expand_degenerate_seq(b"NN"):
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
        for ref in expand_degenerate_seq(b"NNN"):
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
                     expects: list[tuple[bytes, bytes, bytes, int, int]]):
        """ Assert that the actual and expected outputs match. """
        for relvec, expect in zip(relvecs, expects, strict=True):
            with self.subTest(relvec=relvec, expect=expect):
                self.assertEqual(relvec_to_read(ref, np.array(relvec,
                                                              dtype=np.uint8),
                                                MAX_QUAL, MIN_QUAL),
                                 expect)

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
        ref = DNA(b"ACGT")
        relvecs = [[MATCH, MATCH, MATCH, MATCH]]
        expects = [(b"ACGT", b"IIII", b"4=", 1, 4)]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_valid(self):
        """ Test when the read does not cover one or both ends of the
        reference. """
        ref = DNA(b"ACGT")
        relvecs = [
            [NOCOV, MATCH, MATCH, MATCH],
            [MATCH, MATCH, MATCH, NOCOV],
            [NOCOV, MATCH, MATCH, NOCOV],
            [NOCOV, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, NOCOV],
        ]
        expects = [
            (b"CGT", b"III", b"3=", 2, 4),
            (b"ACG", b"III", b"3=", 1, 3),
            (b"CG", b"II", b"2=", 2, 3),
            (b"GT", b"II", b"2=", 3, 4),
            (b"AC", b"II", b"2=", 1, 2),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_middle_invalid(self):
        """ Test when the read does not cover a middle position. """
        ref = DNA(b"ACGT")
        relvecs = [
            [MATCH, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, MATCH],
            [MATCH, NOCOV, NOCOV, MATCH],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Expected [0-9]+ base calls")

    def test_nocov_all_invalid(self):
        """ Test when the read does not cover any positions. """
        ref = DNA(b"ACGT")
        relvecs = [[NOCOV, NOCOV, NOCOV, NOCOV]]
        self.assert_raise(ref, relvecs, ValueError,
                          "Relation vector is blank")

    def test_low_qual_valid(self):
        """ Test when the read has a low-quality base. """
        ref = DNA(b"ACGT")
        relvecs = [
            [ANY_N - SUB_A, MATCH, MATCH, MATCH],
            [MATCH, ANY_N - SUB_C, MATCH, MATCH],
            [MATCH, MATCH, ANY_N - SUB_G, MATCH],
            [MATCH, MATCH, MATCH, ANY_N - SUB_T],
        ]
        expects = [
            (b"NCGT", b"!III", b"1M3=", 1, 4),
            (b"ANGT", b"I!II", b"1=1M2=", 1, 4),
            (b"ACNT", b"II!I", b"2=1M1=", 1, 4),
            (b"ACGN", b"III!", b"3=1M", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_low_qual_invalid(self):
        """ Test when the read has an invalid low-quality base. """
        ref = DNA(b"ACGT")
        relvecs = [
            [ANY_N, MATCH, MATCH, MATCH],
            [MATCH, ANY_N, MATCH, MATCH],
            [MATCH, MATCH, ANY_N, MATCH],
            [MATCH, MATCH, MATCH, ANY_N],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          f"Invalid relation {ANY_N}")

    def test_subst_valid(self):
        """ Test when the read has a substitution. """
        ref = DNA(b"ACGT")
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
            (b"CCGT", b"IIII", b"1X3=", 1, 4),
            (b"GCGT", b"IIII", b"1X3=", 1, 4),
            (b"TCGT", b"IIII", b"1X3=", 1, 4),
            (b"AAGT", b"IIII", b"1=1X2=", 1, 4),
            (b"AGGT", b"IIII", b"1=1X2=", 1, 4),
            (b"ATGT", b"IIII", b"1=1X2=", 1, 4),
            (b"ACAT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACCT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACTT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACGA", b"IIII", b"3=1X", 1, 4),
            (b"ACGC", b"IIII", b"3=1X", 1, 4),
            (b"ACGG", b"IIII", b"3=1X", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_invalid(self):
        """ Test when the read has an invalid substitution. """
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
        relvecs = [
            # 1 deletion
            [MATCH, DELET, MATCH, MATCH],
            [MATCH, MATCH, DELET, MATCH],
            # 2 deletions
            [MATCH, DELET, DELET, MATCH],
        ]
        expects = [
            # 1 deletion
            (b"AGT", b"III", b"1=1D2=", 1, 4),
            (b"ACT", b"III", b"2=1D1=", 1, 4),
            # 2 deletions
            (b"AT", b"II", b"1=2D1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_delete_invalid(self):
        """ Test when the read has a deletion at either end. """
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
            (b"ANCGT", b"IIIII", b"1=1I3=", 1, 4),
            (b"ACNGT", b"IIIII", b"2=1I2=", 1, 4),
            (b"ACGNT", b"IIIII", b"3=1I1=", 1, 4),
            (b"ANCG", b"IIII", b"1=1I2=", 1, 3),
            (b"ACNG", b"IIII", b"2=1I1=", 1, 3),
            (b"CNGT", b"IIII", b"1=1I2=", 2, 4),
            (b"CGNT", b"IIII", b"2=1I1=", 2, 4),
            (b"CNG", b"III", b"1=1I1=", 2, 3),
            # 2 insertions, 1 base apart
            (b"ANCNGT", b"IIIIII", b"1=1I1=1I2=", 1, 4),
            (b"ACNGNT", b"IIIIII", b"2=1I1=1I1=", 1, 4),
            # 2 insertions, 2 bases apart
            (b"ANCGNT", b"IIIIII", b"1=1I2=1I1=", 1, 4),
            # 3 insertions, 1 base apart
            (b"ANCNGNT", b"IIIIIII", b"1=1I1=1I1=1I1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert_end5_invalid(self):
        """ Test when the read has an insertion at the 5' end. """
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
        ref = DNA(b"ACGT")
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
            (b"CNCGT", b"IIIII", b"1X1I3=", 1, 4),
            (b"CCNGT", b"IIIII", b"1X1=1I2=", 1, 4),
            (b"ANTGT", b"IIIII", b"1=1I1X2=", 1, 4),
            (b"ATNGT", b"IIIII", b"1=1X1I2=", 1, 4),
            (b"ATGNT", b"IIIII", b"1=1X1=1I1=", 1, 4),
            # 1 insertion next to 1 low-quality base call
            (b"NNCGT", b"!IIII", b"1M1I3=", 1, 4),
            (b"NCNGT", b"!IIII", b"1M1=1I2=", 1, 4),
            (b"ANNGT", b"II!II", b"1=1I1M2=", 1, 4),
            (b"ANNGT", b"I!III", b"1=1M1I2=", 1, 4),
            (b"ANGNT", b"I!III", b"1=1M1=1I1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)


class TestAsSam(ut.TestCase):
    """ Test function `as_sam`. """

    def test_line_in_sam_format(self):
        line = as_sam(b"FS10000136:77:BPG61616-2310:1:1101:1000:1300", 99,
                      "SARS2_FSE", 1, 42, b"151=", "=", 133, 283,
                      DNA(b"CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTG"
                          b"GAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTC"
                          b"AGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC"),
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:")
        expect = (b"FS10000136:77:BPG61616-2310:1:1101:1000:1300\t99\tSARS2_FSE"
                  b"\t1\t42\t151=\t=\t133\t283\t"
                  b"CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTT"
                  b"ATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCG"
                  b"TTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC\t"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:\n")
        self.assertEqual(line, expect)
