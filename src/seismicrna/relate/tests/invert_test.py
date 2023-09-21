import unittest as ut

import numpy as np

from ..invert import find_relvec_ends, inverse_relate
from ...core.qual import HI_QUAL, LO_QUAL
from ...core.rel import (MATCH, DELET, INS_5, INS_3, INS_8,
                         NOCOV, ANY_N, MINS5, MINS3, ANY_8,
                         SUB_A, SUB_C, SUB_G, SUB_T, NP_TYPE)
from ...core.seq import DNA


class TestFindRelVecEnds(ut.TestCase):
    """ Test function `find_relvec_ends`. """

    def test_matches(self):
        """ Test that a relation vector of all matches is valid. """
        self.assertIsNotNone(find_relvec_ends(np.ones(8, NP_TYPE)))

    def test_nocov_margin(self):
        """ Test that a relation vector with non-covered positions on
        its margins is valid. """
        relvec = np.full(8, NOCOV, NP_TYPE)
        relvec[4] = 1
        self.assertIsNotNone(find_relvec_ends(relvec))

    def test_nocov_middle(self):
        """ Test that a relation vector with non-covered positions in
        the middle is invalid. """
        relvec = np.full(8, MATCH, NP_TYPE)
        relvec[4] = NOCOV
        self.assertRaises(ValueError, find_relvec_ends, relvec)

    def test_blank(self):
        """ Test that an entirely blank relation vector is invalid. """
        self.assertRaises(ValueError, find_relvec_ends,
                          np.full(8, NOCOV, NP_TYPE))

    def test_not_array(self):
        """ Test that a non-array relation vector is invalid. """
        self.assertRaises(TypeError, find_relvec_ends,
                          np.ones(8, NP_TYPE).tolist())

    def test_not_uint8(self):
        """ Test that a non-uint8 relation vector is invalid. """
        self.assertRaises(TypeError, find_relvec_ends, np.ones(8, np.int8))


class TestInverseRelate(ut.TestCase):
    """ Test function `inverse_relate`. """

    def assert_equal(self, ref: DNA,
                     relvecs: list[list[int]],
                     expects: list[tuple[str, str, str, int, int]]):
        """ Assert that the actual and expected outputs match. """
        for relvec, expect in zip(relvecs, expects, strict=True):
            with self.subTest(relvec=relvec, expect=expect):
                self.assertEqual(inverse_relate(ref, np.array(relvec,
                                                              dtype=np.uint8),
                                                HI_QUAL, LO_QUAL),
                                 (DNA(expect[0]),) + expect[1:])

    def assert_raise(self, ref: DNA,
                     relvecs: list[list[int]],
                     error: type[Exception],
                     regex: str):
        """ Assert that the relation vectors raise an exception. """
        for relvec in relvecs:
            with self.subTest(relvec=relvec):
                self.assertRaisesRegex(error, regex, inverse_relate,
                                       ref, np.array(relvec,
                                                     dtype=np.uint8),
                                       HI_QUAL, LO_QUAL)

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

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
