import unittest as ut

from ..encode import encode_match, encode_relate
from ....core.ngs import LO_QUAL, HI_QUAL
from ....core.rel import (MATCH,
                          SUB_A,
                          SUB_C,
                          SUB_G,
                          SUB_T,
                          SUB_N,
                          ANY_B,
                          ANY_D,
                          ANY_H,
                          ANY_V,
                          ANY_N)


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
                code = encode_relate(ref, read, HI_QUAL, HI_QUAL)
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
                self.assertEqual(encode_relate(ref, read, LO_QUAL, HI_QUAL),
                                 sub)


class TestEncodeMatch(ut.TestCase):
    """ Test function `encode_match`. """

    def test_encode_match_hi_qual(self):
        """ Test when the quality is at least the minimum. """
        for read in "ACGT":
            self.assertEqual(encode_match(read, HI_QUAL, HI_QUAL), MATCH)
        self.assertEqual(encode_match('N', HI_QUAL, HI_QUAL), ANY_N)

    def test_encode_match_lo_qual(self):
        """ Test when the quality is less than the minimum. """
        for read, sub in zip("ACGTN", [ANY_B, ANY_D, ANY_H, ANY_V, ANY_N]):
            self.assertEqual(encode_match(read, LO_QUAL, HI_QUAL), sub)

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
