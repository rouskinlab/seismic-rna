import unittest as ut

from ..encode import encode_match, encode_relate
from ...core.qual import LO_QUAL, HI_QUAL
from ...core.rel import (MATCH,
                         SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                         ANY_B, ANY_D, ANY_H, ANY_V, ANY_N)


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
