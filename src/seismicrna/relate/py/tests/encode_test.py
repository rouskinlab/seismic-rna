import unittest as ut
from string import printable

from ..encode import SUBS_ENCODINGS, SUBS_DECODINGS, is_acgt, encode_relate
from ....core.ngs import LO_QUAL, HI_QUAL
from ....core.rel import MATCH, ANY_N
from ....core.seq import DNA


class TestSubsCoding(ut.TestCase):

    def test_subs_encodings(self):
        self.assertEqual(SUBS_ENCODINGS,
                         {"A": 16, "C": 32, "G": 64, "T": 128})

    def test_subs_decodings(self):
        self.assertEqual(SUBS_DECODINGS,
                         {16: "A", 32: "C", 64: "G", 128: "T"})


class TestIsACGT(ut.TestCase):

    def test_is_acgt(self):
        for char in "ACGT":
            self.assertTrue(is_acgt(char))
        for char in "acgtnN":
            self.assertFalse(is_acgt(char))
        for char in printable:
            self.assertEqual(is_acgt(char), char in DNA.four())


class TestEncodeRelate(ut.TestCase):

    def test_encode_relate_hi_qual(self):
        """ Test when the quality is at least the minimum. """
        for ref in "ACGTN":
            for read in "ACGTN":
                code = encode_relate(ref, read, HI_QUAL, HI_QUAL)
                if ref == "N":
                    self.assertEqual(code, ANY_N)
                elif read == "N":
                    self.assertEqual(code, ANY_N ^ SUBS_ENCODINGS[ref])
                elif read == ref:
                    self.assertEqual(code, MATCH)
                else:
                    self.assertEqual(code, SUBS_ENCODINGS[read])

    def test_encode_relate_lo_qual(self):
        """ Test when the quality is less than the minimum. """
        for ref in "ACGTN":
            for read in "ACGTN":
                code = encode_relate(ref, read, LO_QUAL, HI_QUAL)
                if ref == "N":
                    self.assertEqual(code, ANY_N)
                else:
                    self.assertEqual(code, ANY_N ^ SUBS_ENCODINGS[ref])
