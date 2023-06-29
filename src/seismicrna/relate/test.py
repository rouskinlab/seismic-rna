import unittest as ut
from sys import byteorder

import pandas as pd

from .relate import relate_line
from .seqpos import format_seq_pos
from ..core.rel import as_sam, iter_alignments, NOCOV, MED_QUAL
from ..core.seq import DNA


class TestSeqposFormatSeqPos(ut.TestCase):
    """ Test function `seqpos.format_seq_pos`. """

    def test_acgt_index_1_acgt(self):
        """ Test with ACGT, 1-indexed. """
        index = format_seq_pos(DNA(b"ACGT"), [1, 2, 3, 4], 1)
        expect = pd.Index(["A1", "C2", "G3", "T4"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_1_cg(self):
        """ Test with ACGT, 1-indexed. """
        index = format_seq_pos(DNA(b"ACGT"), [2, 3], 1)
        expect = pd.Index(["C2", "G3"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_58_acgt(self):
        """ Test with ACGT, 58-indexed. """
        index = format_seq_pos(DNA(b"ACGT"), [58, 59, 60, 61], 58)
        expect = pd.Index(["A58", "C59", "G60", "T61"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_58_cg(self):
        """ Test with ACGT, 58-indexed. """
        index = format_seq_pos(DNA(b"ACGT"), [59, 60], 58)
        expect = pd.Index(["C59", "G60"])
        self.assertTrue(index.equals(expect))


class TestRelateRelateLineAmbrel(ut.TestCase):
    """ Test function `relate.relate_line`. """

    @staticmethod
    def relate(ref: str, refseq: DNA, read, qual, cigar, end5):
        """ Generate a SAM line from the given information, and use it
        to compute a relation vector. """
        sam_line = as_sam(b"read", 99, ref, end5, 0, cigar, "=", 1, len(read),
                          read, qual)
        muts = bytearray(NOCOV.to_bytes(1, byteorder) * len(refseq))
        relate_line(sam_line, muts, refseq, len(refseq), ref, MED_QUAL, True)
        return muts

    def iter_cases(self, refseq: DNA, max_ins: int = 2):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, relvec in iter_alignments(refseq,
                                                                     max_ins,
                                                                     max_ins,
                                                                     max_ins):
            result = self.relate("ref", refseq, read, qual, cigar, end5)
            with self.subTest(relvec=relvec, result=result):
                self.assertEqual(relvec, result)

    def test_aaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAA. """
        self.iter_cases(DNA(b"AAAA"), 0)

    @ut.skip("Takes a long time to run")
    def test_aaaaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAAAA. """
        self.iter_cases(DNA(b"AAAAAA"), 0)

    def test_aacc_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from AACC. """
        self.iter_cases(DNA(b"AACC"), 1)

    def test_acgt_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from ACGT. """
        self.iter_cases(DNA(b"ACGT"), 1)
