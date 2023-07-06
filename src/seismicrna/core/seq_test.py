"""
Core -- Sequence Testing Module
========================================================================
Auth: Matty

Unit tests for `core.seq`.
"""

from itertools import combinations, product
from string import printable
import unittest as ut

import numpy as np
from .seq import BASES, BASES_ARR, RBASE, DNA, RNA, expand_degenerate_seq


class TestConstants(ut.TestCase):
    """ Test constants of `seq` module. """

    def test_bases(self):
        """ Test that `BASES` contains the four DNA letters. """
        self.assertEqual(BASES, b"ACGT")

    def test_rbase(self):
        """ Test that `RBASE` contains the four RNA letters. """
        self.assertEqual(RBASE, b"ACGU")

    def test_bases_arr(self):
        """ Test that `BASES` contains the four DNA ASCII codes. """
        self.assertTrue(isinstance(BASES_ARR, np.ndarray))
        self.assertIs(BASES_ARR.dtype.type, np.uint8)
        self.assertEqual(BASES_ARR.tolist(), [65, 67, 71, 84])


class TestDna(ut.TestCase):
    """ Test class `DNA`. """

    def test_valid(self):
        """ Test whether valid DNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGT"] * length)):
                dna = DNA("".join(bases).encode())
                self.assertEqual(len(dna), length)
                self.assertEqual(dna.decode(), "".join(bases))

    def test_slice(self):
        """ Test slicing DNA sequences. """
        dnaseq = "GATTACA"
        dna = DNA(dnaseq.encode())
        for i, j in combinations(range(len(dna) + 1), r=2):
            subseq = dna[i: j]
            self.assertTrue(isinstance(subseq, DNA))
            self.assertEqual(subseq.decode(), dnaseq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing DNA sequences. """
        seqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        recs = ["ACGT", "CATGCATGCAGCTGAC", "AGTATGATGATGTCCCCCCACTTTA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = DNA(seq.encode())
                rev = DNA(rec.encode())
                self.assertTrue(isinstance(fwd.rc, DNA))
                self.assertTrue(isinstance(rev.rc, DNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_transcribe(self):
        """ Test transcribing DNA sequences. """
        dseqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        rseqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        for dna, rna in zip(dseqs, rseqs, strict=True):
            with self.subTest(dna=dna, rna=rna):
                tr = DNA(dna.encode()).tr()
                self.assertTrue(isinstance(tr, RNA))
                self.assertEqual(tr, RNA(rna.encode()))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGT":
                self.assertRaises(ValueError, DNA, char.encode())

    def test_zero(self):
        """ Test that zero-length DNA sequences raise ValueError. """
        self.assertRaises(ValueError, DNA, b"")


class TestRna(ut.TestCase):
    """ Test class `RNA`. """

    def test_valid(self):
        """ Test whether valid RNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGU"] * length)):
                rna = RNA("".join(bases).encode())
                self.assertEqual(len(rna), length)
                self.assertEqual(rna.decode(), "".join(bases))

    def test_slice(self):
        """ Test slicing RNA sequences. """
        rnaseq = "GAUUACA"
        rna = RNA(rnaseq.encode())
        for i, j in combinations(range(len(rna) + 1), r=2):
            subseq = rna[i: j]
            self.assertTrue(isinstance(subseq, RNA))
            self.assertEqual(subseq.decode(), rnaseq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing RNA sequences. """
        seqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        recs = ["ACGU", "CAUGCAUGCAGCUGAC", "AGUAUGAUGAUGUCCCCCCACUUUA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = RNA(seq.encode())
                rev = RNA(rec.encode())
                self.assertTrue(isinstance(fwd.rc, RNA))
                self.assertTrue(isinstance(rev.rc, RNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_reverse_transcribe(self):
        """ Test reverse transcribing RNA sequences. """
        rseqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        dseqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        for rna, dna in zip(rseqs, dseqs, strict=True):
            with self.subTest(rna=rna, dna=dna):
                rt = RNA(rna.encode()).rt()
                self.assertTrue(isinstance(rt, DNA))
                self.assertEqual(rt, DNA(dna.encode()))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGU":
                self.assertRaises(ValueError, RNA, char.encode())

    def test_zero(self):
        """ Test that zero-length RNA sequences raise ValueError. """
        self.assertRaises(ValueError, RNA, b"")


class TestExpandDegenerateSeq(ut.TestCase):
    """ Test function `expand_degenerate_seq`. """

    def test_zero_degenerate(self):
        """ Test that the original sequence is returned. """
        self.assertEqual(list(expand_degenerate_seq(b"ACGT")),
                         [b"ACGT"])

    def test_one_degenerate(self):
        """ Test that one sequence is returned for each DNA base. """
        self.assertEqual(list(expand_degenerate_seq(b"ACNT")),
                         [b"ACAT", b"ACCT", b"ACGT", b"ACTT"])

    def test_two_degenerate(self):
        """ Test that one sequence is returned for every combination of
        two DNA bases. """
        self.assertEqual(list(expand_degenerate_seq(b"NCGN")),
                         [b"ACGA", b"ACGC", b"ACGG", b"ACGT",
                          b"CCGA", b"CCGC", b"CCGG", b"CCGT",
                          b"GCGA", b"GCGC", b"GCGG", b"GCGT",
                          b"TCGA", b"TCGC", b"TCGG", b"TCGT"])
