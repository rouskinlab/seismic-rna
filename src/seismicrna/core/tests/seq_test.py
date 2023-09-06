"""

Tests for Sequence Core Module.

========================================================================

"""

from itertools import combinations, product
from string import printable
import unittest as ut

import numpy as np
from ..seq import Seq, DNA, RNA, expand_degenerate_seq


class TestDna(ut.TestCase):
    """ Test class `DNA`. """

    def test_alph(self):
        self.assertEqual(DNA.alph, ('A', 'C', 'N', 'G', 'T'))

    def test_get_comp(self):
        self.assertEqual(DNA.get_comp(), ('T', 'G', 'N', 'C', 'A'))

    def test_get_comptrans(self):
        self.assertEqual(DNA.get_comptrans(),
                         {65: 'T', 67: 'G', 71: 'C', 78: 'N', 84: 'A'})

    def test_get_alphaset(self):
        self.assertEqual(DNA.get_alphaset(), {'A', 'C', 'G', 'T', 'N'})

    def test_valid(self):
        """ Test whether valid DNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGTN"] * length)):
                dna = DNA("".join(bases))
                self.assertEqual(len(dna), length)
                self.assertEqual(str(dna), "".join(bases))

    def test_random(self):
        """ Test whether random DNA sequences can be created. """
        for length in range(1, 5):
            dna = DNA.random(length)
            self.assertTrue(isinstance(dna, DNA))
            self.assertEqual(len(dna), length)

    def test_to_array(self):
        """ Test generating NumPy arrays from DNA sequences. """
        array = DNA("CTANG").to_array()
        self.assertEqual(array.dtype, np.dtype("<U1"))
        self.assertTrue(np.all(array == np.array(['C', 'T', 'A', 'N', 'G'])))

    def test_slice(self):
        """ Test slicing DNA sequences. """
        for length in range(2, 10):
            dna = DNA.random(length)
            seq = str(dna)
            for i, j in combinations(range(length + 1), r=2):
                subseq = dna[i: j]
                self.assertTrue(isinstance(subseq, DNA))
                self.assertEqual(str(subseq), seq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing DNA sequences. """
        seqs = ["ACGTN", "GTCAGCTGCANTGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        recs = ["NACGT", "CATGCANTGCAGCTGAC", "AGTATGATGATGTCCCCCCACTTTA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = DNA(seq)
                rev = DNA(rec)
                self.assertTrue(isinstance(fwd.rc, DNA))
                self.assertTrue(isinstance(rev.rc, DNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_transcribe(self):
        """ Test transcribing DNA sequences. """
        dseqs = ["ACGTN", "GTCAGCTGCANTGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        rseqs = ["ACGUN", "GUCAGCUGCANUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        for dna, rna in zip(dseqs, rseqs, strict=True):
            with self.subTest(dna=dna, rna=rna):
                tr = DNA(dna).tr()
                self.assertTrue(isinstance(tr, RNA))
                self.assertEqual(tr, RNA(rna))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGTN":
                self.assertRaises(ValueError, DNA, char)

    def test_bool(self):
        """ Test that only zero-length DNA sequences are falsy. """
        self.assertFalse(DNA(""))
        for length in range(1, 10):
            self.assertTrue(DNA.random(length))


class TestRna(ut.TestCase):
    """ Test class `RNA`. """

    def test_alph(self):
        self.assertEqual(RNA.alph, ('A', 'C', 'N', 'G', 'U'))

    def test_get_comp(self):
        self.assertEqual(RNA.get_comp(), ('U', 'G', 'N', 'C', 'A'))

    def test_get_comptrans(self):
        self.assertEqual(RNA.get_comptrans(),
                         {65: 'U', 67: 'G', 71: 'C', 78: 'N', 85: 'A'})

    def test_get_alphaset(self):
        self.assertEqual(RNA.get_alphaset(), {'A', 'C', 'G', 'U', 'N'})

    def test_valid(self):
        """ Test whether valid RNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGUN"] * length)):
                rna = RNA("".join(bases))
                self.assertEqual(len(rna), length)
                self.assertEqual(str(rna), "".join(bases))

    def test_random(self):
        """ Test whether random RNA sequences can be created. """
        for length in range(1, 5):
            rna = RNA.random(length)
            self.assertTrue(isinstance(rna, RNA))
            self.assertEqual(len(rna), length)

    def test_to_array(self):
        """ Test generating NumPy arrays from RNA sequences. """
        array = RNA("NGAUC").to_array()
        self.assertEqual(array.dtype, np.dtype("<U1"))
        self.assertTrue(np.all(array == np.array(['N', 'G', 'A', 'U', 'C'])))

    def test_slice(self):
        """ Test slicing RNA sequences. """
        for length in range(2, 10):
            rna = RNA.random(length)
            seq = str(rna)
            for i, j in combinations(range(length + 1), r=2):
                subseq = rna[i: j]
                self.assertTrue(isinstance(subseq, RNA))
                self.assertEqual(str(subseq), seq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing RNA sequences. """
        seqs = ["ACGUN", "GUCAGCUGCANUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        recs = ["NACGU", "CAUGCANUGCAGCUGAC", "AGUAUGAUGAUGUCCCCCCACUUUA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = RNA(seq)
                rev = RNA(rec)
                self.assertTrue(isinstance(fwd.rc, RNA))
                self.assertTrue(isinstance(rev.rc, RNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_reverse_transcribe(self):
        """ Test reverse transcribing RNA sequences. """
        rseqs = ["ACGUN", "GUCAGCUGCANUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        dseqs = ["ACGTN", "GTCAGCTGCANTGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        for rna, dna in zip(rseqs, dseqs, strict=True):
            with self.subTest(rna=rna, dna=dna):
                rt = RNA(rna).rt()
                self.assertTrue(isinstance(rt, DNA))
                self.assertEqual(rt, DNA(dna))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGUN":
                self.assertRaises(ValueError, RNA, char)

    def test_bool(self):
        """ Test that only zero-length RNA sequences are falsy. """
        self.assertFalse(RNA(""))
        for length in range(1, 10):
            self.assertTrue(RNA.random(length))


class TestSeq(ut.TestCase):
    """ Test basic properties of the Seq class. """

    def test_abstract_base_class(self):
        """ Test that instantiating a Seq raises AttributeError. """
        self.assertRaises(AttributeError, Seq, "ACGN")

    def test_equal_dna_dna(self):
        """ Test that DNA instances with the same sequences compare as
        equal. """
        seq = "ACGTN"
        self.assertEqual(DNA(seq), DNA(seq))

    def test_equal_rna_rna(self):
        """ Test that RNA instances with the same sequences compare as
        equal. """
        seq = "ACGUN"
        self.assertEqual(RNA(seq), RNA(seq))

    def test_not_equal_dna_str(self):
        """ Test that DNA and str instances with the same sequences
        compare as not equal. """
        seq = "ACGTN"
        self.assertNotEqual(seq, DNA(seq))
        self.assertNotEqual(DNA(seq), seq)

    def test_not_equal_rna_str(self):
        """ Test that RNA and str instances with the same sequences
        compare as not equal. """
        seq = "ACGUN"
        self.assertNotEqual(seq, RNA(seq))
        self.assertNotEqual(RNA(seq), seq)

    def test_not_equal_dna_rna(self):
        """ Test that DNA and RNA instances with the same sequences
        compare as not equal. """
        seq = "ACGN"
        self.assertNotEqual(DNA(seq), RNA(seq))
        self.assertNotEqual(RNA(seq), DNA(seq))

    def test_hashable_dna(self):
        """ Test that DNA instances are hashable. """
        self.assertTrue(isinstance(hash(DNA("ACGTN")), int))

    def test_hashable_rna(self):
        """ Test that RNA instances are hashable. """
        self.assertTrue(isinstance(hash(RNA("ACGUN")), int))

    def test_set_str_dna_rna(self):
        """ Test that instances of str, DNA, and RNA with identical
        sequences can all be used together in a set. """
        seq = "ACGN"
        seqs = {seq, DNA(seq), RNA(seq)}
        self.assertEqual(len(seqs), 3)
        self.assertTrue(seq in seqs)
        self.assertTrue(DNA(seq) in seqs)
        self.assertTrue(RNA(seq) in seqs)

    def test_dict_str_dna_rna(self):
        """ Test that instances of str, DNA, and RNA with identical
        sequences can all be used together as dict keys. """
        seq = "ACGN"
        seqs = {seq: None, DNA(seq): None, RNA(seq): None}
        self.assertEqual(len(seqs), 3)
        self.assertTrue(seq in seqs)
        self.assertTrue(DNA(seq) in seqs)
        self.assertTrue(RNA(seq) in seqs)


class TestExpandDegenerateSeq(ut.TestCase):
    """ Test function `expand_degenerate_seq`. """

    def test_zero_degenerate(self):
        """ Test that the original sequence is returned. """
        self.assertEqual(list(expand_degenerate_seq(DNA("ACGT"))),
                         list(map(DNA, ["ACGT"])))

    def test_one_degenerate(self):
        """ Test that one sequence is returned for each DNA base. """
        self.assertEqual(list(expand_degenerate_seq(DNA("ACNT"))),
                         list(map(DNA, ["ACAT", "ACCT", "ACGT", "ACTT"])))

    def test_two_degenerate(self):
        """ Test that one sequence is returned for every combination of
        two DNA bases. """
        self.assertEqual(list(expand_degenerate_seq(DNA("NCGN"))),
                         list(map(DNA, ["ACGA", "ACGC", "ACGG", "ACGT",
                                        "CCGA", "CCGC", "CCGG", "CCGT",
                                        "GCGA", "GCGC", "GCGG", "GCGT",
                                        "TCGA", "TCGC", "TCGG", "TCGT"])))
