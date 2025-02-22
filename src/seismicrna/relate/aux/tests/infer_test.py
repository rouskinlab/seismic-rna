import unittest as ut
from typing import Iterable

from seismicrna.relate.aux.infer import infer_read
from seismicrna.core.rel import (DELET,
                                 INS_5,
                                 INS_3,
                                 ANY_N,
                                 SUB_A,
                                 SUB_C,
                                 SUB_G,
                                 SUB_T)
from seismicrna.core.seq import DNA


class TestInferRead(ut.TestCase):

    def assert_equal(self,
                     refseq: DNA,
                     relvecs: Iterable[tuple[int, int, dict[int, int]]],
                     expects: Iterable[tuple[str, str, str]]):
        """ Assert that the actual and expected outputs match. """
        for relvec, expect in zip(relvecs, expects, strict=True):
            with self.subTest(relvec=relvec, expect=expect):
                end5, end3, muts = relvec
                read, qual, cigar = expect
                self.assertEqual(infer_read(refseq, end5, end3, muts),
                                 (DNA(read), qual, cigar))

    def assert_raise(self,
                     ref: DNA,
                     relvecs: Iterable[tuple[int, int, dict[int, int]]],
                     error: type[Exception],
                     regex: str):
        for relvec in relvecs:
            with self.subTest(relvec=relvec):
                self.assertRaisesRegex(error, regex, infer_read, ref, *relvec)

    def test_all_match(self):
        """ Test when the read has four matching bases. """
        ref = DNA("ACGT")
        relvecs = [(1, 4, {})]
        expects = [("ACGT", "IIII", "4=")]
        self.assert_equal(ref, relvecs, expects)

    def test_all_match_n(self):
        """ Test when the read has four matching bases and an ambiguous
        base. """
        ref = DNA("ACNGT")
        relvecs = [(1, 5, {3: ANY_N})]
        expects = [("ACNGT", "II!II", "2=1M2=")]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_valid(self):
        """ Test when the read does not cover one or both ends of the
        reference. """
        ref = DNA("ACGT")
        relvecs = [
            (2, 4, {}),
            (1, 3, {}),
            (2, 3, {}),
            (3, 4, {}),
            (1, 2, {}),
        ]
        expects = [
            ("CGT", "III", "3="),
            ("ACG", "III", "3="),
            ("CG", "II", "2="),
            ("GT", "II", "2="),
            ("AC", "II", "2="),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_low_qual_valid(self):
        """ Test when the read has a low-quality base. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {1: ANY_N - SUB_A}),
            (1, 4, {2: ANY_N - SUB_C}),
            (1, 4, {3: ANY_N - SUB_G}),
            (1, 4, {4: ANY_N - SUB_T}),
        ]
        expects = [
            ("NCGT", "!III", "1M3="),
            ("ANGT", "I!II", "1=1M2="),
            ("ACNT", "II!I", "2=1M1="),
            ("ACGN", "III!", "3=1M"),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_valid(self):
        """ Test when the read has a substitution. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {1: SUB_C}),
            (1, 4, {1: SUB_G}),
            (1, 4, {1: SUB_T}),
            (1, 4, {2: SUB_A}),
            (1, 4, {2: SUB_G}),
            (1, 4, {2: SUB_T}),
            (1, 4, {3: SUB_A}),
            (1, 4, {3: SUB_C}),
            (1, 4, {3: SUB_T}),
            (1, 4, {4: SUB_A}),
            (1, 4, {4: SUB_C}),
            (1, 4, {4: SUB_G}),
        ]
        expects = [
            ("CCGT", "IIII", "1X3="),
            ("GCGT", "IIII", "1X3="),
            ("TCGT", "IIII", "1X3="),
            ("AAGT", "IIII", "1=1X2="),
            ("AGGT", "IIII", "1=1X2="),
            ("ATGT", "IIII", "1=1X2="),
            ("ACAT", "IIII", "2=1X1="),
            ("ACCT", "IIII", "2=1X1="),
            ("ACTT", "IIII", "2=1X1="),
            ("ACGA", "IIII", "3=1X"),
            ("ACGC", "IIII", "3=1X"),
            ("ACGG", "IIII", "3=1X"),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_invalid(self):
        """ Test when the read has an invalid substitution. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {1: SUB_A}),
            (1, 4, {2: SUB_C}),
            (1, 4, {3: SUB_G}),
            (1, 4, {4: SUB_T}),
        ]
        self.assert_raise(ref,
                          relvecs,
                          ValueError,
                          "Cannot substitute [ACGT] to itself")

    def test_delete_valid(self):
        """ Test when the read has deletions. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 deletion
            (1, 4, {2: DELET}),
            (1, 4, {3: DELET}),
            # 2 deletions
            (1, 4, {2: DELET, 3: DELET}),
        ]
        expects = [
            # 1 deletion
            ("AGT", "III", "1=1D2="),
            ("ACT", "III", "2=1D1="),
            # 2 deletions
            ("AT", "II", "1=2D1="),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_delete_invalid(self):
        """ Test when the read has a deletion at either end. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {1: DELET}),
            (2, 4, {2: DELET}),
            (3, 4, {3: DELET}),
            (4, 4, {4: DELET}),
            (1, 4, {1: DELET, 2: DELET, 3: DELET, 4: DELET}),
            (1, 4, {1: DELET, 4: DELET}),
            (1, 4, {4: DELET}),
            (1, 3, {3: DELET}),
            (1, 2, {2: DELET}),
            (1, 1, {1: DELET}),
        ]
        self.assert_raise(ref,
                          relvecs,
                          ValueError,
                          "Deletion cannot be at position [0-9]+ in .+")

    def test_insert5_valid(self):
        """ Test when the read has insertions. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion
            (1, 4, {1: INS_5}),
            (1, 4, {2: INS_5}),
            (1, 4, {3: INS_5}),
            (1, 3, {1: INS_5}),
            (1, 3, {2: INS_5}),
            (2, 4, {2: INS_5}),
            (2, 4, {3: INS_5}),
            (2, 3, {2: INS_5}),
            # 2 insertions, 1 base apart
            (1, 4, {1: INS_5, 2: INS_5}),
            (1, 4, {2: INS_5, 3: INS_5}),
            # 2 insertions, 2 bases apart
            (1, 4, {1: INS_5, 3: INS_5}),
            # 3 insertions, 1 base apart
            (1, 4, {1: INS_5, 2: INS_5, 3: INS_5}),
        ]
        expects = [
            # 1 insertion
            ("ANCGT", "IIIII", "1=1I3="),
            ("ACNGT", "IIIII", "2=1I2="),
            ("ACGNT", "IIIII", "3=1I1="),
            ("ANCG", "IIII", "1=1I2="),
            ("ACNG", "IIII", "2=1I1="),
            ("CNGT", "IIII", "1=1I2="),
            ("CGNT", "IIII", "2=1I1="),
            ("CNG", "III", "1=1I1="),
            # 2 insertions, 1 base apart
            ("ANCNGT", "IIIIII", "1=1I1=1I2="),
            ("ACNGNT", "IIIIII", "2=1I1=1I1="),
            # 2 insertions, 2 bases apart
            ("ANCGNT", "IIIIII", "1=1I2=1I1="),
            # 3 insertions, 1 base apart
            ("ANCNGNT", "IIIIIII", "1=1I1=1I1=1I1="),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert3_valid(self):
        """ Test when the read has insertions. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion
            (1, 4, {2: INS_3}),
            (1, 4, {3: INS_3}),
            (1, 4, {4: INS_3}),
            (1, 3, {2: INS_3}),
            (1, 3, {3: INS_3}),
            (2, 4, {3: INS_3}),
            (2, 4, {4: INS_3}),
            (2, 3, {3: INS_3}),
            # 2 insertions, 1 base apart
            (1, 4, {2: INS_3, 3: INS_3}),
            (1, 4, {3: INS_3, 4: INS_3}),
            # 2 insertions, 2 bases apart
            (1, 4, {2: INS_3, 4: INS_3}),
            # 3 insertions, 1 base apart
            (1, 4, {2: INS_3, 3: INS_3, 4: INS_3}),
        ]
        expects = [
            # 1 insertion
            ("ANCGT", "IIIII", "1=1I3="),
            ("ACNGT", "IIIII", "2=1I2="),
            ("ACGNT", "IIIII", "3=1I1="),
            ("ANCG", "IIII", "1=1I2="),
            ("ACNG", "IIII", "2=1I1="),
            ("CNGT", "IIII", "1=1I2="),
            ("CGNT", "IIII", "2=1I1="),
            ("CNG", "III", "1=1I1="),
            # 2 insertions, 1 base apart
            ("ANCNGT", "IIIIII", "1=1I1=1I2="),
            ("ACNGNT", "IIIIII", "2=1I1=1I1="),
            # 2 insertions, 2 bases apart
            ("ANCGNT", "IIIIII", "1=1I2=1I1="),
            # 3 insertions, 1 base apart
            ("ANCNGNT", "IIIIIII", "1=1I1=1I1=1I1="),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert5_end5_invalid(self):
        """ Test when the read has an insertion at the 5' end. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {4: INS_5}),
            (1, 4, {3: INS_5, 4: INS_5}),
            (1, 4, {2: INS_5, 3: INS_5, 4: INS_5}),
            (1, 4, {1: INS_5, 2: INS_5, 3: INS_5, 4: INS_5}),
            (2, 3, {3: INS_5}),
            (2, 3, {2: INS_5, 3: INS_5}),
        ]
        self.assert_raise(ref,
                          relvecs,
                          ValueError,
                          "Position [0-9]+ in .+ cannot be 5' of an insertion")

    def test_insert3_end3_invalid(self):
        """ Test when the read has an insertion at the 3' end. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {1: INS_3}),
            (1, 4, {1: INS_3, 2: INS_3}),
            (1, 4, {1: INS_3, 2: INS_3, 3: INS_3}),
            (1, 4, {1: INS_3, 2: INS_3, 3: INS_3, 4: INS_3}),
            (2, 3, {2: INS_3}),
            (2, 3, {2: INS_3, 3: INS_3}),
        ]
        self.assert_raise(ref,
                          relvecs,
                          ValueError,
                          "Position [0-9]+ in .+ cannot be 3' of an insertion")

    def test_insert_deletion_invalid(self):
        """ Test when the read has an insertion next to a deletion. """
        ref = DNA("ACGT")
        relvecs = [
            (1, 4, {2: INS_3 + DELET}),
            (1, 4, {3: INS_3 + DELET}),
            (1, 4, {4: INS_3 + DELET}),
            (1, 4, {1: DELET + INS_5}),
            (1, 4, {2: DELET + INS_5}),
            (1, 4, {3: DELET + INS_5}),
        ]
        self.assert_raise(ref,
                          relvecs,
                          ValueError,
                          "Position .+ is del and ins")

    def test_insert5_non_match_valid(self):
        """ Test when the read has insertions next to substitutions or
        low-quality base calls. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion next to 1 substitution
            (1, 4, {1: SUB_C + INS_5}),
            (1, 4, {1: SUB_C, 2: INS_5}),
            (1, 4, {2: INS_3 + SUB_T}),
            (1, 4, {2: SUB_T + INS_5}),
            (1, 4, {2: SUB_T, 3: INS_5}),
            # 1 insertion next to 1 low-quality base call
            (1, 4, {1: ANY_N - SUB_A + INS_5}),
            (1, 4, {1: ANY_N - SUB_A, 2: INS_5}),
            (1, 4, {1: INS_5, 2: INS_3 + ANY_N - SUB_C}),
            (1, 4, {2: ANY_N - SUB_C + INS_5}),
            (1, 4, {2: ANY_N - SUB_C, 3: INS_5}),
        ]
        expects = [
            # 1 insertion next to 1 substitution
            ("CNCGT", "IIIII", "1X1I3="),
            ("CCNGT", "IIIII", "1X1=1I2="),
            ("ANTGT", "IIIII", "1=1I1X2="),
            ("ATNGT", "IIIII", "1=1X1I2="),
            ("ATGNT", "IIIII", "1=1X1=1I1="),
            # 1 insertion next to 1 low-quality base call
            ("NNCGT", "!IIII", "1M1I3="),
            ("NCNGT", "!IIII", "1M1=1I2="),
            ("ANNGT", "II!II", "1=1I1M2="),
            ("ANNGT", "I!III", "1=1M1I2="),
            ("ANGNT", "I!III", "1=1M1=1I1="),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert3_non_match_valid(self):
        """ Test when the read has insertions next to substitutions or
        low-quality base calls. """
        ref = DNA("ACGT")
        relvecs = [
            # 1 insertion next to 1 substitution
            (1, 4, {1: SUB_C, 2: INS_3}),
            (1, 4, {1: SUB_C, 3: INS_3}),
            (1, 4, {2: INS_3 + SUB_T}),
            (1, 4, {2: SUB_T, 3: INS_3}),
            (1, 4, {2: SUB_T, 4: INS_3}),
            # 1 insertion next to 1 low-quality base call
            (1, 4, {1: ANY_N - SUB_A, 2: INS_3}),
            (1, 4, {1: ANY_N - SUB_A, 3: INS_3}),
            (1, 4, {2: INS_3 + ANY_N - SUB_C}),
            (1, 4, {2: ANY_N - SUB_C, 3: INS_3}),
            (1, 4, {2: ANY_N - SUB_C, 4: INS_3}),
        ]
        expects = [
            # 1 insertion next to 1 substitution
            ("CNCGT", "IIIII", "1X1I3="),
            ("CCNGT", "IIIII", "1X1=1I2="),
            ("ANTGT", "IIIII", "1=1I1X2="),
            ("ATNGT", "IIIII", "1=1X1I2="),
            ("ATGNT", "IIIII", "1=1X1=1I1="),
            # 1 insertion next to 1 low-quality base call
            ("NNCGT", "!IIII", "1M1I3="),
            ("NCNGT", "!IIII", "1M1=1I2="),
            ("ANNGT", "II!II", "1=1I1M2="),
            ("ANNGT", "I!III", "1=1M1I2="),
            ("ANGNT", "I!III", "1=1M1=1I1="),
        ]
        self.assert_equal(ref, relvecs, expects)


if __name__ == "__main__":
    ut.main(verbosity=2)
