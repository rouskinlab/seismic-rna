import unittest as ut
from itertools import chain
from typing import Generator, Sequence

import numpy as np

from ..iterrelv import iter_relvecs_q53, iter_relvecs_all
from ...core.rel import (MATCH, DELET, INS_5, INS_3,
                         SUB_A, SUB_C, SUB_G, SUB_T,
                         ANY_B, ANY_D, ANY_H, ANY_V,
                         ANY_N, NOCOV)
from ...core.seq import DNA, expand_degenerate_seq


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

########################################################################
#                                                                      #
# ©2023, the Rouskin Lab.                                              #
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
