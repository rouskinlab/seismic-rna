import unittest as ut

import numpy as np

from ..code import (IRREC,
                    INDEL,
                    NOCOV,
                    MATCH,
                    DELET,
                    INS_5,
                    INS_3,
                    INSRT,
                    SUB_A,
                    SUB_C,
                    SUB_G,
                    SUB_T,
                    SUB_N,
                    ANY_B,
                    ANY_D,
                    ANY_H,
                    ANY_V,
                    ANY_N,
                    REL_SIZE,
                    REL_TYPE)


class TestConstants(ut.TestCase):

    def test_rel_size_type(self):
        self.assertIs(REL_SIZE, 1)
        self.assertIs(REL_TYPE, np.uint8)

    def test_primary_codes(self):
        """ Test the primary relation codes. """
        for exp, code in enumerate([MATCH, DELET, INS_5, INS_3,
                                    SUB_A, SUB_C, SUB_G, SUB_T]):
            self.assertIsInstance(code, int)
            self.assertEqual(code, 2 ** exp)

    def test_derived_codes(self):
        """ Test the derived relation codes. """
        self.assertEqual(IRREC, 0)
        self.assertEqual(INSRT, 12)
        self.assertEqual(INDEL, 14)
        self.assertEqual(SUB_N, 240)
        self.assertEqual(ANY_N, 241)
        self.assertEqual(NOCOV, 255)
        self.assertEqual(ANY_B, ANY_N - SUB_A)
        self.assertEqual(ANY_D, ANY_N - SUB_C)
        self.assertEqual(ANY_H, ANY_N - SUB_G)
        self.assertEqual(ANY_V, ANY_N - SUB_T)
