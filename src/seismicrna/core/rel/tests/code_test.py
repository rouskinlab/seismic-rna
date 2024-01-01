"""

Tests for the Relate Core Module

========================================================================

"""

import unittest as ut

import numpy as np

from ..code import (IRREC,
                    INDEL,
                    NOCOV,
                    MATCH,
                    DELET,
                    INS_5,
                    INS_3,
                    INS_8,
                    MINS5,
                    MINS3,
                    ANY_8,
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
                    REL_TYPE)


class TestConstants(ut.TestCase):
    """ Test constants of `rel` module. """

    def test_np_type(self):
        self.assertIs(REL_TYPE, np.uint8)

    def test_primary_codes(self):
        """ Test the primary relation codes. """
        for exp, code in enumerate([MATCH, DELET, INS_5, INS_3,
                                    SUB_A, SUB_C, SUB_G, SUB_T]):
            self.assertEqual(code, 2 ** exp)

    def test_derived_codes(self):
        """ Test the derived relation codes. """
        self.assertEqual(IRREC, 0)
        self.assertEqual(MINS5, 5)
        self.assertEqual(MINS3, 9)
        self.assertEqual(INS_8, 12)
        self.assertEqual(ANY_8, 13)
        self.assertEqual(INDEL, 14)
        self.assertEqual(SUB_N, 240)
        self.assertEqual(ANY_N, 241)
        self.assertEqual(NOCOV, 255)
        self.assertEqual(ANY_B, ANY_N - SUB_A)
        self.assertEqual(ANY_D, ANY_N - SUB_C)
        self.assertEqual(ANY_H, ANY_N - SUB_G)
        self.assertEqual(ANY_V, ANY_N - SUB_T)

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
