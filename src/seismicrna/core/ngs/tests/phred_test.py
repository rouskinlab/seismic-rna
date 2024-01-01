"""

Tests for the Quality Core Module

========================================================================

"""

import unittest as ut

from ..phred import (LO_QUAL,
                     OK_QUAL,
                     HI_QUAL,
                     LO_PHRED,
                     OK_PHRED,
                     HI_PHRED,
                     decode_phred,
                     encode_phred)


class TestConstants(ut.TestCase):

    def test_quals(self):
        self.assertLess(LO_QUAL, OK_QUAL)
        self.assertLess(OK_QUAL, HI_QUAL)

    def test_phreds(self):
        self.assertLess(LO_PHRED, OK_PHRED)
        self.assertLess(OK_PHRED, HI_PHRED)


class TestDecode(ut.TestCase):

    def test_decode(self):
        self.assertEqual(decode_phred('I', 33), 40)
        self.assertEqual(decode_phred('!', 33), 0)


class TestEncode(ut.TestCase):

    def test_encode(self):
        self.assertEqual(encode_phred(40, 33), 'I')
        self.assertEqual(encode_phred(0, 33), '!')

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
