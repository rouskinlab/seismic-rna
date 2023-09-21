import unittest as ut

from ..sam import _find_blank_range


class TestFindBlankRange(ut.TestCase):
    """ Test function `_find_blank_range`. """

    def test_side5_zero_length(self):
        self.assertEqual(_find_blank_range(False, 0, 1, 10),
                         (10, 10))

    def test_side5_under_length(self):
        self.assertEqual(_find_blank_range(False, 3, 1, 10),
                         (3, 10))

    def test_side5_equal_length(self):
        self.assertEqual(_find_blank_range(False, 10, 1, 10),
                         (10, 10))

    def test_side5_over_length(self):
        self.assertEqual(_find_blank_range(False, 11, 1, 10),
                         (10, 10))

    def test_side5_neg_length(self):
        self.assertRaisesRegex(ValueError, "Length of read must be ≥ 1",
                               _find_blank_range, False, -1, 1, 10)

    def test_side3_zero_length(self):
        self.assertEqual(_find_blank_range(True, 0, 1, 10),
                         (0, 0))

    def test_side3_under_length(self):
        self.assertEqual(_find_blank_range(True, 3, 1, 10),
                         (0, 7))

    def test_side3_equal_length(self):
        self.assertEqual(_find_blank_range(True, 10, 1, 10),
                         (0, 0))

    def test_side3_over_length(self):
        self.assertEqual(_find_blank_range(True, 11, 1, 10),
                         (0, 0))

    def test_side3_neg_length(self):
        self.assertRaisesRegex(ValueError, "Length of read must be ≥ 1",
                               _find_blank_range, True, -1, 1, 10)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                              #
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
