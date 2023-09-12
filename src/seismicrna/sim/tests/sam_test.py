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
