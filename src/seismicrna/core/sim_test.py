import unittest as ut

from .sim import rand_dna
from .seq import DNA


class TestRandDna(ut.TestCase):
    """ Test function `rand_dna`. """

    def test_type(self):
        """ Test that the type of the return value is DNA. """
        self.assertIs(type(rand_dna(1)), DNA)

    def test_length(self):
        """ Test that the length of the DNA sequence is as expected. """
        for length in range(1, 10):
            with self.subTest(length=length):
                self.assertEqual(len(rand_dna(length)), length)

    def test_invalid_length(self):
        """ Test that lengths ≤ 0 raise ValueError. """
        for length in range(0, -10, -1):
            with self.subTest(length=length):
                self.assertRaises(ValueError, rand_dna, length)
