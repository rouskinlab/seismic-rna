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
