import unittest as ut

import numpy as np
import pandas as pd

from .sect import encode_primer, index_to_pos, index_to_seq, seq_pos_to_index
from .seq import DNA


class TestEncodePrimer(ut.TestCase):
    """ Test function `encode_primer`. """

    def test_valid_primers(self):
        """ Test that valid primers can be encoded. """
        ref = "myref"
        primer1 = "ACTACGACGTGACTAGCT"
        primer2 = "TCTACTCCATTTTCAATACT"
        result = encode_primer(ref, primer1, primer2)
        expected_value = ref, DNA(primer1.encode()), DNA(primer2.encode())
        expected_type = str, DNA, DNA
        self.assertEqual(result, expected_value)
        self.assertEqual(tuple(map(type, result)), expected_type)


class TestIndexToPos(ut.TestCase):
    """ Test function `index_to_pos`. """

    def test_valid_full(self):
        """ Test with a valid full sequence. """
        seq = DNA(b"ACGT")
        pos = [1, 2, 3, 4]
        start = 1
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))

    def test_valid_slice(self):
        """ Test with a valid slice of a sequence. """
        seq = DNA(b"ACAGCCTAG")
        pos = [7, 8, 9, 10, 11]
        start = 6
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))

    def test_valid_noncontig(self):
        """ Test with non-contiguous sequence. """
        seq = DNA(b"ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))


class TestIndexToSeq(ut.TestCase):
    """ Test function `index_to_seq`. """

    def test_valid_full(self):
        """ Test with a valid full sequence. """
        seq = DNA(b"ACGT")
        pos = [1, 2, 3, 4]
        start = 1
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(seq, result)

    def test_valid_slice(self):
        """ Test with a valid slice of a sequence. """
        seq = DNA(b"ACAGCCTAG")
        pos = [7, 8, 9, 10, 11]
        start = 6
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(DNA(b"CAGCC"), result)

    def test_valid_noncontig(self):
        """ Test with non-contiguous sequence, allowing gaps. """
        seq = DNA(b"ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        result = index_to_seq(seq_pos_to_index(seq, pos, start),
                              allow_gaps=True)
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(DNA(b"AGCA"), result)

    def test_invalid_noncontig(self):
        """ Test with non-contiguous sequence, forbidding gaps. """
        seq = DNA(b"ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        self.assertRaisesRegex(ValueError,
                               ("A sequence cannot be assembled from an index "
                                "with missing positions"),
                               index_to_seq,
                               seq_pos_to_index(seq, pos, start))


class TestSeqPosToIndex(ut.TestCase):
    """ Test function `seq_pos_to_index`. """

    def test_valid_full_1(self):
        """ Test with every position in the sequence, starting at 1. """
        seq = DNA(b"ACGT")
        pos = [1, 2, 3, 4]
        start = 1
        expected = pd.MultiIndex.from_arrays([pos, ["A", "C", "G", "T"]])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_full_9(self):
        """ Test with every position in the sequence, starting at 9. """
        seq = DNA(b"ACGT")
        pos = [9, 10, 11, 12]
        start = 9
        expected = pd.MultiIndex.from_arrays([pos, ["A", "C", "G", "T"]])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_slice_6(self):
        """ Test with a slice of the sequence, starting at 6. """
        seq = DNA(b"ACAGCCTAG")
        pos = [7, 8, 9, 10, 11]
        start = 6
        expected = pd.MultiIndex.from_arrays([pos, ["C", "A", "G", "C", "C"]])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_noncontig_2(self):
        """ Test with non-contiguous sequence, starting at 2. """
        seq = DNA(b"ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        expected = pd.MultiIndex.from_arrays([pos, ["A", "G", "C", "A"]])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_empty_1(self):
        """ Test with no positions, starting at 1. """
        seq = DNA(b"ACGT")
        pos = []
        start = 1
        expected = pd.MultiIndex.from_arrays([[], []])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_invalid_full_0(self):
        """ Test with every position in the sequence, starting at 0. """
        seq = DNA(b"ACGT")
        pos = [1, 2, 3, 4]
        start = 0
        self.assertRaisesRegex(ValueError,
                               "The start position must be ≥ 1, but got 0",
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_less_start_2(self):
        """ Test with a position less than start (= 2). """
        seq = DNA(b"ACGT")
        pos = [1, 2, 3, 4]
        start = 2
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≥ start .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_greater_end_9(self):
        """ Test with a position greater than end, starting at 9. """
        seq = DNA(b"ACGT")
        pos = [9, 10, 11, 13]
        start = 9
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≤ end .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_dup_1(self):
        """ Test with duplicated positions, starting at 1. """
        seq = DNA(b"ACGT")
        pos = [1, 2, 2, 4]
        start = 1
        self.assertRaisesRegex(ValueError,
                               "Duplicated positions: .*",
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_unsort_1(self):
        """ Test with unsorted positions, starting at 1. """
        seq = DNA(b"ACGT")
        pos = [4, 3, 2, 1]
        start = 1
        self.assertRaisesRegex(ValueError,
                               "Unsorted positions: .*",
                               seq_pos_to_index,
                               seq, pos, start)
