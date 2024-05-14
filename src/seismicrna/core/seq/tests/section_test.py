import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.seq.section import (BASE_NAME,
                                         FULL_NAME,
                                         POS_NAME,
                                         SEQ_INDEX_NAMES,
                                         Section,
                                         get_shared_index,
                                         iter_windows,
                                         hyphenate_ends,
                                         index_to_pos,
                                         index_to_seq,
                                         intersect,
                                         unite,
                                         seq_pos_to_index,
                                         window_to_margins)
from seismicrna.core.seq.xna import DNA


class TestConstants(ut.TestCase):

    def test_full_name(self):
        self.assertEqual(FULL_NAME, "full")

    def test_pos_name(self):
        self.assertEqual(POS_NAME, "Position")

    def test_base_name(self):
        self.assertEqual(BASE_NAME, "Base")

    def test_seq_index_names(self):
        self.assertEqual(SEQ_INDEX_NAMES, (POS_NAME, BASE_NAME))


class TestHyphenateEnds(ut.TestCase):
    """ Test function `hyphenate_ends`. """

    def test_hyphenate(self):
        self.assertEqual(hyphenate_ends(8, 91), "8-91")


class TestSectionInit(ut.TestCase):

    def test_full(self):
        """ Create a section from a full reference sequence. """
        for length in [0, 20, 100]:
            with self.subTest(length=length):
                seq = DNA.random(length)
                section = Section("myref", seq)
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, FULL_NAME)
                self.assertTrue(section.full)

    def test_full_end5(self):
        """ Create a section with a hyphenated name from a full
        reference sequence. """
        length = 20
        seq = DNA.random(length)
        section = Section("myref", seq, end5=1)
        self.assertEqual(section.ref, "myref")
        self.assertEqual(section.seq, seq)
        self.assertEqual(section.end5, 1)
        self.assertEqual(section.end3, length)
        self.assertEqual(section.name, f"1-{length}")
        self.assertTrue(section.full)

    def test_full_end3(self):
        """ Create a section with a hyphenated name from a full
        reference sequence. """
        length = 20
        seq = DNA.random(length)
        section = Section("myref", seq, end3=length)
        self.assertEqual(section.ref, "myref")
        self.assertEqual(section.seq, seq)
        self.assertEqual(section.end5, 1)
        self.assertEqual(section.end3, length)
        self.assertEqual(section.name, f"1-{length}")
        self.assertTrue(section.full)

    def test_full_blank_name(self):
        """ Create a section with a hyphenated name from a full
        reference sequence. """
        for length in [0, 20, 100]:
            with self.subTest(length=length):
                seq = DNA.random(length)
                section = Section("myref", seq, name="")
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, f"1-{length}")
                self.assertTrue(section.full)

    def test_full_given_name(self):
        """ Create a section with a given name from a full reference
        sequence. """
        for length in [0, 20, 100]:
            with self.subTest(length=length):
                seq = DNA.random(length)
                section = Section("myref", seq, name="mysection")
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, "mysection")
                self.assertTrue(section.full)

    def test_slice_end5_greater(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [20, 100]:
            seq = DNA.random(length)
            for end5 in range(2, length + 1):
                with self.subTest(length=length, end5=end5):
                    section = Section("myref", seq, end5=end5)
                    self.assertEqual(section.ref, "myref")
                    self.assertEqual(section.seq, seq[end5 - 1:])
                    self.assertEqual(section.end5, end5)
                    self.assertEqual(section.end3, length)
                    self.assertEqual(section.name, f"{end5}-{length}")
                    self.assertFalse(section.full)

    def test_slice_end5_equal(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [0, 20, 100]:
            seq = DNA.random(length)
            with self.subTest(length=length):
                section = Section("myref", seq, end5=1)
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, f"1-{length}")
                self.assertTrue(section.full)

    def test_slice_end5_less(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [0, 20, 100]:
            seq = DNA.random(length)
            for end5 in range(-length, 1):
                with self.subTest(length=length, end5=end5):
                    self.assertRaisesRegex(ValueError,
                                           "Need end5 ≥ seq5",
                                           Section,
                                           "myref",
                                           seq,
                                           end5=end5)

    def test_slice_end3_less(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [20, 100]:
            seq = DNA.random(length)
            for end3 in range(1, length):
                with self.subTest(length=length, end3=end3):
                    section = Section("myref", seq, end3=end3)
                    self.assertEqual(section.ref, "myref")
                    self.assertEqual(section.seq, seq[: end3])
                    self.assertEqual(section.end5, 1)
                    self.assertEqual(section.end3, end3)
                    self.assertEqual(section.name, f"1-{end3}")
                    self.assertFalse(section.full)

    def test_slice_end3_equal(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [0, 20, 100]:
            seq = DNA.random(length)
            with self.subTest(length=length):
                section = Section("myref", seq, end3=length)
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, f"1-{length}")
                self.assertTrue(section.full)

    def test_slice_end3_greater(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [20, 100]:
            seq = DNA.random(length)
            for end3 in range(length + 1, 2 * length + 1):
                with self.subTest(length=length, end3=end3):
                    self.assertRaisesRegex(ValueError,
                                           "Need end3 ≤ seq3",
                                           Section,
                                           "myref",
                                           seq,
                                           end3=end3)

    def test_slice_end5_end3(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [20, 100]:
            seq = DNA.random(length)
            for end5 in range(2, length + 1):
                for end3 in range(end5 - 1, length):
                    with self.subTest(length=length, end5=end5, end3=end3):
                        section = Section("myref", seq, end5=end5, end3=end3)
                        self.assertEqual(section.ref, "myref")
                        self.assertEqual(section.seq, seq[end5 - 1: end3])
                        self.assertEqual(section.end5, end5)
                        self.assertEqual(section.end3, end3)
                        self.assertEqual(section.name, f"{end5}-{end3}")
                        self.assertFalse(section.full)

    def test_slice_end5_end3_invalid(self):
        """ Create a section from a slice of a reference sequence. """
        for length in [20, 100]:
            seq = DNA.random(length)
            for end5 in range(2, length + 1):
                for end3 in range(1, end5 - 1):
                    with self.subTest(length=length, end5=end5, end3=end3):
                        self.assertRaisesRegex(ValueError,
                                               "Need end5 ≤ end3",
                                               Section,
                                               "myref",
                                               seq,
                                               end5=end5,
                                               end3=end3)

    def test_partial_seq5_greater(self):
        """ Create a section from part of a reference sequence. """
        for length in [20, 40]:
            seq = DNA.random(length)
            for seq5 in range(2, 2 * length + 1):
                with self.subTest(length=length, seq5=seq5):
                    section = Section("myref", seq, seq5=seq5)
                    self.assertEqual(section.ref, "myref")
                    self.assertEqual(section.seq, seq)
                    self.assertEqual(section.end5, seq5)
                    self.assertEqual(section.end3, seq5 + length - 1)
                    self.assertEqual(section.name,
                                     f"{seq5}-{seq5 + length - 1}")
                    self.assertFalse(section.full)

    def test_partial_seq5_equal(self):
        """ Create a section from part of a reference sequence. """
        for length in [0, 20, 40]:
            seq = DNA.random(length)
            with self.subTest(length=length):
                section = Section("myref", seq, seq5=1)
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, FULL_NAME)
                self.assertTrue(section.full)

    def test_partial_seq5_less(self):
        """ Create a section from part of a reference sequence. """
        for length in [0, 20, 40]:
            seq = DNA.random(length)
            for seq5 in range(-length, 1):
                with self.subTest(length=length, seq5=seq5):
                    self.assertRaisesRegex(ValueError,
                                           "seq5 must be ≥ 1",
                                           Section,
                                           "myref",
                                           seq,
                                           seq5=seq5)

    def test_partial_reflen_greater(self):
        """ Create a section from part of a reference sequence. """
        for length in [20, 40]:
            seq = DNA.random(length)
            for reflen in range(length + 1, 2 * length + 1):
                with self.subTest(length=length, reflen=reflen):
                    section = Section("myref", seq, reflen=reflen)
                    self.assertEqual(section.ref, "myref")
                    self.assertEqual(section.seq, seq)
                    self.assertEqual(section.end5, 1)
                    self.assertEqual(section.end3, length)
                    self.assertEqual(section.name, f"1-{length}")
                    self.assertFalse(section.full)

    def test_partial_reflen_equal(self):
        """ Create a section from part of a reference sequence. """
        for length in [20, 40]:
            seq = DNA.random(length)
            with self.subTest(length=length):
                section = Section("myref", seq, reflen=length)
                self.assertEqual(section.ref, "myref")
                self.assertEqual(section.seq, seq)
                self.assertEqual(section.end5, 1)
                self.assertEqual(section.end3, length)
                self.assertEqual(section.name, FULL_NAME)
                self.assertTrue(section.full)

    def test_partial_reflen_less(self):
        """ Create a section from part of a reference sequence. """
        for length in [20, 40]:
            seq = DNA.random(length)
            for reflen in range(0, length):
                with self.subTest(length=length, reflen=reflen):
                    self.assertRaisesRegex(ValueError,
                                           ("The 3' end of the given sequence "
                                            "is [0-9]+, but the full reference "
                                            "is only [0-9]+ nt"),
                                           Section,
                                           "myref",
                                           seq,
                                           reflen=reflen)

    def test_partial_seq5_reflen(self):
        """ Create a section from part of a reference sequence. """
        for length in [20, 40]:
            seq = DNA.random(length)
            for seq5 in range(2, 2 * length + 1):
                for reflen in range(length + seq5 - 1, 2 * length + 1):
                    with self.subTest(length=length, seq5=seq5, reflen=reflen):
                        section = Section("ref", seq, seq5=seq5, reflen=reflen)
                        self.assertEqual(section.ref, "ref")
                        self.assertEqual(section.seq, seq)
                        self.assertEqual(section.end5, seq5)
                        self.assertEqual(section.end3, seq5 + length - 1)
                        self.assertEqual(section.name,
                                         f"{seq5}-{seq5 + length - 1}")
                        self.assertFalse(section.full)

    def test_partial_slice(self):
        """ Section from part of a slice of a reference sequence. """
        seq = DNA("ACGCAGAAACGGCACTTTCTCACGAGAAGTTGAGCTGAGAATGTGTGGCGGGCGGTCTA")
        seq5 = 38
        rlen = 107
        end5 = 45
        end3 = 78
        sect = Section("ref", seq, seq5=seq5, reflen=rlen, end5=end5, end3=end3)
        self.assertEqual(sect.ref, "ref")
        self.assertEqual(sect.seq, DNA("AACGGCACTTTCTCACGAGAAGTTGAGCTGAGAA"))
        self.assertEqual(sect.end5, end5)
        self.assertEqual(sect.end3, end3)
        self.assertEqual(sect.name, f"{end5}-{end3}")
        self.assertFalse(sect.full)

    def test_partial_slice_invalid_end5(self):
        """ Section from part of a slice of a reference sequence. """
        seq = DNA("ACGCAGAAACGGCACTTTCTCACGAGAAGTTGAGCTGAGAATGTGTGGCGGGCGGTCTA")
        seq5 = 38
        rlen = 107
        end5 = 37
        end3 = 78
        self.assertRaisesRegex(ValueError,
                               "Need end5 ≥ seq5, but got 37 < 38",
                               Section,
                               "myref",
                               seq,
                               seq5=seq5,
                               reflen=rlen,
                               end5=end5,
                               end3=end3)

    def test_partial_slice_invalid_end3(self):
        """ Section from part of a slice of a reference sequence. """
        seq = DNA("ACGCAGAAACGGCACTTTCTCACGAGAAGTTGAGCTGAGAATGTGTGGCGGGCGGTCTA")
        seq5 = 38
        rlen = 107
        end5 = 45
        end3 = 97
        self.assertRaisesRegex(ValueError,
                               "Need end3 ≤ seq3, but got 97 > 96",
                               Section,
                               "myref",
                               seq,
                               seq5=seq5,
                               reflen=rlen,
                               end5=end5,
                               end3=end3)

    def test_partial_slice_invalid_seq5(self):
        """ Section from part of a slice of a reference sequence. """
        seq = DNA("ACGCAGAAACGGCACTTTCTCACGAGAAGTTGAGCTGAGAATGTGTGGCGGGCGGTCTA")
        seq5 = 0
        rlen = 107
        end5 = 45
        end3 = 78
        self.assertRaisesRegex(ValueError,
                               "seq5 must be ≥ 1, but got 0",
                               Section,
                               "myref",
                               seq,
                               seq5=seq5,
                               reflen=rlen,
                               end5=end5,
                               end3=end3)

    def test_partial_slice_invalid_reflen(self):
        """ Section from part of a slice of a reference sequence. """
        seq = DNA("ACGCAGAAACGGCACTTTCTCACGAGAAGTTGAGCTGAGAATGTGTGGCGGGCGGTCTA")
        seq5 = 38
        rlen = 95
        end5 = 45
        end3 = 78
        self.assertRaisesRegex(ValueError,
                               ("The 3' end of the given sequence is 96, "
                                "but the full reference is only 95 nt"),
                               Section,
                               "myref",
                               seq,
                               seq5=seq5,
                               reflen=rlen,
                               end5=end5,
                               end3=end3)


class TestSectionEqual(ut.TestCase):

    def test_equal_full(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        self.assertEqual(section1, section2)

    def test_equal_part(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        self.assertEqual(section1, section2)

    def test_equal_mask(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section1.add_mask("mask", [9])
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section2.add_mask("mask", [9])
        self.assertEqual(section1, section2)

    def test_diff_ref(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref1",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        section2 = Section(ref="ref2",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_seq(self):
        section1 = Section(ref="ref",
                           seq=DNA.random(30),
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=DNA.random(30),
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_seq5(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=2,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=3,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_end5(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=15,
                           end3=32,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_end3(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=33,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_name(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection1")
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           reflen=40,
                           end5=14,
                           end3=32,
                           name="mysection2")
        self.assertNotEqual(section1, section2)

    def test_diff_full(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           reflen=40,
                           name="mysection")
        section2 = Section(ref="ref",
                           seq=seq,
                           name="mysection")
        self.assertNotEqual(section1, section2)

    def test_diff_mask_name(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section1.add_mask("mask1", [9])
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section2.add_mask("mask2", [9])
        self.assertNotEqual(section1, section2)

    def test_diff_mask_positions(self):
        seq = DNA.random(30)
        section1 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section1.add_mask("mask", [9])
        section2 = Section(ref="ref",
                           seq=seq,
                           seq5=7,
                           name="mysection")
        section2.add_mask("mask", [10])
        self.assertNotEqual(section1, section2)


class TestSectionLength(ut.TestCase):

    def test_full_lengths(self):
        for length in [0, 1, 10, 100]:
            with self.subTest(length=length):
                seq = DNA.random(length)
                section = Section("myref", seq)
                self.assertEqual(section.length, length)
                self.assertEqual(section.length, len(section.seq))

    def test_slice_lengths(self):
        for length in [0, 1, 10, 100]:
            seq = DNA.random(length)
            for end5 in range(1, length + 1):
                for end3 in range(end5 - 1, length + 1):
                    with self.subTest(length=length, end5=end5, end3=end3):
                        section = Section("myref", seq, end5=end5, end3=end3)
                        self.assertEqual(section.length, end3 - end5 + 1)
                        self.assertEqual(section.length, len(section.seq))

    def test_partial_lengths(self):
        for length in [0, 1, 10, 100]:
            seq = DNA.random(length)
            for seq5 in range(1, length + 1):
                for rlen in range(length + seq5 - 1, 2 * length + seq5 - 1):
                    with self.subTest(length=length, seq5=seq5, reflen=rlen):
                        section = Section("myref", seq, seq5=seq5, reflen=rlen)
                        self.assertEqual(section.length, length)
                        self.assertEqual(section.length, len(section.seq))


class TestSectionRange(ut.TestCase):

    def test_range(self):
        seq = DNA("ACATGGGGGGACACGTAACTGTTTTCAGAGGAATTGTAGGCATGGATTAACATTCCGTG")
        seq5 = 19
        end5 = 26
        end3 = 77
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        expect = pd.MultiIndex.from_arrays([list(range(end5, end3 + 1)),
                                            list(str(seq[end5 - seq5:
                                                         end3 - seq5 + 1]))],
                                           names=[POS_NAME, BASE_NAME])
        self.assertIsInstance(section.range, pd.MultiIndex)
        self.assertTrue(section.range.equals(expect))

    def test_range_int(self):
        seq = DNA("ACATGGGGGGACACGTAACTGTTTTCAGAGGAATTGTAGGCATGGATTAACATTCCGTG")
        seq5 = 19
        end5 = 26
        end3 = 77
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        expect = np.array(list(range(end5, end3 + 1)))
        self.assertIsInstance(section.range_int, np.ndarray)
        self.assertTrue(section.range_int.dtype is np.dtype('int64'))
        self.assertTrue(np.array_equal(section.range_int, expect))

    def test_range_one(self):
        seq = DNA("ACATGGGGGGACACGTAACTGTTTTCAGAGGAATTGTAGGCATGGATTAACATTCCGTG")
        seq5 = 19
        end5 = 26
        end3 = 77
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        expect = np.array(list(range(1, end3 - end5 + 2)))
        self.assertIsInstance(section.range_one, np.ndarray)
        self.assertIs(section.range_one.dtype, np.dtype('int64'))
        self.assertTrue(np.array_equal(section.range_one, expect))


class TestSectionCopy(ut.TestCase):

    def test_copy(self):
        for length in [2, 10, 100]:
            for reflen in [length, length + 1]:
                for end5 in [1, 2]:
                    for name in [None, "section"]:
                        for masks in [{}, {"mask": [end5]}]:
                            section = Section(ref="ref",
                                              seq=DNA.random(length),
                                              reflen=reflen,
                                              end5=end5,
                                              name=name)
                            for mask, values in masks.items():
                                section.add_mask(mask, values)
                            copied = section.copy()
                            self.assertIsNot(copied, section)
                            self.assertEqual(copied, section)


class TestSectionAddMask(ut.TestCase):

    def test_add_mask(self):
        seq = DNA("GAACCGTGACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        self.assertEqual(list(section._masks.keys()), list())
        section.add_mask("mymask", [18, 20, 25])
        self.assertEqual(list(section._masks.keys()), ["mymask"])
        mymask = section._masks["mymask"]
        self.assertIsInstance(mymask, np.ndarray)
        self.assertIs(mymask.dtype, np.dtype('int64'))
        self.assertTrue(np.array_equal(mymask, np.array([18, 20, 25])))

    def test_add_unordered_mask(self):
        seq = DNA("GAACCGTGACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask", [20, 25, 18, 20])
        self.assertEqual(list(section._masks.keys()), ["mymask"])
        self.assertTrue(np.array_equal(section._masks["mymask"],
                                       np.array([18, 20, 25])))

    def test_add_mask_invert(self):
        seq = DNA("GAACCGTGACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask", [18, 20, 25], complement=True)
        self.assertEqual(list(section._masks.keys()), ["mymask"])
        self.assertTrue(np.array_equal(section._masks["mymask"],
                                       np.array([17, 19, 21, 22, 23, 24, 26])))

    def test_add_overlapping_masks(self):
        seq = DNA("GAACCGTGACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [18, 20, 25])
        section.add_mask("mymask2", [20, 22, 25])
        self.assertEqual(list(section._masks.keys()), ["mymask1", "mymask2"])
        self.assertTrue(np.array_equal(section._masks["mymask1"],
                                       np.array([18, 20, 25])))
        self.assertTrue(np.array_equal(section._masks["mymask2"],
                                       np.array([22])))

    def test_add_mask_invalid(self):
        seq = DNA("GAACCGTGACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        self.assertRaisesRegex(ValueError,
                               "Got positions to mask outside",
                               section.add_mask,
                               "mymask",
                               [16, 20, 25])


class TestSectionMaskNames(ut.TestCase):

    def test_mask_names(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        section = Section("myref", seq)
        names = list()
        for i in range(5):
            name = f"mask{i}"
            names.append(name)
            masked = np.random.choice(section.range_int, size=3, replace=False)
            section.add_mask(name, masked)
        self.assertEqual(section.mask_names, names)

    def test_get_mask(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        section = Section("myref", seq)
        masks = dict()
        for i in range(5):
            name = f"mask{i}"
            masked = np.random.choice(section.unmasked_int, size=3, replace=False)
            masks[name] = np.sort(masked)
            section.add_mask(name, masked)
        for name, masked in masks.items():
            self.assertTrue(np.array_equal(section.get_mask(name), masked))


class TestSectionMaskList(ut.TestCase):

    def test_mask_list(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.mask_list([19, 24, 23, 19])
        self.assertTrue(np.array_equal(section._masks[section.MASK_LIST],
                                       np.array([19, 23, 24])))


class TestSectionMaskGU(ut.TestCase):

    def test_find_gu(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        gu_pos = section._find_gu()
        self.assertIsInstance(gu_pos, np.ndarray)
        self.assertIs(gu_pos.dtype, np.dtype('int64'))
        self.assertTrue(np.array_equal(gu_pos, np.array([18, 19, 21, 23, 25])))

    def test_mask_gu(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.mask_gu()
        self.assertTrue(np.array_equal(section._masks[section.MASK_GU],
                                       section._find_gu()))


class TestSectionMaskPolyA(ut.TestCase):

    def test_find_polya(self):
        seq = DNA("AAAATCGTTCAAGAAAAACAAACA")
        seq5 = 8
        section = Section("myref", seq, seq5=seq5)
        expects = {0: np.array([], dtype=int),
                   1: np.array([8, 9, 10, 11, 18, 19, 21, 22, 23, 24, 25, 27,
                                28, 29, 31]),
                   2: np.array([8, 9, 10, 11, 18, 19, 21, 22, 23, 24, 25, 27,
                                28, 29]),
                   3: np.array([8, 9, 10, 11, 21, 22, 23, 24, 25, 27, 28, 29]),
                   4: np.array([8, 9, 10, 11, 21, 22, 23, 24, 25]),
                   5: np.array([21, 22, 23, 24, 25])}
        for min_length, expect in expects.items():
            with self.subTest(min_length=min_length):
                self.assertTrue(np.array_equal(section._find_polya(min_length),
                                               expect))

    def test_mask_polya(self):
        seq = DNA("AAAATCGTTCAAGAAAAACAAACA")
        seq5 = 8
        for min_length in range(6):
            section = Section("myref", seq, seq5=seq5)
            section.mask_polya(min_length)
            self.assertTrue(np.array_equal(section._masks[section.MASK_POLYA],
                                           section._find_polya(min_length)))


class TestSectionMasked(ut.TestCase):

    def test_masked_bool(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1], dtype=bool)
        self.assertTrue(np.array_equal(section.masked_bool, expect))

    def test_masked_int(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([49, 52, 54, 55, 58])
        self.assertTrue(np.array_equal(section.masked_int, expect))

    def test_masked_zero(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([3, 6, 8, 9, 12])
        self.assertTrue(np.array_equal(section.masked_zero, expect))


class TestSectionUnmasked(ut.TestCase):

    def test_unmasked_bool(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0], dtype=bool)
        self.assertTrue(np.array_equal(section.unmasked_bool, expect))

    def test_unmasked_int(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([46, 47, 48, 50, 51, 53, 56, 57])
        self.assertTrue(np.array_equal(section.unmasked_int, expect))

    def test_unmasked_zero(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = np.array([0, 1, 2, 4, 5, 7, 10, 11])
        self.assertTrue(np.array_equal(section.unmasked_zero, expect))

    def test_unmasked(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        expect = pd.MultiIndex.from_arrays([[46, 47, 48, 50, 51, 53, 56, 57],
                                            list("CCGCCACT")],
                                           names=[POS_NAME, BASE_NAME])
        self.assertTrue(np.array_equal(section.unmasked, expect))

    def test_size(self):
        seq = DNA("CCCGCATCCCGACCAACACTAAGA")
        seq5 = 38
        end5 = 46
        end3 = 58
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.add_mask("mymask1", [49, 54, 58])
        section.add_mask("mymask2", [55, 49, 52])
        self.assertEqual(section.size, section.unmasked_int.size)


class TestSubsection(ut.TestCase):

    def test_full_section_full_sub(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        subsection = section.subsection()
        self.assertIsNot(subsection, section)
        self.assertEqual(subsection, section)

    def test_full_section_full_sub_name(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        subsection = section.subsection(name="mysubsection")
        self.assertEqual(subsection.seq, seq)
        self.assertEqual(subsection.name, "mysubsection")

    def test_full_section_full_length_sub(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        subsection = section.subsection(end5=1, end3=len(seq))
        self.assertEqual(subsection.seq, seq)
        self.assertEqual(subsection.name, f"1-{len(seq)}")

    def test_full_section_full_sub_end5(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        end5 = 10
        subsection = section.subsection(end5=end5)
        self.assertEqual(subsection.seq, seq[end5 - 1:])
        self.assertEqual(subsection.end5, end5)
        self.assertEqual(subsection.end3, len(seq))
        self.assertEqual(subsection.name, f"{end5}-{len(seq)}")

    def test_full_section_full_sub_end3(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        end3 = 23
        subsection = section.subsection(end3=end3)
        self.assertEqual(subsection.seq, seq[: end3])
        self.assertEqual(subsection.end5, 1)
        self.assertEqual(subsection.end3, end3)
        self.assertEqual(subsection.name, f"1-{end3}")

    def test_trunc_section_full_sub(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        end5, end3 = 4, 29
        section = Section("myref", seq, end5=end5, end3=end3)
        subsection = section.subsection()
        self.assertIsNot(subsection, section)
        self.assertEqual(subsection, section)

    def test_trunc_section_trunc_sub(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq, end5=4, end3=29)
        end5, end3 = 9, 21
        subsection = section.subsection(end5=end5, end3=end3)
        self.assertEqual(subsection.seq, seq[end5 - 1: end3])
        self.assertEqual(subsection.end5, end5)
        self.assertEqual(subsection.end3, end3)
        self.assertEqual(subsection.name, f"{end5}-{end3}")

    def test_partial_trunc_section_trunc_sub(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        seq5 = 3
        section = Section("myref", seq, seq5=seq5, end5=4, end3=29)
        end5, end3 = 9, 21
        subsection = section.subsection(end5=end5, end3=end3)
        self.assertEqual(subsection.seq, seq[end5 - seq5: end3 - seq5 + 1])
        self.assertEqual(subsection.end5, end5)
        self.assertEqual(subsection.end3, end3)
        self.assertEqual(subsection.name, f"{end5}-{end3}")


class TestIntersect(ut.TestCase):

    def test_empty_invalid(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot intersect zero sections",
                               intersect)

    def test_one_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        inter = intersect(section)
        self.assertIsNot(inter, section)
        self.assertEqual(inter, section)

    def test_one_full_named(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        inter = intersect(section, name="myintersection")
        self.assertEqual(inter.seq, section.seq)
        self.assertEqual(inter.end5, section.end5)
        self.assertEqual(inter.end3, section.end3)
        self.assertTrue(inter.full)
        self.assertEqual(inter.name, "myintersection")

    def test_one_slice(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq, end5=13, end3=33)
        inter = intersect(section)
        self.assertIsNot(inter, section)
        self.assertEqual(inter, section)

    def test_two_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq)
        section2 = Section("myref", seq)
        inter = intersect(section1, section2)
        self.assertIsNot(inter, section1)
        self.assertIsNot(inter, section2)
        self.assertEqual(inter, section1)
        self.assertEqual(inter, section2)

    def test_two_overlapping(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        end5, end3 = 15, 30
        section1 = Section("myref", seq, end5=end5)
        section2 = Section("myref", seq, end3=end3)
        inter = intersect(section1, section2)
        self.assertEqual(inter.seq, seq[end5 - 1: end3])
        self.assertEqual(inter.end5, end5)
        self.assertEqual(inter.end3, end3)
        self.assertEqual(inter.name, f"{end5}-{end3}")

    def test_two_disjoint(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=6, end3=14)
        section2 = Section("myref", seq, end5=19, end3=35)
        inter = intersect(section1, section2)
        self.assertEqual(inter.seq, DNA(""))
        self.assertEqual(inter.end5, 19)
        self.assertEqual(inter.end3, 18)
        self.assertEqual(inter.name, f"19-18")

    def test_three_overlapping(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        end5, end3 = 15, 30
        section1 = Section("myref", seq, end5=end5)
        section2 = Section("myref", seq, end3=end3)
        section3 = Section("myref", seq, end5=(end5 - 2), end3=(end3 + 2))
        inter = intersect(section1, section2, section3)
        self.assertEqual(inter.seq, seq[end5 - 1: end3])
        self.assertEqual(inter.end5, end5)
        self.assertEqual(inter.end3, end3)
        self.assertEqual(inter.name, f"{end5}-{end3}")

    def test_one_masked(self):
        seq = DNA.random(30)
        section = Section("myref", seq)
        section.add_mask("mask1", [19, 13, 26])
        inter = intersect(section)
        self.assertEqual(inter.masked_int.tolist(), [13, 19, 26])

    def test_two_masked(self):
        seq = DNA.random(30)
        section1 = Section("myref", seq, end5=3, end3=21)
        section2 = Section("myref", seq, end5=9, end3=28)
        section1.add_mask("mask1", [5, 6, 7, 8, 9, 10, 20])
        section2.add_mask("mask2", [10, 15, 20, 25])
        inter = intersect(section1, section2)
        self.assertEqual(inter.end5, 9)
        self.assertEqual(inter.end3, 21)
        self.assertEqual(inter.masked_int.tolist(), [9, 10, 15, 20])

    def test_diff_refs(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCT")
        section1 = Section("ref1", seq)
        section2 = Section("ref2", seq)
        self.assertRaisesRegex(ValueError,
                               "Expected exactly one reference",
                               intersect,
                               section1,
                               section2)

    def test_diff_seqs(self):
        seq1 = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCT")
        seq2 = DNA("CATTCTCGTAGTCAACTATCGCGTCTCTCT")
        section1 = Section("myref", seq1)
        section2 = Section("myref", seq2)
        self.assertRaisesRegex(ValueError,
                               "Expected exactly one sequence",
                               intersect,
                               section1,
                               section2)
        section1 = Section("myref", seq1, end3=10)
        section2 = Section("myref", seq2, end3=10)
        self.assertEqual(intersect(section1, section2).seq,
                         DNA("CATTCTCGTA"))


class TestUnite(ut.TestCase):

    def test_empty_invalid(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot unite zero sections",
                               unite)

    def test_one_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        union = unite(section)
        self.assertIsNot(union, section)
        self.assertEqual(union, section)

    def test_one_full_named(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        union = unite(section, name="myunion")
        self.assertEqual(union.seq, section.seq)
        self.assertEqual(union.end5, section.end5)
        self.assertEqual(union.end3, section.end3)
        self.assertTrue(union.full)
        self.assertEqual(union.name, "myunion")

    def test_one_slice(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq, end5=13, end3=33)
        union = unite(section)
        self.assertIsNot(union, section)
        self.assertEqual(union, section)

    def test_two_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq)
        section2 = Section("myref", seq)
        union = unite(section1, section2)
        self.assertIsNot(union, section1)
        self.assertIsNot(union, section2)
        self.assertEqual(union, section1)
        self.assertEqual(union, section2)

    def test_two_overlapping(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=12, end3=31)
        section2 = Section("myref", seq, end5=4, end3=23)
        union = unite(section1, section2)
        self.assertEqual(union.seq, seq[4 - 1: 31])
        self.assertEqual(union.end5, 4)
        self.assertEqual(union.end3, 31)
        self.assertEqual(union.name, "4-31")
        self.assertEqual(union.unmasked_int.tolist(),
                         list(range(4, 31 + 1)))

    def test_two_overlapping_refseq(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=12, end3=31)
        section2 = Section("myref", seq, end5=4, end3=23)
        self.assertEqual(unite(section1, section2, refseq=seq),
                         unite(section1, section2))

    def test_two_disjoint(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=23, end3=31)
        section2 = Section("myref", seq, end5=4, end3=12)
        union = unite(section1, section2)
        self.assertEqual(union.seq, DNA("TCTCGTAGTNNNNNNNNNNGTCTCTCTA"))
        self.assertEqual(union.end5, 4)
        self.assertEqual(union.end3, 31)
        self.assertEqual(union.name, f"4-31")
        self.assertEqual(union.unmasked_int.tolist(),
                         list(range(4, 12 + 1)) + list(range(23, 31 + 1)))

    def test_two_disjoint_refseq(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=23, end3=31)
        section2 = Section("myref", seq, end5=4, end3=12)
        union = unite(section1, section2, refseq=seq)
        self.assertEqual(union.seq, seq[4 - 1: 31])
        self.assertEqual(union.end5, 4)
        self.assertEqual(union.end3, 31)
        self.assertEqual(union.name, f"4-31")
        self.assertEqual(union.unmasked_int.tolist(),
                         list(range(4, 12 + 1)) + list(range(23, 31 + 1)))

    def test_two_disjoint_wrong_refseq(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=23, end3=31)
        section2 = Section("myref", seq, end5=4, end3=12)
        union = unite(section1,
                      section2,
                      refseq=DNA("CATTCTCGTAGTCAACTATCGCGTCTCTCTACTCTT"))
        self.assertEqual(union.seq, DNA("TCTCGTAGTCAACTATCGCGTCTCTCTA"))
        self.assertRaisesRegex(ValueError,
                               "Expected exactly one sequence",
                               unite,
                               section1,
                               section2,
                               refseq=DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTGCTCTT"))

    def test_one_masked(self):
        seq = DNA.random(30)
        section = Section("myref", seq)
        section.add_mask("mask1", [19, 13, 26])
        union = unite(section)
        self.assertEqual(union.masked_int.tolist(), [13, 19, 26])

    def test_two_masked(self):
        seq = DNA.random(30)
        section1 = Section("myref", seq, end5=3, end3=21)
        section2 = Section("myref", seq, end5=9, end3=28)
        section1.add_mask("mask1", [5, 6, 7, 8, 9, 10, 20])
        section2.add_mask("mask2", [10, 15, 20, 25])
        union = unite(section1, section2)
        self.assertEqual(union.end5, 3)
        self.assertEqual(union.end3, 28)
        self.assertEqual(union.masked_int.tolist(), [5, 6, 7, 8, 10, 20, 25])

    def test_diff_refs(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCT")
        section1 = Section("ref1", seq)
        section2 = Section("ref2", seq)
        self.assertRaisesRegex(ValueError,
                               "Expected exactly one reference",
                               unite,
                               section1,
                               section2)

    def test_diff_seqs(self):
        seq1 = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCT")
        seq2 = DNA("CATTCTCGTAGTCAACTATCGCGTCTCTCT")
        section1 = Section("myref", seq1)
        section2 = Section("myref", seq2)
        self.assertRaisesRegex(ValueError,
                               "Sequences differ",
                               unite,
                               section1,
                               section2)
        section1 = Section("myref", seq1, end3=10)
        section2 = Section("myref", seq2, end3=10)
        self.assertEqual(unite(section1, section2).seq,
                         DNA("CATTCTCGTA"))


class TestIndexToPos(ut.TestCase):
    """ Test function `index_to_pos`. """

    def test_valid_full(self):
        """ Test with a valid full sequence. """
        seq = DNA("ACGT")
        start = 1
        pos = list(range(start, len(seq) + start))
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))

    def test_valid_slice(self):
        """ Test with a valid slice of a sequence. """
        seq = DNA("ACAGCCTAG")
        pos = list(range(7, 11 + 1))
        start = 6
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))

    def test_valid_noncontig(self):
        """ Test with non-contiguous sequence. """
        seq = DNA("ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        result = index_to_pos(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, np.ndarray))
        self.assertTrue(np.array_equal(pos, result))


class TestIndexToSeq(ut.TestCase):
    """ Test function `index_to_seq`. """

    def test_valid_full(self):
        """ Test with a valid full sequence. """
        seq = DNA("ACGT")
        start = 1
        pos = list(range(start, len(seq) + start))
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(seq, result)

    def test_valid_no_pos(self):
        """ Test with a valid sequence and no positions. """
        seq = DNA("ACGT")
        pos = []
        start = 1
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(DNA(""), result)

    def test_valid_empty_seq(self):
        """ Test with an empty sequence and no positions. """
        seq = DNA("")
        pos = []
        start = 1
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(seq, result)

    def test_valid_slice(self):
        """ Test with a valid slice of a sequence. """
        seq = DNA("ACAGCCTAG")
        pos = list(range(7, 11 + 1))
        start = 6
        result = index_to_seq(seq_pos_to_index(seq, pos, start))
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(DNA("CAGCC"), result)

    def test_valid_noncontig(self):
        """ Test with non-contiguous sequence, allowing gaps. """
        seq = DNA("ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        result = index_to_seq(seq_pos_to_index(seq, pos, start),
                              allow_gaps=True)
        self.assertTrue(isinstance(result, DNA))
        self.assertEqual(DNA("AGCA"), result)

    def test_invalid_noncontig(self):
        """ Test with non-contiguous sequence, forbidding gaps. """
        seq = DNA("ACAGCCTAG")
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
        seq = DNA("ACGT")
        start = 1
        pos = list(range(start, len(seq) + start))
        expected = pd.MultiIndex.from_arrays([pos, ['A', 'C', 'G', 'T']])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_full_9(self):
        """ Test with every position in the sequence, starting at 9. """
        seq = DNA("ACGT")
        start = 9
        pos = list(range(start, len(seq) + start))
        expected = pd.MultiIndex.from_arrays([pos, ['A', 'C', 'G', 'T']])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_slice_6(self):
        """ Test with a slice of the sequence, starting at 6. """
        seq = DNA("ACAGCCTAG")
        pos = list(range(7, 11 + 1))
        start = 6
        expected = pd.MultiIndex.from_arrays([pos, ['C', 'A', 'G', 'C', 'C']])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_noncontig_2(self):
        """ Test with non-contiguous sequence, starting at 2. """
        seq = DNA("ACAGCCTAG")
        pos = [4, 5, 7, 9]
        start = 2
        expected = pd.MultiIndex.from_arrays([pos, ['A', 'G', 'C', 'A']])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_empty_1(self):
        """ Test with no positions, starting at 1. """
        seq = DNA("ACGT")
        pos = []
        start = 1
        expected = pd.MultiIndex.from_arrays([[], []])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_valid_empty_seq(self):
        """ Test with an empty sequence and no positions. """
        seq = DNA("")
        pos = []
        start = 1
        expected = pd.MultiIndex.from_arrays([[], []])
        self.assertTrue(expected.equals(seq_pos_to_index(seq, pos, start)))

    def test_invalid_empty_seq_1(self):
        """ Test with an empty sequence and start position 0. """
        seq = DNA("")
        pos = [0]
        start = 0
        self.assertRaisesRegex(ValueError,
                               "The start position must be ≥ 1, but got 0",
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_empty_seq_2(self):
        """ Test with an empty sequence and position 0. """
        seq = DNA("")
        pos = [0]
        start = 1
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≥ start .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_empty_seq_3(self):
        """ Test with an empty sequence and one position. """
        seq = DNA("")
        pos = [1]
        start = 1
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≤ end .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_full_0(self):
        """ Test with every position in the sequence, starting at 0. """
        seq = DNA("ACGT")
        pos = list(range(1, len(seq) + 1))
        start = 0
        self.assertRaisesRegex(ValueError,
                               "The start position must be ≥ 1, but got 0",
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_less_start_2(self):
        """ Test with a position less than start (= 2). """
        seq = DNA("ACGT")
        pos = list(range(1, len(seq) + 1))
        start = 2
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≥ start .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_greater_end_9(self):
        """ Test with a position greater than end, starting at 9. """
        seq = DNA("ACGT")
        pos = [10, 11, 12, 13]
        start = 9
        self.assertRaisesRegex(ValueError,
                               ("All positions must be ≤ end .*, "
                                "but got .*"),
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_dup_1(self):
        """ Test with duplicated positions, starting at 1. """
        seq = DNA("ACGT")
        pos = [1, 2, 2, 4]
        start = 1
        self.assertRaisesRegex(ValueError,
                               "Duplicated positions: .*",
                               seq_pos_to_index,
                               seq, pos, start)

    def test_invalid_unsort_1(self):
        """ Test with unsorted positions, starting at 1. """
        seq = DNA("ACGT")
        pos = [4, 3, 2, 1]
        start = 1
        self.assertRaisesRegex(ValueError,
                               "Unsorted positions: .*",
                               seq_pos_to_index,
                               seq, pos, start)


class TestGetSharedIndex(ut.TestCase):
    """ Test function `get_shared_index`. """

    def test_zero_valid(self):
        index = get_shared_index([], empty_ok=True)
        self.assertIsInstance(index, pd.MultiIndex)
        self.assertEqual(tuple(index.names), SEQ_INDEX_NAMES)

    def test_zero_invalid(self):
        self.assertRaisesRegex(ValueError,
                               "No indexes were given",
                               get_shared_index,
                               [])

    def test_one_valid(self):
        index = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                           ["A", "C", "G", "T", "N"]],
                                          names=SEQ_INDEX_NAMES)
        self.assertIs(get_shared_index([index]), index)

    def test_one_invalid(self):
        index = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                           ["A", "C", "G", "T", "N"]],
                                          names=[BASE_NAME, POS_NAME])
        self.assertRaisesRegex(ValueError,
                               "Expected levels named",
                               get_shared_index,
                               [index])

    def test_two_valid(self):
        index0 = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                            ["A", "C", "G", "T", "N"]],
                                           names=SEQ_INDEX_NAMES)
        index1 = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                            ["A", "C", "G", "T", "N"]],
                                           names=SEQ_INDEX_NAMES)
        self.assertIs(get_shared_index([index0, index1]), index0)

    def test_two_invalid(self):
        index0 = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                            ["A", "C", "G", "T", "N"]],
                                           names=SEQ_INDEX_NAMES)
        index1 = pd.MultiIndex.from_arrays([[4, 6, 10, 11, 12],
                                            ["A", "C", "G", "T", "A"]],
                                           names=SEQ_INDEX_NAMES)
        self.assertRaisesRegex(ValueError,
                               "Indexes 0 and 1 differ",
                               get_shared_index,
                               [index0, index1])


class TestWindowToMargins(ut.TestCase):
    """ Test function `window_to_margins`. """

    def test_window_margins(self):
        expected = {1: (0, 0),
                    2: (0, 1),
                    3: (1, 1),
                    4: (1, 2),
                    5: (2, 2),
                    6: (2, 3),
                    7: (3, 3)}
        for window, margins in expected.items():
            self.assertEqual(window_to_margins(window), margins)

    def test_window_length(self):
        for window in range(1, 101):
            self.assertEqual(sum(window_to_margins(window)) + 1, window)

    def test_window_0(self):
        self.assertRaisesRegex(ValueError,
                               "window must be ≥ 1, but got 0",
                               window_to_margins,
                               0)

    def test_window_float(self):
        self.assertRaisesRegex(TypeError,
                               "window must be int, but got float",
                               window_to_margins,
                               1.)


class TestGetWindows(ut.TestCase):
    """ Test function `get_windows`. """

    def compare(self, expected: list, *series: pd.Series, **kwargs):
        windows = list(iter_windows(*series, **kwargs))
        for (wcenter, wseries), (ecenter, eseries) in zip(windows,
                                                          expected,
                                                          strict=True):
            self.assertEqual(wcenter, ecenter)
            for window, expect in zip(wseries, eseries, strict=True):
                self.assertTrue(np.allclose(window, expect, equal_nan=True))

    def test_empty(self):
        self.assertEqual(list(iter_windows(size=1)), [])

    def test_1_series_size_1_min_1_excl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3.],)),
                    (2, ([1.],)),
                    (4, ([6.],)),
                    (6, ([8.],)),
                    (8, ([4.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=1)

    def test_1_series_size_1_min_1_incl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3.],)),
                    (2, ([1.],)),
                    (4, ([6.],)),
                    (5, ([np.nan],)),
                    (6, ([8.],)),
                    (7, ([np.nan],)),
                    (8, ([4.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=1, include_nan=True)

    def test_1_series_size_1_min_0_excl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3.],)),
                    (2, ([1.],)),
                    (4, ([6.],)),
                    (5, ([],)),
                    (6, ([8.],)),
                    (7, ([],)),
                    (8, ([4.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=1, min_count=0)

    def test_1_series_size_1_min_0_incl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3.],)),
                    (2, ([1.],)),
                    (4, ([6.],)),
                    (5, ([np.nan],)),
                    (6, ([8.],)),
                    (7, ([np.nan],)),
                    (8, ([4.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=1, min_count=0, include_nan=True)

    def test_1_series_size_2_min_1_excl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3., 1.],)),
                    (2, ([1.],)),
                    (4, ([6.],)),
                    (5, ([8.],)),
                    (6, ([8.],)),
                    (7, ([4.],)),
                    (8, ([4., 5.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=2)

    def test_1_series_size_2_min_1_incl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3., 1.],)),
                    (2, ([1.],)),
                    (4, ([6., np.nan],)),
                    (5, ([np.nan, 8.],)),
                    (6, ([8., np.nan],)),
                    (7, ([np.nan, 4.],)),
                    (8, ([4., 5.],)),
                    (9, ([5.],))]
        self.compare(expected, series, size=2, include_nan=True)

    def test_1_series_size_2_min_2_excl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3., 1.],)),
                    (8, ([4., 5.],))]
        self.compare(expected, series, size=2, min_count=2)

    def test_1_series_size_2_min_2_incl_nan(self):
        series = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.],
                           seq_pos_to_index(DNA("ACGTNGATC"),
                                            [1, 2, 4, 5, 6, 7, 8, 9],
                                            start=1))
        expected = [(1, ([3., 1.],)),
                    (4, ([6., np.nan],)),
                    (5, ([np.nan, 8.],)),
                    (6, ([8., np.nan],)),
                    (7, ([np.nan, 4.],)),
                    (8, ([4., 5.],))]
        self.compare(expected, series, size=2, min_count=2, include_nan=True)

    def test_2_series_size_1_min_1_excl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3.], [2.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (6, ([8.], [5.])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=1)

    def test_2_series_size_1_min_1_incl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3.], [2.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (5, ([np.nan], [np.nan])),
                    (6, ([8.], [5.])),
                    (7, ([np.nan], [9.])),
                    (8, ([4.], [np.nan])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=1, include_nan=True)

    def test_2_series_size_1_min_0_excl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3.], [2.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (5, ([], [])),
                    (6, ([8.], [5.])),
                    (7, ([], [])),
                    (8, ([], [])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=1, min_count=0)

    def test_2_series_size_1_min_0_incl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3.], [2.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (5, ([np.nan], [np.nan])),
                    (6, ([8.], [5.])),
                    (7, ([np.nan], [9.])),
                    (8, ([4.], [np.nan])),
                    (9, ([5.], [6.]))]
        self.compare(expected,
                     series1,
                     series2,
                     size=1,
                     min_count=0,
                     include_nan=True)

    def test_2_series_size_2_min_1_excl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (5, ([8.], [5.])),
                    (6, ([8.], [5.])),
                    (8, ([5.], [6.])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=2)

    def test_2_series_size_2_min_1_incl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.])),
                    (2, ([1.], [7.])),
                    (4, ([6., np.nan], [0., np.nan])),
                    (5, ([np.nan, 8.], [np.nan, 5.])),
                    (6, ([8., np.nan], [5., 9.])),
                    (7, ([np.nan, 4.], [9., np.nan])),
                    (8, ([4., 5.], [np.nan, 6.])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=2, include_nan=True)

    def test_2_series_size_2_min_2_excl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.]))]
        self.compare(expected, series1, series2, size=2, min_count=2)

    def test_2_series_size_2_min_2_incl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.])),
                    (4, ([6., np.nan], [0., np.nan])),
                    (5, ([np.nan, 8.], [np.nan, 5.])),
                    (6, ([8., np.nan], [5., 9.])),
                    (7, ([np.nan, 4.], [9., np.nan])),
                    (8, ([4., 5.], [np.nan, 6.]))]
        self.compare(expected,
                     series1,
                     series2,
                     size=2,
                     min_count=2,
                     include_nan=True)

    def test_2_series_size_2_min_0_excl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.])),
                    (2, ([1.], [7.])),
                    (4, ([6.], [0.])),
                    (5, ([8.], [5.])),
                    (6, ([8.], [5.])),
                    (7, ([], [])),
                    (8, ([5.], [6.])),
                    (9, ([5.], [6.]))]
        self.compare(expected, series1, series2, size=2, min_count=0)

    def test_2_series_size_2_min_0_incl_nan(self):
        index = seq_pos_to_index(DNA("ACGTNGATC"),
                                 [1, 2, 4, 5, 6, 7, 8, 9],
                                 start=1)
        series1 = pd.Series([3., 1., 6., np.nan, 8., np.nan, 4., 5.], index)
        series2 = pd.Series([2., 7., 0., np.nan, 5., 9., np.nan, 6.], index)
        expected = [(1, ([3., 1.], [2., 7.])),
                    (2, ([1.], [7.])),
                    (4, ([6., np.nan], [0., np.nan])),
                    (5, ([np.nan, 8.], [np.nan, 5.])),
                    (6, ([8., np.nan], [5., 9.])),
                    (7, ([np.nan, 4.], [9., np.nan])),
                    (8, ([4., 5.], [np.nan, 6.])),
                    (9, ([5.], [6.]))]
        self.compare(expected,
                     series1,
                     series2,
                     size=2,
                     min_count=0,
                     include_nan=True)


if __name__ == "__main__":
    ut.main()

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
