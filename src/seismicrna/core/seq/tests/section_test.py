import unittest as ut

import numpy as np
import pandas as pd

from ..section import (BASE_NAME,
                       FULL_NAME,
                       POS_INDEX,
                       POS_NAME,
                       Section,
                       index_to_pos,
                       index_to_seq,
                       intersection,
                       seq_pos_to_index)
from ..xna import DNA


class TestConstants(ut.TestCase):

    def test_full_name(self):
        self.assertEqual(FULL_NAME, "full")

    def test_pos_index(self):
        """ Test that sequence positions are 1-indexed. """
        self.assertEqual(POS_INDEX, 1)

    def test_pos_name(self):
        self.assertEqual(POS_NAME, "Position")

    def test_base_name(self):
        self.assertEqual(BASE_NAME, "Base")


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


class TestSectionLength(ut.TestCase):

    def test_full_lengths(self):
        for length in [0, 20, 100]:
            with self.subTest(length=length):
                seq = DNA.random(length)
                section = Section("myref", seq)
                self.assertEqual(section.length, length)
                self.assertEqual(section.length, len(section.seq))

    def test_slice_lengths(self):
        for length in [0, 20, 100]:
            seq = DNA.random(length)
            for end5 in range(1, length + 1):
                for end3 in range(end5 - 1, length + 1):
                    with self.subTest(length=length, end5=end5, end3=end3):
                        section = Section("myref", seq, end5=end5, end3=end3)
                        self.assertEqual(section.length, end3 - end5 + 1)
                        self.assertEqual(section.length, len(section.seq))

    def test_partial_lengths(self):
        for length in [0, 20, 100]:
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
        section.add_mask("mymask", [18, 20, 25], invert=True)
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


class TestSectionMaskPos(ut.TestCase):

    def test_mask_pos(self):
        seq = DNA("GAACCGTNACGGATCTCGCAATGT")
        seq5 = 8
        end5 = 17
        end3 = 26
        section = Section("myref", seq, seq5=seq5, end5=end5, end3=end3)
        section.mask_pos([19, 24, 23, 19])
        self.assertTrue(np.array_equal(section._masks[section.MASK_POS],
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


class TestIntersection(ut.TestCase):

    def test_empty_invalid(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot intersect zero sections",
                               intersection)

    def test_one_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        inter = intersection(section)
        self.assertIsNot(inter, section)
        self.assertEqual(inter, section)

    def test_one_full_named(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq)
        inter = intersection(section, name="myintersection")
        self.assertEqual(inter.seq, section.seq)
        self.assertEqual(inter.end5, section.end5)
        self.assertEqual(inter.end3, section.end3)
        self.assertTrue(inter.full)
        self.assertEqual(inter.name, f"myintersection")

    def test_one_slice(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section = Section("myref", seq, end5=13, end3=33)
        inter = intersection(section)
        self.assertIsNot(inter, section)
        self.assertEqual(inter, section)

    def test_two_full(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq)
        section2 = Section("myref", seq)
        inter = intersection(section1, section2)
        self.assertIsNot(inter, section1)
        self.assertIsNot(inter, section2)
        self.assertEqual(inter, section1)
        self.assertEqual(inter, section2)

    def test_two_overlapping(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        end5, end3 = 15, 30
        section1 = Section("myref", seq, end5=end5)
        section2 = Section("myref", seq, end3=end3)
        inter = intersection(section1, section2)
        self.assertEqual(inter.seq, seq[end5 - 1: end3])
        self.assertEqual(inter.end5, end5)
        self.assertEqual(inter.end3, end3)
        self.assertEqual(inter.name, f"{end5}-{end3}")

    def test_two_disjoint(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        section1 = Section("myref", seq, end5=6, end3=17)
        section2 = Section("myref", seq, end5=19, end3=35)
        inter = intersection(section1, section2)
        self.assertEqual(inter.seq, DNA(''))
        self.assertEqual(inter.end5, 19)
        self.assertEqual(inter.end3, 18)
        self.assertEqual(inter.name, f"19-18")

    def test_three_overlapping(self):
        seq = DNA("CATTCTCGTAGTCAACTTTCGCGTCTCTCTACTCTT")
        end5, end3 = 15, 30
        section1 = Section("myref", seq, end5=end5)
        section2 = Section("myref", seq, end3=end3)
        section3 = Section("myref", seq, end5=(end5 - 2), end3=(end3 + 2))
        inter = intersection(section1, section2, section3)
        self.assertEqual(inter.seq, seq[end5 - 1: end3])
        self.assertEqual(inter.end5, end5)
        self.assertEqual(inter.end3, end3)
        self.assertEqual(inter.name, f"{end5}-{end3}")


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

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
