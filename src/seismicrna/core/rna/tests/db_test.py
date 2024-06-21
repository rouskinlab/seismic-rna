import unittest as ut

from seismicrna.core.rna.db import (PAIRED_MARKS,
                                    format_db_structure,
                                    parse_db_structure)


class TestPairedMarks(ut.TestCase):

    def test_paired_marks(self):
        self.assertEqual(list(PAIRED_MARKS.items()),
                         list(map(tuple, [")(", "><", "][", "}{"])))


class TestParseDbStructure(ut.TestCase):

    def test_no_pairs(self):
        expect = []
        for length in range(5):
            struct = "." * length
            for seq5 in range(1, 5):
                with self.subTest(length=length, seq5=seq5):
                    self.assertEqual(parse_db_structure(struct, seq5),
                                     expect)

    def test_deep_pairs(self):
        for cmark, omark in PAIRED_MARKS.items():
            for length in range(5):
                struct = omark * length + cmark * length
                for seq5 in range(1, 5):
                    with self.subTest(cmark=cmark,
                                      omark=omark,
                                      length=length,
                                      seq5=seq5):
                        expect = [(seq5 + i, 2 * length + seq5 - (i + 1))
                                  for i in range(length)]
                        self.assertEqual(parse_db_structure(struct, seq5),
                                         expect)

    def test_shallow_pairs(self):
        for cmark, omark in PAIRED_MARKS.items():
            for length in range(5):
                struct = (omark + cmark) * length
                for seq5 in range(1, 5):
                    with self.subTest(cmark=cmark,
                                      omark=omark,
                                      length=length,
                                      seq5=seq5):
                        expect = [(seq5 + 2 * i, seq5 + 2 * i + 1)
                                  for i in range(length)]
                        self.assertEqual(parse_db_structure(struct, seq5),
                                         expect)

    def test_multi_marks(self):
        struct = ".(.[<.{..}>.]{.})."
        expect = [(2, 17), (4, 13), (5, 11), (7, 10), (14, 16)]
        self.assertEqual(parse_db_structure(struct), expect)

    def test_pseudoknot(self):
        struct = ".((<).>)"
        expect = [(2, 8), (3, 5), (4, 7)]
        self.assertEqual(parse_db_structure(struct), expect)

    def test_multi_pseudoknot(self):
        struct = "(<{[)>}]"
        expect = [(1, 5), (2, 6), (3, 7), (4, 8)]
        self.assertEqual(parse_db_structure(struct), expect)

    def test_dangling_opener(self):
        struct = ".{.(.)."
        self.assertRaisesRegex(ValueError,
                               "Position 2 has an unmatched '{'",
                               parse_db_structure,
                               struct)

    def test_dangling_closer(self):
        struct = ".[.].>."
        self.assertRaisesRegex(ValueError,
                               "Position 13 has an unmatched '>'",
                               parse_db_structure,
                               struct, 8)


class TestFormatDbStructure(ut.TestCase):

    def test_no_pairs(self):
        for length in range(5):
            pairs = []
            expect = "." * length
            self.assertEqual(format_db_structure(pairs, length), expect)

    def test_one_pair(self):
        for length in range(11, 20):
            for seq5 in range(1, 5):
                for position in range(seq5, seq5 + 5):
                    for gap in range(3):
                        with self.subTest(length=length,
                                          seq5=seq5,
                                          position=position,
                                          gap=gap):
                            pairs = [(position, position + gap + 1)]
                            expect = "".join([
                                "." * (position - seq5),
                                "(",
                                "." * gap,
                                ")",
                                "." * (length - (position - seq5 + gap + 2))
                            ])
                            self.assertEqual(format_db_structure(pairs,
                                                                 length,
                                                                 seq5),
                                             expect)

    def test_deep_pairs(self):
        for length in range(5):
            expect = "(" * length + ")" * length
            for seq5 in range(1, 5):
                pairs = [(seq5 + i, 2 * length + seq5 - (i + 1))
                         for i in range(length)]
                with self.subTest(length=length, seq5=seq5):
                    self.assertEqual(format_db_structure(pairs,
                                                         length * 2,
                                                         seq5),
                                     expect)

    def test_shallow_pairs(self):
        for length in range(5):
            expect = "()" * length
            for seq5 in range(1, 5):
                with self.subTest(length=length, seq5=seq5):
                    pairs = [(seq5 + 2 * i, seq5 + 2 * i + 1)
                             for i in range(length)]
                    self.assertEqual(format_db_structure(pairs,
                                                         length * 2,
                                                         seq5),
                                     expect)

    def test_pseudoknot(self):
        pairs = [(2, 8), (3, 5), (4, 7)]
        expect = ".((<).>)"
        self.assertEqual(format_db_structure(pairs, len(expect)), expect)

    def test_multi_pseudoknot(self):
        pairs = [(1, 5), (2, 6), (3, 7), (4, 8)]
        expect = "(<[{)>]}"
        self.assertEqual(format_db_structure(pairs, len(expect)), expect)

    def test_excess_pseudoknot(self):
        pairs = [(1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]
        self.assertRaisesRegex(ValueError,
                               "Cannot write base-pair",
                               format_db_structure,
                               pairs, 10)

    def test_invalid_pos5(self):
        pairs = [(2, 8), (0, 5), (4, 7)]
        for seq5 in range(1, 5):
            self.assertRaisesRegex(ValueError,
                                   f"5' partner must be ≥ {seq5}, but got 0",
                                   format_db_structure,
                                   pairs, 10, seq5)

    def test_invalid_pair(self):
        pairs = [(2, 8), (5, 3), (4, 7)]
        self.assertRaisesRegex(ValueError,
                               f"5' partner must be less than 3' partner, "
                               f"but got 5 ≥ 3",
                               format_db_structure,
                               pairs, 10)

    def test_invalid_pos3(self):
        pairs = [(2, 8), (3, 5), (4, 7)]
        self.assertRaisesRegex(ValueError,
                               f"3' partner must be ≤ 7, but got 8",
                               format_db_structure,
                               pairs, 7, 1)
        self.assertRaisesRegex(ValueError,
                               f"3' partner must be ≤ 7, but got 8",
                               format_db_structure,
                               pairs, 6, 2)

    def test_repeat_pos5(self):
        pairs = [(2, 8), (2, 5), (4, 7)]
        self.assertRaisesRegex(ValueError,
                               f"Got >1 base pair for position 2",
                               format_db_structure,
                               pairs, 8)

    def test_repeat_pos3(self):
        pairs = [(2, 8), (3, 4), (4, 7)]
        self.assertRaisesRegex(ValueError,
                               f"Got >1 base pair for position 4",
                               format_db_structure,
                               pairs, 8)


if __name__ == "__main__":
    ut.main()
