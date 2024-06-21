import unittest as ut

from seismicrna.core.rna.db import (PAIRED_MARKS,
                                    _parse_db_structure)


class TestPairedMarks(ut.TestCase):

    def test_paired_marks(self):
        self.assertEqual(sorted(PAIRED_MARKS.items()),
                         sorted(map(tuple, [")(", "][", "}{", "><"])))


class TestParseDbStructure(ut.TestCase):

    def test_no_pairs(self):
        expect = []
        for length in range(5):
            struct = "." * length
            for seq5 in range(1, 5):
                self.assertEqual(_parse_db_structure(struct, seq5), expect)

    def test_deep_pairs(self):
        for cmark, omark in PAIRED_MARKS.items():
            for length in range(5):
                struct = omark * length + cmark * length
                for seq5 in range(1, 5):
                    expect = [(seq5 + i, 2 * length + seq5 - (i + 1))
                              for i in range(length)]
                    self.assertEqual(_parse_db_structure(struct, seq5), expect)

    def test_shallow_pairs(self):
        for cmark, omark in PAIRED_MARKS.items():
            for length in range(5):
                struct = (omark + cmark) * length
                for seq5 in range(1, 5):
                    expect = [(seq5 + 2 * i, seq5 + 2 * i + 1)
                              for i in range(length)]
                    self.assertEqual(_parse_db_structure(struct, seq5), expect)

    def test_pseudoknot(self):
        struct = ".((<).>)"
        expect = [(2, 8), (3, 5), (4, 7)]
        self.assertEqual(_parse_db_structure(struct, 1), expect)


if __name__ == "__main__":
    ut.main()
