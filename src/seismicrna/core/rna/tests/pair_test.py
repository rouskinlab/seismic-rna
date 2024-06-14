import unittest as ut

import pandas as pd

from seismicrna.core.rna.pair import (UNPAIRED,
                                      pairs_to_dict,
                                      pairs_to_table,
                                      dict_to_pairs,
                                      dict_to_table,
                                      table_to_dict,
                                      table_to_pairs,
                                      find_enclosing_pairs)
from seismicrna.core.seq import FIELD_END5, FIELD_END3, DNA, Section


class TestConstants(ut.TestCase):

    def test_unpaired(self):
        """ Partner given to unpaired bases in CT format. """
        self.assertEqual(UNPAIRED, 0)


class TestPairsToDict(ut.TestCase):

    def test_empty(self):
        pairs = iter([])
        expect = dict()
        self.assertEqual(pairs_to_dict(pairs), expect)

    def test_normal(self):
        pairs = iter([(2, 5), (7, 4), (10, 13)])
        expect = {2: 5, 4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        self.assertEqual(pairs_to_dict(pairs), expect)

    def test_invalid_position(self):
        pairs = iter([(2, 5), (0, 4), (10, 13)])
        self.assertRaisesRegex(ValueError,
                               "Position must be ≥ 1, but got 0",
                               pairs_to_dict,
                               pairs)

    def test_conflicting_pair_1(self):
        pairs = iter([(2, 5), (7, 4), (10, 7)])
        self.assertRaisesRegex(ValueError,
                               "Position 7 was given pairs with both 4 and 10",
                               pairs_to_dict,
                               pairs)

    def test_conflicting_pair_2(self):
        pairs = iter([(7, 10), (2, 5), (7, 4)])
        self.assertRaisesRegex(ValueError,
                               "Position 7 was given pairs with both 10 and 4",
                               pairs_to_dict,
                               pairs)


class TestDictToPairs(ut.TestCase):

    def test_empty(self):
        pairs = dict()
        expect = list()
        self.assertEqual(list(dict_to_pairs(pairs)), expect)

    def test_normal(self):
        pairs = {2: 5, 4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        expect = [(2, 5), (4, 7), (10, 13)]
        self.assertEqual(list(dict_to_pairs(pairs)), expect)

    def test_missing_reverse_1(self):
        pairs = {2: 5, 4: 7, 7: 4, 10: 13, 13: 10}
        self.assertRaisesRegex(ValueError,
                               r"Pair \(2, 5\) is missing its reverse \(5, 2\)",
                               lambda x: list(dict_to_pairs(x)),
                               pairs)

    def test_missing_reverse_2(self):
        pairs = {4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        self.assertRaisesRegex(ValueError,
                               r"Pair \(5, 2\) is missing its reverse \(2, 5\)",
                               lambda x: list(dict_to_pairs(x)),
                               pairs)

    def test_invalid_position(self):
        pairs = {0: 5, 4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        self.assertRaisesRegex(ValueError,
                               "Position must be ≥ 1, but got 0",
                               lambda x: list(dict_to_pairs(x)),
                               pairs)

    def test_invalid_pair(self):
        pairs = {2: 5, 4: 4, 5: 2, 10: 13, 13: 10}
        self.assertRaisesRegex(ValueError,
                               "Position 4 is paired with itself",
                               lambda x: list(dict_to_pairs(x)),
                               pairs)


class TestPairsToTable(ut.TestCase):

    def test_no_pairs(self):
        pairs = list()
        for n in range(5):
            section = Section("myref", DNA.random(n))
            expect = pd.Series(UNPAIRED, section.range)
            result = pairs_to_table(pairs, section)
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.equals(expect))

    def test_normal(self):
        section = Section("myref", DNA.random(13))
        pairs = iter([(2, 5), (7, 4), (10, 13)])
        expect = pd.Series([0, 5, 0, 7, 2, 0, 4, 0, 0, 13, 0, 0, 10],
                           section.range)
        result = pairs_to_table(pairs, section)
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))

    def test_invalid_position_1(self):
        section = Section("myref", DNA.random(13))
        pairs = iter([(2, 5), (7, 4), (10, 14)])
        self.assertRaisesRegex(ValueError,
                               "Position 14 is not in",
                               pairs_to_table,
                               pairs, section)

    def test_invalid_position_2(self):
        section = Section("myref", DNA.random(13))
        pairs = iter([(2, 5), (7, 4), (14, 10)])
        self.assertRaisesRegex(ValueError,
                               "Position 14 is not in",
                               pairs_to_table,
                               pairs, section)

    def test_conflicting_pair_1(self):
        section = Section("myref", DNA.random(13))
        pairs = iter([(2, 5), (7, 4), (10, 7)])
        self.assertRaisesRegex(ValueError,
                               "Position 7 was given pairs with both 4 and 10",
                               pairs_to_table,
                               pairs, section)

    def test_conflicting_pair_2(self):
        section = Section("myref", DNA.random(13))
        pairs = iter([(7, 10), (2, 5), (7, 4)])
        self.assertRaisesRegex(ValueError,
                               "Position 7 was given pairs with both 10 and 4",
                               pairs_to_table,
                               pairs, section)


class TestDictToTable(ut.TestCase):

    def test_empty(self):
        pairs = dict()
        for n in range(5):
            section = Section("myref", DNA.random(n))
            expect = pd.Series(UNPAIRED, section.range)
            result = dict_to_table(pairs, section)
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.equals(expect))

    def test_normal(self):
        section = Section("myref", DNA.random(13))
        pairs = {2: 5, 4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        expect = pd.Series([0, 5, 0, 7, 2, 0, 4, 0, 0, 13, 0, 0, 10],
                           section.range)
        result = dict_to_table(pairs, section)
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))


class TestTableToDict(ut.TestCase):

    def test_empty(self):
        pairs = dict()
        for n in range(5):
            section = Section("myref", DNA.random(n))
            result = table_to_dict(dict_to_table(pairs, section))
            self.assertEqual(result, pairs)

    def test_normal(self):
        section = Section("myref", DNA.random(13))
        pairs = {2: 5, 4: 7, 5: 2, 7: 4, 10: 13, 13: 10}
        result = table_to_dict(dict_to_table(pairs, section))
        self.assertEqual(result, pairs)


class TestTableToPairs(ut.TestCase):

    def test_empty(self):
        pairs = list()
        for n in range(5):
            section = Section("myref", DNA.random(n))
            result = list(table_to_pairs(pairs_to_table(pairs, section)))
            self.assertEqual(result, pairs)

    def test_normal(self):
        section = Section("myref", DNA.random(13))
        pairs = [(2, 5), (4, 7), (10, 13)]
        result = list(table_to_pairs(pairs_to_table(iter(pairs), section)))
        self.assertEqual(result, pairs)


class TestFindEnclosingPairs(ut.TestCase):

    def test_no_pairs(self):
        for n in range(5):
            section = Section("myref", DNA.random(n))
            pairs = list()
            table = pairs_to_table(pairs, section)
            expect = pd.DataFrame(0, section.range, [FIELD_END5, FIELD_END3])
            result = find_enclosing_pairs(table)
            self.assertIsInstance(result, pd.DataFrame)
            self.assertTrue(result.equals(expect))

    def test_nested(self):
        section = Section("myref", DNA.random(14))
        pairs = [(2, 7), (4, 5), (10, 13)]
        table = pairs_to_table(pairs, section)
        expect = pd.DataFrame([[0, 0],
                               [2, 7],
                               [2, 7],
                               [4, 5],
                               [4, 5],
                               [2, 7],
                               [2, 7],
                               [0, 0],
                               [0, 0],
                               [10, 13],
                               [10, 13],
                               [10, 13],
                               [10, 13],
                               [0, 0]],
                              section.range,
                              [FIELD_END5, FIELD_END3])
        result = find_enclosing_pairs(table)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.equals(expect))

    def test_pseudoknotted(self):
        section = Section("myref", DNA.random(14))
        pairs = [(2, 5), (4, 7), (10, 13)]
        table = pairs_to_table(pairs, section)
        self.assertRaisesRegex(ValueError,
                               r"Pairs \(2, 5\) and \(4, 7\) are not nested",
                               find_enclosing_pairs,
                               table)


if __name__ == "__main__":
    ut.main()
