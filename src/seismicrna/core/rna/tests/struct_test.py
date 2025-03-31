import unittest as ut

import pandas as pd

from seismicrna.core.rna.struct import RNAStructure
from seismicrna.core.seq.region import Region


class TestRNAStructure(ut.TestCase):

    def test_init_success(self):
        region = Region("myref", "ACCGT")
        structure = RNAStructure(title="mystructure",
                                 branch="mybranch",
                                 region=region,
                                 pairs=[(1, 5), (2, 4)])
        self.assertIs(structure.region, region)
        self.assertEqual(structure.title, "mystructure")
        self.assertEqual(structure.branch, "mybranch")
        self.assertTrue(structure.table.equals(pd.Series([5, 4, 0, 2, 1],
                                                         region.range)))

    def test_init_outside_region(self):
        region = Region("myref", "ACCGT")
        self.assertRaisesRegex(ValueError,
                               "Position 6 is not in Region",
                               RNAStructure,
                               title="mystructure",
                               branch="mybranch",
                               region=region,
                               pairs=[(1, 6), (2, 4)])

    def test_from_db_string(self):
        struct = RNAStructure.from_db_string(".((...).)",
                                             "ACGATACAG",
                                             seq5=11,
                                             ref="myref",
                                             reg="myreg",
                                             title="mystruct",
                                             branch="mybranch")
        self.assertEqual(struct.title, "mystruct")
        self.assertEqual(struct.branch, "mybranch")
        self.assertEqual(struct.region,
                         Region("myref", "ACGATACAG", seq5=11, name="myreg"))
        self.assertSetEqual(struct.pairs, {(12, 19), (13, 17)})

    def test_is_paired(self):
        # Test case 1: Simple structure with some paired and unpaired bases
        region = Region("myref", "ACCGT")
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(1, 5), (2, 4)])
        expected = pd.Series([True, True, False, True, True], region.range)
        self.assertTrue(structure.is_paired.equals(expected))

        # Test case 2: All bases unpaired
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[])
        expected = pd.Series([False, False, False, False, False], region.range)
        self.assertTrue(structure.is_paired.equals(expected))

        # Test case 4: Longer sequence with nested pairs
        region = Region("myref", "ACCGTACCGT")
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]
        )
        expected = pd.Series([True, True, True, True, True,
                              True, True, True, True, True],
                             region.range)
        self.assertTrue(structure.is_paired.equals(expected))

        # Test case 5: Pseudoknot structure
        region = Region("myref", "ACCGTACCGT")
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)]
        )
        expected = pd.Series([True, True, True, True, True,
                              True, True, True, True, True],
                             region.range)
        self.assertTrue(structure.is_paired.equals(expected))

    def test_is_middle(self):
        # Simple structure with three middle positions (3,8), (4,7) and (5,6) 
        # in stack (2,9)-(3,8)-(4,7)-(5,6)
        region = Region("myref", "ACCGTACCGT")
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 9), (3, 8), (4, 7), (5, 6)])
        expected = pd.Series([False, False, True, True, True,
                              True, True, True, False, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # Simple structure with two middle positions (2,9) and (3,8) 
        # in stack (1,10)-(2,9)-(3,8)-(4,7)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(1, 10), (2, 9), (3, 8), (4, 7)])
        expected = pd.Series([False, True, True, False, False,
                              False, False, True, True, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # Simple structure with one middle position (3,7) in stack 
        # (2,8)-(3,7)-(4,6)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 8), (3, 7), (4, 6)])
        expected = pd.Series([False, False, True, False, False,
                              False, True, False, False, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # No middle positions (all unpaired)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[])
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # No middle positions (only two consecutive pairs)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 8), (3, 7)])
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # Multiple middle positions in a longer sequence
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]
        )
        expected = pd.Series([False, True, True, True, True,
                              True, True, True, True, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))

        # Pseudoknot structure with no middle positions
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)]
        )
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_middle.equals(expected))


if __name__ == "__main__":
    ut.main()
