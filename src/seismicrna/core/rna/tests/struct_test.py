import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.rna.struct import RNAStructure, calc_wfmi
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
        # Simple structure with some paired and unpaired bases
        region = Region("myref", "ACCGT")
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(1, 5), (2, 4)])
        expected = pd.Series([True, True, False, True, True], region.range)
        self.assertTrue(structure.is_paired.equals(expected))

        # All bases unpaired
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[])
        expected = pd.Series([False, False, False, False, False], region.range)
        self.assertTrue(structure.is_paired.equals(expected))

        # Longer sequence with nested pairs
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

        # Pseudoknot structure
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

    def test_is_paired_internally(self):
        # Simple structure with two internal pairs (3,8) and (4,7)
        # in stack (2,9)-(3,8)-(4,7)-(5,6)
        region = Region("myref", "ACCGTACCGT")
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 9), (3, 8), (4, 7), (5, 6)])
        expected = pd.Series([False, False, True, True, False,
                              False, True, True, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # Simple structure with two internal pairs (2,9) and (3,8)
        # in stack (1,10)-(2,9)-(3,8)-(4,7)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(1, 10), (2, 9), (3, 8), (4, 7)])
        expected = pd.Series([False, True, True, False, False,
                              False, False, True, True, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # Simple structure with one internal pair (3,7) in stack
        # (2,8)-(3,7)-(4,6)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 8), (3, 7), (4, 6)])
        expected = pd.Series([False, False, True, False, False,
                              False, True, False, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # No internal pairs (all unpaired)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[])
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # No internal pairs (only two consecutive pairs)
        structure = RNAStructure(title="mystructure",
                                 region=region,
                                 pairs=[(2, 8), (3, 7)])
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # Multiple internal pairs in a longer sequence
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]
        )
        expected = pd.Series([False, True, True, True, False,
                              False, True, True, True, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # Pseudoknot structure with no internal pairs
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)]
        )
        expected = pd.Series([False, False, False, False, False,
                              False, False, False, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))

        # Pseudoknot structure with no internal pairs: [((({])))}
        structure = RNAStructure(
            title="mystructure",
            region=region,
            pairs=[(1, 6), (2, 9), (3, 8), (4, 7), (5, 10)]
        )
        expected = pd.Series([False, False, True, False, False,
                              False, False, True, False, False],
                             region.range)
        self.assertTrue(structure.is_paired_internally.equals(expected))


class TestCalcWfmi(ut.TestCase):
    def test_identical_structures(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertEqual(wfmi, 1.0)

    def test_completely_different_structures(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertEqual(wfmi, 0.0)

    def test_all_unpaired(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region, pairs=[])
        struct2 = RNAStructure(title="struct2", region=region, pairs=[])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertEqual(wfmi, 1.0)

    def test_one_empty_one_paired(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region, pairs=[])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertEqual(wfmi, 0.0)

    def test_internal_pairs_only(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(2, 9), (3, 8), (4, 7), (5, 6)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(2, 9), (3, 8), (4, 7), (5, 6)])
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=False)
        self.assertEqual(wfmi, 1.0)

    def test_different_regions(self):
        region1 = Region("myref", "ACCGT")
        region2 = Region("myref", "GCTAC")
        struct1 = RNAStructure(title="struct1", region=region1,
                               pairs=[(1, 5), (2, 4)])
        struct2 = RNAStructure(title="struct2", region=region2,
                               pairs=[(1, 5), (2, 4)])
        with self.assertRaises(ValueError):
            calc_wfmi(struct1, struct2)

    def test_mixed_paired_unpaired(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 10), (2, 9), (4, 7)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertAlmostEqual(wfmi, 0.2 + 0.8 * 2/3, places=3)

    def test_mixed_paired_unpaired_subset(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 10), (2, 9)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertAlmostEqual(wfmi, 0.4 + 0.6 * (2/3)**0.5, places=3)

    def test_pseudoknot_structures(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 5), (2, 6), (3, 7), (4, 8), (9, 10)])
        wfmi = calc_wfmi(struct1, struct2)
        self.assertEqual(wfmi, 1.0)

    def test_all_vs_internal_pairs_1(self):
        region = Region("myref", "ACCGTACCGT")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(3, 8), (4, 7), (5, 6)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(2, 9), (3, 8), (4, 7)])
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=True)
        self.assertAlmostEqual(wfmi, 2/10 + 8/10 * 2/3, places=3)
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=False)
        self.assertAlmostEqual(wfmi, 2/6, places=3)

    def test_all_vs_internal_pairs_2(self):
        region = Region("myref", "ACCGTACCGTAA")
        struct1 = RNAStructure(title="struct1", region=region,
                               pairs=[(2, 9), (3, 8), (4, 7), (5, 6)])
        struct2 = RNAStructure(title="struct2", region=region,
                               pairs=[(1, 10), (2, 9), (3, 8), (4, 7)])
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=True)
        self.assertAlmostEqual(wfmi, 2/12 + 10/12 * 3/4, places=3)
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=False)
        self.assertAlmostEqual(wfmi, 2/8 + 6/8 * 1/2, places=3)

    def test_empty_region(self):
        region = Region("myref", "")
        struct1 = RNAStructure(title="struct1", region=region, pairs=[])
        struct2 = RNAStructure(title="struct2", region=region, pairs=[])
        wfmi = calc_wfmi(struct1, struct2, terminal_pairs=True)
        self.assertTrue(np.isnan(wfmi))


if __name__ == "__main__":
    ut.main()
