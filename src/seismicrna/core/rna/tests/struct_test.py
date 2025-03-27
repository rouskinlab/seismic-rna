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

    def test_is_paired(self):



if __name__ == "__main__":
    ut.main()
