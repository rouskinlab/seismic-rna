import unittest as ut

from ..seq import DNA
from ..xam import as_sam


class TestAsSam(ut.TestCase):
    """ Test function `as_sam`. """

    def test_line_in_sam_format(self):
        line = as_sam("FS10000136:77:BPG61616-2310:1:1101:1000:1300", 99,
                      "SARS2_FSE", 1, 42, "151=", "=", 133, 283,
                      DNA("CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTG"
                          "GAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTC"
                          "AGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC"),
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:")
        expect = ("FS10000136:77:BPG61616-2310:1:1101:1000:1300\t99\tSARS2_FSE"
                  "\t1\t42\t151=\t=\t133\t283\t"
                  "CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTT"
                  "ATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCG"
                  "TTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC\t"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:\n")
        self.assertEqual(line, expect)

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
