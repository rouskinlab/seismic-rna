import unittest as ut

from seismicrna.relate.py.cigar import (CIG_ALIGN,
                                        CIG_MATCH,
                                        CIG_SUBST,
                                        CIG_DELET,
                                        CIG_INSRT,
                                        CIG_SCLIP,
                                        parse_cigar)


class TestParseCigar(ut.TestCase):
    """ Test function `parse_cigar`. """

    def test_cigar_match_subst_valid(self):
        """ Parse a valid CIGAR string with match and subst codes. """
        cigar = "9S23=1X13=1D9=2I56=3S"
        expect = [(CIG_SCLIP, 9),
                  (CIG_MATCH, 23),
                  (CIG_SUBST, 1),
                  (CIG_MATCH, 13),
                  (CIG_DELET, 1),
                  (CIG_MATCH, 9),
                  (CIG_INSRT, 2),
                  (CIG_MATCH, 56),
                  (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)

    def test_cigar_align_valid(self):
        """ Parse a valid CIGAR string with align codes. """
        cigar = "9S37M1D9M2I56M3S"
        expect = [(CIG_SCLIP, 9),
                  (CIG_ALIGN, 37),
                  (CIG_DELET, 1),
                  (CIG_ALIGN, 9),
                  (CIG_INSRT, 2),
                  (CIG_ALIGN, 56),
                  (CIG_SCLIP, 3)]
        self.assertEqual(list(parse_cigar(cigar)), expect)


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
