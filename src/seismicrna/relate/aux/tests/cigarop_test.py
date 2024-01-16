import string
import unittest as ut

from seismicrna.relate.py.cigar import (CIG_ALIGN,
                                        CIG_MATCH,
                                        CIG_SUBST,
                                        CIG_DELET,
                                        CIG_INSRT,
                                        CIG_SCLIP)
from seismicrna.relate.aux.cigarop import (CigarOp,
                                           count_cigar_muts,
                                           find_cigar_op_pos)


class TestCigarOp(ut.TestCase):
    """ Test class `CigarOp`. """

    ops = CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_DELET, CIG_INSRT, CIG_SCLIP

    def test_cigar_init_valid(self):
        """ Test that CigarOp instances can be created. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            self.assertIsInstance(cigar_op, CigarOp)
            self.assertEqual(cigar_op.op, op)
            self.assertEqual(cigar_op._len, 1)

    def test_cigar_init_invalid(self):
        """ Test that CigarOp instances cannot be created with invalid
        operations. """
        for op in string.printable:
            if op not in self.ops:
                self.assertRaisesRegex(ValueError, f"Invalid CIGAR operation: ",
                                       CigarOp, op)

    def test_cigar_lengthen(self):
        """ Test lengthening CIGAR operations. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            for length in range(1, 10):
                self.assertEqual(cigar_op._len, length)
                cigar_op.lengthen()

    def test_cigar_str(self):
        """ Test string representations of CIGAR operations. """
        for op in self.ops:
            cigar_op = CigarOp(op)
            for length in range(1, 10):
                self.assertEqual(str(cigar_op), f"{length}{op}")
                cigar_op.lengthen()


class TestCountCigarMuts(ut.TestCase):
    """ Test function `count_cigar_muts`. """

    def test_cigar_match_subst_valid(self):
        """ Count mutations in a valid CIGAR string. """
        self.assertEqual(count_cigar_muts("9S23=1X13=1D9=2I56=3S"), 4)


class TestFindCigarOpPos(ut.TestCase):
    """ Test function `find_cigar_op_pos`. """

    def test_cigar_xeq_aln_valid(self):
        """ Find aligned positions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_ALIGN)),
                         [])

    def test_cigar_m_aln_valid(self):
        """ Find aligned positions in a CIGAR string with M codes. """
        self.assertEqual(
            list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M", CIG_ALIGN)),
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
             26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 43, 44,
             45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 61, 62,
             63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
             79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 93, 94, 95,
             96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
             110, 111, 112, 113, 114, 115, 116, 117]
        )

    def test_cigar_xeq_mat_valid(self):
        """ Find matches in a CIGAR string with =/X codes. """
        self.assertEqual(
            list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=", CIG_MATCH)),
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
             26, 27, 28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45,
             46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 61, 62, 63,
             64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
             80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 93, 94, 95, 96,
             97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
             111, 112, 113, 114, 115, 116, 117]
        )

    def test_cigar_m_mat_valid(self):
        """ Find matches in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_MATCH)),
                         [])

    def test_cigar_xeq_sub_valid(self):
        """ Find substitutions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_SUBST)),
                         [33])

    def test_cigar_m_sub_valid(self):
        """ Find substitutions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_SUBST)),
                         [])

    def test_cigar_xeq_del_valid(self):
        """ Find deletions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_DELET)),
                         [])

    def test_cigar_m_del_valid(self):
        """ Find deletions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_DELET)),
                         [])

    def test_cigar_xeq_ins_valid(self):
        """ Find insertions in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_INSRT)),
                         [34, 35, 36, 59, 60, 92])

    def test_cigar_m_ins_valid(self):
        """ Find insertions in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_INSRT)),
                         [34, 35, 36, 59, 60, 92])

    def test_cigar_xeq_scl_valid(self):
        """ Find soft clippings in a CIGAR string with =/X codes. """
        self.assertEqual(list(find_cigar_op_pos("9S23=1X3I13=1D9=2I31=1I25=",
                                                CIG_SCLIP)),
                         [1, 2, 3, 4, 5, 6, 7, 8, 9])

    def test_cigar_m_scl_valid(self):
        """ Find soft clippings in a CIGAR string with M codes. """
        self.assertEqual(list(find_cigar_op_pos("9S24M3I13M1D9M2I31M1I25M",
                                                CIG_SCLIP)),
                         [1, 2, 3, 4, 5, 6, 7, 8, 9])


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
