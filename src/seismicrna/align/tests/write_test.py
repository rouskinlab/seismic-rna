import unittest as ut

from seismicrna.align.write import calc_flags_sep_strands
from seismicrna.core.logs import Level, set_config, restore_config
from seismicrna.core.ngs import (FLAG_PAIRED,
                                 FLAG_PROPER,
                                 FLAG_FIRST,
                                 FLAG_SECOND,
                                 FLAG_REVERSE)


class TestCalcFlags(ut.TestCase):

    def test_f1r2_paired_unmixed(self):
        expect = (([FLAG_FIRST | FLAG_PAIRED,
                    FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED],
                   [FLAG_SECOND | FLAG_REVERSE,
                    FLAG_FIRST]),
                  ([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED,
                    FLAG_SECOND | FLAG_PAIRED],
                   [FLAG_SECOND,
                    FLAG_FIRST | FLAG_REVERSE]))
        result = calc_flags_sep_strands(True, True, False)
        self.assertEqual(result, expect)

    @restore_config
    def test_f1r2_paired_mixed(self):
        set_config(verbosity=Level.ERROR)
        expect = (([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                    FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                   [FLAG_SECOND | FLAG_REVERSE,
                    FLAG_FIRST]),
                  ([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                    FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                   [FLAG_SECOND,
                    FLAG_FIRST | FLAG_REVERSE]))
        result = calc_flags_sep_strands(True, True, True)
        self.assertEqual(result, expect)

    def test_f1r2_single(self):
        expect = (([0],
                   [FLAG_REVERSE | FLAG_PAIRED]),
                  ([FLAG_REVERSE],
                   [FLAG_PAIRED]))
        for mixed in [False, True]:
            result = calc_flags_sep_strands(True, False, mixed)
            self.assertEqual(result, expect)

    def test_f2r1_paired_unmixed(self):
        expect = (([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED,
                    FLAG_SECOND | FLAG_PAIRED],
                   [FLAG_SECOND,
                    FLAG_FIRST | FLAG_REVERSE]),
                  ([FLAG_FIRST | FLAG_PAIRED,
                    FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED],
                   [FLAG_SECOND | FLAG_REVERSE,
                    FLAG_FIRST]))
        result = calc_flags_sep_strands(False, True, False)
        self.assertEqual(result, expect)

    @restore_config
    def test_f2r1_paired_mixed(self):
        set_config(verbosity=Level.ERROR)
        expect = (([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                    FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                   [FLAG_SECOND,
                    FLAG_FIRST | FLAG_REVERSE]),
                  ([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                    FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                   [FLAG_SECOND | FLAG_REVERSE,
                    FLAG_FIRST]))
        result = calc_flags_sep_strands(False, True, True)
        self.assertEqual(result, expect)

    def test_f2r1_single(self):
        expect = (([FLAG_REVERSE],
                   [FLAG_PAIRED]),
                  ([0],
                   [FLAG_REVERSE | FLAG_PAIRED]))
        for mixed in [False, True]:
            result = calc_flags_sep_strands(False, False, mixed)
            self.assertEqual(result, expect)


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
