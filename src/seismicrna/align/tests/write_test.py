import unittest as ut

from seismicrna.align.write import calc_flags_sep_strands
from seismicrna.core.logs import Level, set_config, restore_config
from seismicrna.core.ngs import (FLAG_PAIRED,
                                 FLAG_PROPER,
                                 FLAG_FIRST,
                                 FLAG_SECOND,
                                 FLAG_REVERSE)


class TestCalcFlags(ut.TestCase):

    @restore_config
    def test_f1r2_paired(self):
        set_config(verbosity=Level.ERROR)
        for bt2_mixed in [False, True]:
            expect = (([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND | FLAG_REVERSE,
                        FLAG_FIRST]),
                      ([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND,
                        FLAG_FIRST | FLAG_REVERSE]))
            result = calc_flags_sep_strands(True, True, bt2_mixed)
            self.assertEqual(result, expect)

    def test_f1r2_single(self):
        expect = (([0],
                   [FLAG_REVERSE | FLAG_PAIRED]),
                  ([FLAG_REVERSE],
                   [FLAG_PAIRED]))
        for mixed in [False, True]:
            result = calc_flags_sep_strands(True, False, mixed)
            self.assertEqual(result, expect)

    @restore_config
    def test_f2r1_paired(self):
        set_config(verbosity=Level.ERROR)
        for bt2_mixed in [False, True]:
            expect = (([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND,
                        FLAG_FIRST | FLAG_REVERSE]),
                      ([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND | FLAG_REVERSE,
                        FLAG_FIRST]))
            result = calc_flags_sep_strands(False, True, bt2_mixed)
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
    ut.main(verbosity=2)
