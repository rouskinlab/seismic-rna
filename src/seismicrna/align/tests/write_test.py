import unittest as ut

from seismicrna.align.write import calc_flags
from seismicrna.core.ngs import FLAG_FIRST, FLAG_SECOND, FLAG_REVERSE


class TestCalcFlags(ut.TestCase):

    def test_f1r2_paired(self):
        expect = (([FLAG_FIRST, FLAG_SECOND | FLAG_REVERSE],
                   [FLAG_SECOND | FLAG_REVERSE, FLAG_FIRST]),
                  ([FLAG_FIRST | FLAG_REVERSE, FLAG_SECOND],
                   [FLAG_SECOND, FLAG_FIRST | FLAG_REVERSE]))
        result = calc_flags(True, True)
        self.assertEqual(result, expect)

    def test_f1r2_single(self):
        expect = (([0],
                   [FLAG_REVERSE]),
                  ([FLAG_REVERSE],
                   [0]))
        result = calc_flags(True, False)
        self.assertEqual(result, expect)

    def test_f2r1_paired(self):
        expect = (([FLAG_FIRST | FLAG_REVERSE, FLAG_SECOND],
                   [FLAG_SECOND, FLAG_FIRST | FLAG_REVERSE]),
                  ([FLAG_FIRST, FLAG_SECOND | FLAG_REVERSE],
                   [FLAG_SECOND | FLAG_REVERSE, FLAG_FIRST]))
        result = calc_flags(False, True)
        self.assertEqual(result, expect)

    def test_f2r1_single(self):
        expect = (([FLAG_REVERSE],
                   [0]),
                  ([0],
                   [FLAG_REVERSE]))
        result = calc_flags(False, False)
        self.assertEqual(result, expect)


if __name__ == "__main__":
    ut.main()
