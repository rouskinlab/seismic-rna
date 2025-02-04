import unittest as ut

from seismicrna.core.arg.cli import *


class TestMARCD(ut.TestCase):

    def test_default_min_marcd_run_equals_default_max_marcd_join(self):
        default_marcd = 0.0175
        self.assertEqual(opt_min_marcd_run.default, default_marcd)
        self.assertEqual(opt_max_marcd_join.default, default_marcd)


if __name__ == "__main__":
    ut.main(verbosity=2)
