import unittest as ut

import numpy as np
from scipy.stats import beta, dirichlet

from seismicrna.core.stats import (calc_beta_mv,
                                   calc_beta_params,
                                   calc_dirichlet_params,
                                   calc_dirichlet_mv)

rng = np.random.default_rng()


def rand_dirichlet_alpha(n: int):
    """ Simulate `n` alpha parameters for a Dirichlet distribution. """
    # Alpha parameters can be any positive real number.
    return -np.log(rng.random(n))


class TestCalcDirichletMV(ut.TestCase):

    def test_scipy_dirichlet(self):
        for n in range(2, 6):
            with self.subTest(n=n):
                alpha = rand_dirichlet_alpha(n)
                mean, var = calc_dirichlet_mv(alpha)
                self.assertTrue(np.allclose(mean, dirichlet.mean(alpha)))
                self.assertTrue(np.allclose(var, dirichlet.var(alpha)))


class TestCalcDirichletParams(ut.TestCase):

    def test_invert(self):
        for n in range(2, 6):
            with self.subTest(n=n):
                alpha_true = rand_dirichlet_alpha(n)
                alpha = calc_dirichlet_params(*calc_dirichlet_mv(alpha_true))
                self.assertTrue(np.allclose(alpha, alpha_true))


class TestCalcBetaMV(ut.TestCase):

    def test_scipy_beta(self):
        a, b = rand_dirichlet_alpha(2)
        mean, var = calc_beta_mv(a, b)
        self.assertTrue(np.isclose(mean, beta.mean(a, b)))
        self.assertTrue(np.isclose(var, beta.var(a, b)))


class TestCalcBetaParams(ut.TestCase):

    def test_invert(self):
        a_true, b_true = rand_dirichlet_alpha(2)
        a, b = calc_beta_params(*calc_beta_mv(a_true, b_true))
        self.assertTrue(np.isclose(a, a_true))
        self.assertTrue(np.isclose(b, b_true))


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
