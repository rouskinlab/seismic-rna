import unittest as ut

import numpy as np
from scipy.stats import beta, dirichlet

from seismicrna.core.stats import (calc_beta_mv,
                                   calc_beta_params,
                                   calc_dirichlet_params,
                                   calc_dirichlet_mv,
                                   kumaraswamy_pdf)

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


class TestKumaraswamyPDF(ut.TestCase):

    def test_auc(self):
        x = np.linspace(0, 1, 10001)
        for a in [0.01, 0.1, 1.0, 10.]:
            for b in [0.01, 0.1, 1.0, 10.]:
                with np.errstate(divide="ignore"):
                    expect = ((a * b)
                              * (x ** (a - 1))
                              * ((1 - x ** a) ** (b - 1)))
                    result = kumaraswamy_pdf(x, a, b)
                self.assertTrue(np.allclose(expect, result, equal_nan=True))


if __name__ == "__main__":
    ut.main(verbosity=2)
