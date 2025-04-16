import numpy as np

from .array import ensure_same_length, get_length

from .validate import require_greater, require_between


def _validate_alpha(alpha: np.ndarray):
    if (length := get_length(alpha, "alpha")) < 2:
        raise ValueError(
            f"Must have at least 2 alpha parameters, but got {length}"
        )
    if alpha.min() <= 0.:
        raise ValueError(f"All alpha parameters must be > 0, but got {alpha}")
    return length


def _validate_mean_variance(mean: np.ndarray, variance: np.ndarray, n: int = 1):
    if (length := ensure_same_length(mean, variance, "mean", "variance")) < 2:
        raise ValueError(f"Must have at least 2 alpha parameters, but got {n}")
    if mean.min() <= 0.:
        raise ValueError(f"Every mean must be > 0, but got {mean}")
    if mean.max() >= n:
        raise ValueError(f"Every mean must be < {n}, but got {mean}")
    if not np.isclose((mean_sum := mean.sum()), n):
        raise ValueError(f"All means must sum to {n}, but got {mean_sum}")
    if variance.min() <= 0.:
        raise ValueError(f"Every variance must be > 0, but got {variance}")
    return length


def calc_dirichlet_mv(alpha: np.ndarray):
    """ Find the means and variances of a Dirichlet distribution from
    its concentration parameters.

    Parameters
    ----------
    alpha: np.ndarray
        Concentration parameters of the Dirichlet distribution.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Means and variances of the Dirichlet distribution.
    """
    _validate_alpha(alpha)
    conc = alpha.sum()
    mean = alpha / conc
    variance = (mean * (1. - mean)) / (conc + 1.)
    return mean, variance


def calc_dirichlet_params(mean: np.ndarray, variance: np.ndarray):
    """ Find the concentration parameters of a Dirichlet distribution
    from its mean and variance.

    Parameters
    ----------
    mean: np.ndarray
        Means.
    variance: np.ndarray
        Variances.

    Returns
    -------
    np.ndarray
        Concentration parameters.
    """
    _validate_mean_variance(mean, variance)
    concentrations = (mean * (1. - mean)) / variance - 1.
    concentration = concentrations.mean()
    if not np.allclose(concentrations, concentration):
        raise ValueError(f"Incompatible means {mean} and variances {variance}")
    if concentration <= 0.:
        raise ValueError("Every variance must be < mean * (1 - mean), "
                         f"but got {variance} > {mean * (1. - mean)}")
    return concentration * mean


def calc_beta_mv(alpha: float, beta: float):
    """ Find the mean and variance of a beta distribution from its alpha
    and beta parameters.

    Parameters
    ----------
    alpha: float
        Alpha parameter of the beta distribution.
    beta: float
        Beta parameter of the beta distribution.

    Returns
    -------
    tuple[float, float]
        Mean and variance of the beta distribution.
    """
    (mean, _), (variance, _) = calc_dirichlet_mv(np.array([alpha, beta]))
    return float(mean), float(variance)


def calc_beta_params(mean: float, variance: float):
    """ Find the alpha and beta parameters of a beta distribution from
    its mean and variance.

    Parameters
    ----------
    mean: float
        Mean of the beta distribution.
    variance: float
        Variance of the beta distribution.

    Returns
    -------
    tuple[float, float]
        Alpha and beta parameters of the beta distribution.
    """
    alpha, beta = calc_dirichlet_params(np.array([mean, 1. - mean]),
                                        np.array([variance, variance]))
    return float(alpha), float(beta)


def kumaraswamy_pdf(x: np.ndarray, a: float | int, b: float | int):
    """ Kumaraswamy distribution probability density function (PDF).

    Parameters
    ----------
    x: np.ndarray
        Input values; must be in the interval [0, 1].
    a: float | int
        Shape parameter a; must be > 0.
    b: float | int
        Shape parameter b; must be > 0.

    Returns
    -------
    np.ndarray
        Kumaraswamy distribution PDF at input values.
    """
    require_greater("a", a, 0., classes=(float, int))
    require_greater("b", b, 0., classes=(float, int))
    return ((a * b)
            * np.power(x, a - 1.)
            * np.power(1. - np.power(x, a), b - 1.))


def double_kumaraswamy_pdf(x: np.ndarray, 
                           w: float | int, 
                           a1: float | int, 
                           b1: float | int, 
                           a2: float | int, 
                           b2: float | int):
    """ Double Kumaraswamy distribution probability density function
    (PDF).

    Parameters
    ----------
    x: np.ndarray
        Input values; must be in the interval [0, 1].
    w: float | int
        Weight for distribution 1; must be in the interval [0, 1].
    a1: float | int
        Shape parameter a for distribution 1; must be > 0.
    b1: float | int
        Shape parameter b for distribution 1; must be > 0.
    a2: float | int
        Shape parameter a for distribution 2; must be > 0.
    b2: float | int
        Shape parameter b for distribution 2; must be > 0.

    Returns
    -------
    np.ndarray
        Double Kumaraswamy distribution PDF at input values.
    """
    require_between("w", w, 0., 1., classes=(float, int))
    return w * kumaraswamy_pdf(x, a1, b1) + (1 - w) * kumaraswamy_pdf(x, a2, b2)
