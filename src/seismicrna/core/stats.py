import numpy as np

from .array import ensure_same_length, get_length


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
    if (n := get_length(alpha, "alpha")) < 2:
        raise ValueError(f"Must have at least 2 alpha parameters, but got {n}")
    if alpha.min() <= 0.:
        raise ValueError(f"All alpha parameters must be > 0, but got {alpha}")
    concentration = alpha.sum()
    mean = alpha / concentration
    variance = (mean * (1. - mean)) / (concentration + 1.)
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
    if (n := ensure_same_length(mean, variance, "mean", "variance")) < 2:
        raise ValueError(f"Must have at least 2 alpha parameters, but got {n}")
    if mean.min() <= 0.:
        raise ValueError(f"Every mean must be > 0, but got {mean}")
    if mean.max() >= 1.:
        raise ValueError(f"Every mean must be < 0, but got {mean}")
    if not np.isclose((mean_sum := mean.sum()), 1.):
        raise ValueError(f"All means must sum to 1, but got {mean_sum}")
    if variance.min() <= 0.:
        raise ValueError(f"Every variance must be > 0, but got {variance}")
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
