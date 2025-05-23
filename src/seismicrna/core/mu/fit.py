import numpy as np
import pandas as pd

from .nan import auto_remove_nan
from ..logs import logger
from ..validate import require_between

FPAIRED = "fpaired"
PAIRED_ALPHA = "paired_alpha"
PAIRED_BETA = "paired_beta"
UNPAIRED_ALPHA = "unpaired_alpha"
UNPAIRED_BETA = "unpaired_beta"

# Log-normal distributions of beta distribution parameters:
# - mutation rate paired, alpha
# - mutation rate paired, beta
# - mutation rate unpaired, alpha
# - mutation rate unpaired, beta
DEFAULT_MEAN = {
    "A": np.array([+0.987601, +4.307069, +0.958335, +2.893603]),
    "C": np.array([+0.088802, +4.228603, +0.580111, +2.962363]),
    "G": np.zeros(4, dtype=float),
    "T": np.zeros(4, dtype=float),
}
DEFAULT_COV = {
    "A": np.array([[+0.649940, +0.815070, -0.059535, -0.007325],
                   [+0.815070, +1.272768, -0.232291, -0.164021],
                   [-0.059535, -0.232291, +0.632850, +0.793738],
                   [-0.007325, -0.164021, +0.793738, +1.083976]]),
    "C": np.array([[+0.292252, +0.436732, +0.008293, +0.037922],
                   [+0.436732, +0.927412, -0.000711, +0.095237],
                   [+0.008293, -0.000711, +0.478951, +0.541135],
                   [+0.037922, +0.095237, +0.541135, +0.787187]]),
    "G": np.eye(4, dtype=float),
    "T": np.eye(4, dtype=float),
}


@auto_remove_nan
def fit_beta_mixture_model(mus: np.ndarray | pd.Series,
                           ab_params_mean: np.ndarray,
                           ab_params_cov: np.ndarray,
                           fpaired_alpha: float = 1.,
                           fpaired_beta: float = 1.,
                           eps: float = 1.e-6,
                           n_trials: int = 20,
                           maxiter: int = 10000,
                           ftol: float = 1.e-8,
                           gtol: float = 1.e-8):
    """ Fit a two-component beta mixture model to mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray | pd.Series
        Array of mutation rates (values between 0 and 1)
    ab_params_mean: numpy.ndarray
        Mean of the log-normal distribution of the alpha and beta parameters
    ab_params_cov: numpy.ndarray
        Covariance matrix of the log-normal distribution of the alpha and beta parameters
    fpaired_alpha: float
        Alpha parameter for the beta distribution of the fraction of paired bases
    fpaired_beta: float
        Beta parameter for the beta distribution of the fraction of paired bases
    eps: float
        Clip all mus < eps to eps, and all mus > (1 - eps) to (1 - eps).
    n_trials: int
        Number of optimization attempts with different random initial parameters
    maxiter: int
        Maximum number of iterations
    ftol: float
        Tolerance for the function value
    gtol: float
        Tolerance for the gradient
    
    Returns
    -------
    dict[str, float]
        Best fitting parameters found across all trials
    """
    logger.routine("Began fitting beta mixture model")
    from scipy.stats import beta, multivariate_normal
    from scipy.optimize import OptimizeResult, minimize
    # Beta distributions cannot handle 0 and 1, so clip the mutation
    # rates to eps and (1 - eps).
    require_between("eps", eps, 0., 1., inclusive=False, classes=float)
    mus = np.asarray_chkfinite(np.clip(mus, eps, 1. - eps))
    # Log-normal distributions of the alpha and beta parameters
    ab_params_dist = multivariate_normal(ab_params_mean, ab_params_cov)
    # Beta distribution of the fraction of paired bases.
    fpaired_dist = beta(fpaired_alpha, fpaired_beta)

    def calc_log_posterior(x: np.ndarray):
        """ Log posterior probability of the parameters (x) given the
        mutation rates (mus).

        The posterior probability P(x | mus) is given by:
        P(x | mus) = P(mus | x) * P(x) / P(mus)

        P(mus | x) is the likelihood of mus given x, which is given by:
        P(mus | x) = ∏ P(mu_i | x)
                   = ∏ (P(mu_i | pa, pb) * P(paired_i | mu_i, x))
                        +
                        P(mu_i | ua, ub) * (1 - P(paired_i | mu_i, x)))
        where pa, pb, ua, and ub are the alpha and beta parameters for
        paired and unpaired bases, and fp is the fraction paired.
        Note that taking the product over all mu_i assumes independence
        of the mutation rates, which is not strictly true.

        P(paired_i | mu_i, x) is the probability that base i is paired:
        P(paired_i | mu_i, x) = (fp * P(mu_i | pa, pb)
                                 /
                                 (fp * P(mu_i | pa, pb) 
                                  +
                                  (1 - fp) * P(mu_i | ua, ub)))

        P(x) is the prior probability of x, which is given by:
        P(x) = P(fp) * P(pa, pb, ua, ub)

        P(mus) is the probability of the mutation rates, which does not
        depend on the parameters x and is therefore not relevant to the
        optimization of x, so it is not included in the calculation.

        Parameters
        ----------
        x: numpy.ndarray
            Length-5 array of parameters:
            0.  fp: Fraction of paired bases
            1.  pa: Alpha parameter for paired bases
            2.  pb: Beta parameter for paired bases
            3.  ua: Alpha parameter for unpaired bases
            4.  ub: Beta parameter for unpaired bases

        Returns
        -------
        float
            Log posterior probability of the parameters.
        """
        # Calculate the probability that each base is paired.
        paired_pdf = beta.pdf(mus, x[1], x[2])
        unpaired_pdf = beta.pdf(mus, x[3], x[4])
        paired_joint = x[0] * paired_pdf
        unpaired_joint = (1. - x[0]) * unpaired_pdf
        with np.errstate(invalid="ignore"):
            paired_prob = paired_joint / (paired_joint + unpaired_joint)
        unpaired_prob = 1. - paired_prob
        # Calculate the log likelihood.
        log_likelihood = np.sum(np.log(
            paired_prob * paired_pdf + unpaired_prob * unpaired_pdf
        ))
        # Calculate the prior log probability of the parameters x.
        fp_log_prior = fpaired_dist.logpdf(x[0])
        pa_pb_ua_ub_log_prior = ab_params_dist.logpdf(np.log(x[1:]))
        return log_likelihood + fp_log_prior + pa_pb_ua_ub_log_prior

    def succeeded(r: OptimizeResult):
        """ Check if optimization succeeded. """
        if not isinstance(r, OptimizeResult):
            return False
        if not r.success:
            return False
        fpaired_, pa_, pb_, ua_, ub_ = r.x
        assert eps <= fpaired_ <= 1. - eps
        assert eps <= pa_ <= 1. / eps
        assert eps <= pb_ <= 1. / eps
        assert eps <= ua_ <= 1. / eps
        assert eps <= ub_ <= 1. / eps
        # Unpaired bases should have a higher mean.
        mean_paired = pa_ / (pa_ + pb_)
        mean_unpaired = ua_ / (ua_ + ub_)
        if mean_paired >= mean_unpaired:
            return False
        # The pseudoenergy (proporional to the log of unpaired divided
        # by paired PDF) must increase monotonically. A necessary and
        # sufficient condition is that pa ≤ ua and pb ≥ ub.
        if pa_ > ua_ or pb_ < ub_:
            return False
        return True

    # Run the optimization n_trials times with different random initial
    # parameters.
    best_result = None
    best_log_posterior = -np.inf
    for i in range(n_trials):
        # Generate random initial parameters.
        fpaired = fpaired_dist.rvs()
        ab_params = np.exp(ab_params_dist.rvs())
        x0 = np.concatenate([[fpaired], ab_params])
        # Run the optimization.
        result = minimize(lambda x: -calc_log_posterior(x),
                          x0,
                          bounds=[(eps, 1. - eps)] + [(eps, 1. / eps)] * 4,
                          method="L-BFGS-B",
                          options={"maxiter": maxiter,
                                   "ftol": ftol,
                                   "gtol": gtol})
        if succeeded(result):
            # If the optimization was successful, then check if the log
            # posterior is the best so far.
            log_posterior = calc_log_posterior(result.x)
            logger.detail(f"Trial {i} succeeded: {log_posterior}; {result.x}")
            if log_posterior > best_log_posterior:
                best_log_posterior = log_posterior
                best_result = result.x
        else:
            logger.detail(f"Trial {i} failed")

    if best_result is None:
        raise RuntimeError("Failed to fit parameters in any trial")

    fpaired, pa, pb, ua, ub = map(float, best_result)
    params = {FPAIRED: fpaired,
              PAIRED_ALPHA: pa,
              PAIRED_BETA: pb,
              UNPAIRED_ALPHA: ua,
              UNPAIRED_BETA: ub}
    logger.detail(f"Fit beta mixture model parameters: {params}")
    logger.routine("Ended fitting beta mixture model")
    plot_beta_mixture(mus, params)
    return params


def plot_beta_mixture(mus, params):
    """
    Plot the fitted beta mixture model against the data histogram.

    Parameters:
    -----------
    mus : numpy.ndarray
        The data (mutation rates)
    weights : list or numpy.ndarray
        The weights of the mixture components [weight1, weight2]
    params : list or numpy.ndarray
        The parameters of the beta distributions [alpha1, beta1, alpha2, beta2]
    """
    import matplotlib.pyplot as plt
    from scipy.stats import beta

    weight1, alpha1, beta1, alpha2, beta2 = params.values()
    weight2 = 1. - weight1

    # Create grid for plotting
    x_grid = np.linspace(0, 1, 10000)

    # Compute component densities and mixture density
    pdf1 = beta.pdf(x_grid, alpha1, beta1)
    pdf2 = beta.pdf(x_grid, alpha2, beta2)
    mixture_pdf = weight1 * pdf1 + weight2 * pdf2
    energy = -0.61597 * np.log(pdf1 / pdf2)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot histogram of data
    plt.hist(mus, bins=30, density=True, alpha=0.5, label='Data')

    # Plot component densities and mixture density
    plt.plot(x_grid, weight1 * pdf1, 'r--', label=f'Component 1 (α={alpha1:.2f}, β={beta1:.2f}, p={weight1:.2f})')
    plt.plot(x_grid, weight2 * pdf2, 'g--', label=f'Component 2 (α={alpha2:.2f}, β={beta2:.2f}, p={weight2:.2f})')
    plt.plot(x_grid, mixture_pdf, 'b-', label='Mixture')
    plt.plot(x_grid, energy, "k-", label="Pseudoenergy")

    # Add annotations
    mean1 = alpha1 / (alpha1 + beta1)
    mean2 = alpha2 / (alpha2 + beta2)
    plt.axvline(x=mean1, color='r', linestyle=':', alpha=0.7)
    plt.axvline(x=mean2, color='g', linestyle=':', alpha=0.7)

    plt.title('Beta Mixture Model Fit')
    plt.xlabel('Mutation Rate')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()
