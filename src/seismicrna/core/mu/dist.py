from typing import Callable

import numpy as np
import pandas as pd

from ..seq import BASE_NAME, BASEA, BASEC, Region
from ..stats import kumaraswamy_pdf, double_kumaraswamy_pdf
from ..validate import require_isinstance, require_atleast, require_between


# Gas constant
R = 1.98720425864083e-3  # kcal/(mol*K)


# Reference parameters for each base.
REF_PARAMS = {
    BASEA: ((1.30298999, 131.5804808),
            (0.557323482, 1.388569622, 11.16626786, 2.778348856, 618.6712217)),
    BASEC: ((0.774939875, 44.15809053),
            (0.374710709, 1.332311855, 13.4833211, 1.411382563, 54.77131265)),
}


def get_scaled_pdf_func(pdf_func: Callable):
    """ Return a function gives the PDF of a random variable Y = cX,
    where the PDF of X is given by pdf_func. """
    
    def scaled_pdf_func(y: np.ndarray,
                        c: float,
                        *params: float):
        return pdf_func(y / c, *params) / c
    
    return scaled_pdf_func


scaled_kumaraswamy_pdf = get_scaled_pdf_func(kumaraswamy_pdf)
scaled_double_kumaraswamy_pdf = get_scaled_pdf_func(double_kumaraswamy_pdf)


def get_mus_bases(mus: pd.Series | pd.DataFrame):
    """ Get the mutation rates for each base. """
    require_isinstance("mus", mus, (pd.Series, pd.DataFrame))
    bases = mus.index.get_level_values(BASE_NAME)
    invalid_bases = np.setdiff1d(bases, list(REF_PARAMS)).tolist()
    if invalid_bases:
        raise ValueError(
            f"Cannot yet normalize mutation rates for bases: {invalid_bases}"
        )
    return {base: mus.loc[bases == base] for base in REF_PARAMS}


def calc_pdfs(mus_bases: dict[str, pd.Series | pd.DataFrame],
              scale_factor: float):
    """ Calculate the PDFs for each base. """
    pdfs = dict()
    for base, mus_base in mus_bases.items():
        # Calculate pseudoenergies using the method in Cordero et al.
        # Biochemistry 2012, 51, 36, 7037–7039; except use Kumaraswamy
        # distributions with custom empirical parameters, scaled by the
        # scale factor that is optimal for this dataset.
        paired_params, unpaired_params = REF_PARAMS[base]
        paired_pdf = scaled_kumaraswamy_pdf(
            mus_base.values, scale_factor, *paired_params
        )
        unpaired_pdf = scaled_double_kumaraswamy_pdf(
            mus_base.values, scale_factor, *unpaired_params
        )
        pdfs[base] = paired_pdf, unpaired_pdf
    return pdfs


def calc_scale_factor(mus: pd.Series | pd.DataFrame,
                      f_paired: float | int,
                      eps: float = 1.e-6,
                      **kwargs):
    """ Calculate the scale factor for the parameters to fit the given
    mutation rates. """
    require_between("f_paired", f_paired, 0., 1., classes=(float, int))
    require_between("eps", eps, 0., 1., inclusive=False, classes=float)
    mus_bases = get_mus_bases(mus)

    def objective(scale_factor: float):
        loglike = 0.
        pdfs = calc_pdfs(mus_bases, scale_factor)
        for paired_pdf, unpaired_pdf in pdfs.values():
            mixed_pdf = f_paired * paired_pdf + (1. - f_paired) * unpaired_pdf
            loglike += np.sum(np.log(mixed_pdf))
        return -loglike

    from scipy.optimize import minimize_scalar
    result = minimize_scalar(objective, bounds=(eps, 1. / eps), **kwargs)
    if not result.success:
        raise RuntimeError(
            f"Failed to optimize scale factor: {result.message}"
        )
    return float(result.x)


def calc_pseudoenergies(mus: pd.Series | pd.DataFrame,
                        temperature: float | int,
                        f_paired: float | int,
                        **kwargs):
    """ Calculate the pseudoenergy of each base to predict structures,
    in kcal/mol. """
    require_atleast("temperature", temperature, 0., classes=(float, int))
    scale_factor = calc_scale_factor(mus, f_paired, **kwargs)
    # Calculate the pseudoenergies for each type of base via the method
    # in Cordero et al. Biochemistry 2012, 51, 36, 7037–7039;
    # except use Kumaraswamy distributions with custom empirical
    # parameters, scaled by the scale factor that is optimal for this
    # dataset.
    pseudoenergies = mus.copy(deep=True)
    mus_bases = get_mus_bases(mus)
    pdfs = calc_pdfs(mus_bases, scale_factor)
    for base, (paired_pdf, unpaired_pdf) in pdfs.items():
        pseudoenergies.loc[mus_bases[base].index] = (
            (R * temperature) * np.log(unpaired_pdf / paired_pdf)
        )
    return pseudoenergies


def simulate_and_plot_distributions(n_samples: int = 1000,
                                    f_paired: float = 0.5,
                                    base: str = BASEA,
                                    seed: int = 42,
                                    figsize: tuple = (10, 6),
                                    **kwargs):
    """ Simulate mutation rates and plot the original and scaled distributions.
    
    Parameters
    ----------
    n_samples: int
        Number of mutation rates to simulate.
    f_paired: float
        Fraction of paired bases.
    base: str
        Base type to simulate (A, C, G, or U).
    seed: int
        Random seed for reproducibility.
    figsize: tuple
        Figure size (width, height) in inches.
    **kwargs:
        Additional arguments to pass to calc_scale_factor.
        
    Returns
    -------
    tuple
        (scale_factor, figure) - The optimized scale factor and the matplotlib figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy.random import RandomState
    
    # Set random seed for reproducibility
    rng = RandomState(seed)
    
    # Get the parameters for the specified base
    params = {"paired": (2, 48), "unpaired": (10, 10)}
    
    # Simulate mutation rates from the paired and unpaired distributions
    paired_samples = rng.beta(params['paired'][0], params['paired'][1], n_samples)
    unpaired_samples = rng.beta(params['unpaired'][0], params['unpaired'][1], n_samples)
    
    # Mix the samples according to f_paired
    mixed_samples = np.zeros(n_samples)
    paired_mask = rng.random(n_samples) < f_paired
    mixed_samples[paired_mask] = paired_samples[paired_mask]
    mixed_samples[~paired_mask] = unpaired_samples[~paired_mask]
    
    # Create a DataFrame with the simulated mutation rates
    mus = pd.Series(mixed_samples, index=Region("x", "A" * n_samples).range)
    
    # Calculate the optimal scale factor
    scale_factor = calc_scale_factor(mus, f_paired, **kwargs)
    
    # Create a range of x values for plotting the PDFs
    order = np.argsort(mus)
    x = mus.values[order]
    
    # Calculate the PDFs for the original distributions (scale_factor = 1)
    original_pdfs = calc_pdfs({base: mus}, 1.0)
    original_paired_pdf, original_unpaired_pdf = original_pdfs[base]
    original_mixed_pdf = f_paired * original_paired_pdf + (1 - f_paired) * original_unpaired_pdf
    
    # Calculate the PDFs for the scaled distributions
    scaled_pdfs = calc_pdfs({base: mus}, scale_factor)
    scaled_paired_pdf, scaled_unpaired_pdf = scaled_pdfs[base]
    scaled_mixed_pdf = f_paired * scaled_paired_pdf + (1 - f_paired) * scaled_unpaired_pdf
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=figsize)
    
    # Plot the original distributions
    ax1.plot(x, original_paired_pdf[order], 'b-', label='Paired')
    ax1.plot(x, original_unpaired_pdf[order], 'r-', label='Unpaired')
    #ax1.plot(x, original_mixed_pdf, 'g-', label='Mixed')
    ax1.hist(mixed_samples, bins=30, density=True, alpha=0.3, color='gray')
    ax1.set_title(f'Original Distributions (scale_factor = 1.0)')
    ax1.set_xlabel('Mutation Rate')
    ax1.set_ylabel('Density')
    ax1.legend()
    
    # Plot the scaled distributions
    ax2.plot(x, scaled_paired_pdf[order], 'b-', label='Paired')
    ax2.plot(x, scaled_unpaired_pdf[order], 'r-', label='Unpaired')
    #ax2.plot(x, scaled_mixed_pdf, 'g-', label='Mixed')
    ax2.hist(mixed_samples, bins=30, density=True, alpha=0.3, color='gray')
    ax2.set_title(f'Scaled Distributions (scale_factor = {scale_factor:.4f})')
    ax2.set_xlabel('Mutation Rate')
    ax2.set_ylabel('Density')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()
    return scale_factor, fig
