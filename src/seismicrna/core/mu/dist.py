from typing import Callable

import numpy as np
import pandas as pd

from ..seq.region import BASE_NAME
from ..seq.xna import BASEA, BASEC
from ..stats import kumaraswamy_pdf, double_kumaraswamy_pdf
from ..validate import require_isinstance, require_atleast, require_between


# Gas constant
R = 1.98720425864083e-3  # kcal/(mol*K)


# Reference parameters for each base.
REF_PARAMS = {
    BASEA: (
        (1.30298999, 131.5804808),
        (0.557323482, 1.388569622, 11.16626786, 2.778348856, 618.6712217),
    ),
    BASEC: (
        (0.774939875, 44.15809053),
        (0.374710709, 1.332311855, 13.4833211, 1.411382563, 54.77131265),
    ),
}


def get_scaled_pdf_func(pdf_func: Callable):
    """Return a function gives the PDF of a random variable Y = cX,
    where the PDF of X is given by pdf_func."""

    def scaled_pdf_func(y: np.ndarray, c: float, *params: float):
        return pdf_func(y / c, *params) / c

    return scaled_pdf_func


scaled_kumaraswamy_pdf = get_scaled_pdf_func(kumaraswamy_pdf)
scaled_double_kumaraswamy_pdf = get_scaled_pdf_func(double_kumaraswamy_pdf)


def get_mus_bases(mus: pd.Series | pd.DataFrame):
    """Get the mutation rates for each base."""
    require_isinstance("mus", mus, (pd.Series, pd.DataFrame))
    bases = mus.index.get_level_values(BASE_NAME)
    invalid_bases = np.setdiff1d(bases, list(REF_PARAMS)).tolist()
    if invalid_bases:
        raise ValueError(
            f"Cannot yet normalize mutation rates for bases: {invalid_bases}"
        )
    return {base: mus.loc[bases == base] for base in REF_PARAMS}


def calc_pdfs(
    mus_bases: dict[str, pd.Series | pd.DataFrame],
    scale_factor: float,
    eps: float = 1.0e-6,
):
    """Calculate the PDFs for each base."""
    pdfs = dict()
    for base, mus_base in mus_bases.items():
        # Calculate pseudoenergies using the method in Cordero et al.
        # Biochemistry 2012, 51, 36, 7037–7039; except use Kumaraswamy
        # distributions with custom empirical parameters, scaled by the
        # scale factor that is optimal for this dataset.
        paired_params, unpaired_params = REF_PARAMS[base]
        # The Kumaraswamy distribution has support on the open interval
        # (0, 1), so clip reactivities off the boundary: an exact 0 (or 1)
        # makes the density diverge (x**(a - 1) with a < 1), which turns
        # the log-likelihood into NaN and defeats the scale-factor fit.
        mus = np.clip(mus_base.values, eps, 1.0 - eps)
        # A scale factor smaller than a reactivity pushes it past the
        # distribution's support (mus / scale_factor > 1), giving a NaN
        # density; the optimizer simply rejects such scale factors, so
        # ignore the expected out-of-support warnings here.
        with np.errstate(invalid="ignore", divide="ignore"):
            paired_pdf = scaled_kumaraswamy_pdf(mus, scale_factor, *paired_params)
            unpaired_pdf = scaled_double_kumaraswamy_pdf(
                mus, scale_factor, *unpaired_params
            )
        pdfs[base] = paired_pdf, unpaired_pdf
    return pdfs


def calc_scale_factor(
    mus: pd.Series | pd.DataFrame, f_paired: float | int, eps: float = 1.0e-6, **kwargs
):
    """Calculate the scale factor for the parameters to fit the given
    mutation rates."""
    require_between("f_paired", f_paired, 0.0, 1.0, classes=(float, int))
    require_between("eps", eps, 0.0, 1.0, inclusive=False, classes=float)
    mus_bases = get_mus_bases(mus)

    def objective(scale_factor: float):
        loglike = 0.0
        pdfs = calc_pdfs(mus_bases, scale_factor, eps)
        for paired_pdf, unpaired_pdf in pdfs.values():
            mixed_pdf = f_paired * paired_pdf + (1.0 - f_paired) * unpaired_pdf
            loglike += np.sum(np.log(mixed_pdf))
        return -loglike

    from scipy.optimize import minimize_scalar

    result = minimize_scalar(objective, bounds=(eps, 1.0 / eps), **kwargs)
    if not result.success:
        raise RuntimeError(f"Failed to optimize scale factor: {result.message}")
    return float(result.x)


def calc_pseudoenergies(
    mus: pd.Series | pd.DataFrame,
    temperature: float | int,
    f_paired: float | int,
    eps: float = 1.0e-6,
    **kwargs,
):
    """Calculate the pseudoenergy of each base to predict structures,
    in kcal/mol."""
    require_atleast("temperature", temperature, 0.0, classes=(float, int))
    scale_factor = calc_scale_factor(mus, f_paired, eps=eps, **kwargs)
    # Calculate the pseudoenergies for each type of base via the method
    # in Cordero et al. Biochemistry 2012, 51, 36, 7037–7039;
    # except use Kumaraswamy distributions with custom empirical
    # parameters, scaled by the scale factor that is optimal for this
    # dataset.
    pseudoenergies = mus.copy(deep=True)
    mus_bases = get_mus_bases(mus)
    pdfs = calc_pdfs(mus_bases, scale_factor, eps)
    for base, (paired_pdf, unpaired_pdf) in pdfs.items():
        pseudoenergies.loc[mus_bases[base].index] = (R * temperature) * np.log(
            unpaired_pdf / paired_pdf
        )
    return pseudoenergies
