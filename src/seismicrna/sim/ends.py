import os
from pathlib import Path

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (opt_ct_file,
                        opt_center_fmean,
                        opt_center_fvar,
                        opt_length_fmean,
                        opt_length_fvar,
                        opt_force,
                        opt_max_procs)
from ..core.batch import END5_COORD, END3_COORD
from ..core.rna import find_ct_region
from ..core.run import run_func
from ..core.stats import calc_beta_params
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]

rng = np.random.default_rng()

PROPORTION = "Proportion"


def _validate_int(n: int,
                  what: str,
                  minimum: int | None = None,
                  maximum: int | None = None):
    if not isinstance(n, int):
        raise ValueError(f"{what} must be int, but got {type(n).__name__}")
    if minimum is not None:
        _validate_int(minimum, "minimum")
        if n < minimum:
            raise ValueError(f"{what} must be ≥ {minimum}, but got {n}")
    if maximum is not None:
        _validate_int(maximum, "maximum")
        if n > maximum:
            raise ValueError(f"{what} must be ≤ {maximum}, but got {n}")


def _validate_fraction(f: float | int, what: str):
    if not isinstance(f, (float, int)):
        raise TypeError(f"{what} must be float/int, but got {type(f).__name__}")
    if not 0. <= f <= 1.:
        raise ValueError(f"{what} must be ≥ 0 and ≤ 1, but got {f}")


def _calc_p_bins(nbins: int, mean: float, fvar: float):
    """ Calculate the probability of each bin.

    Parameters
    ----------
    nbins: int
        Number of bins to calculate the probability of; must be ≥ 1.
    mean: float
        Mean value; must be ≥ 0 and ≤ nbins - 1.
    fvar: float
        Variance as a fraction of its maximum, from 0 to 1.
    """
    _validate_int(nbins, "nbins", minimum=1)
    if nbins == 1:
        # There is only one bin, which must have probability 1.
        return np.ones(1)
    # Calculate the mean as a fraction of the degrees of freedom.
    dof = nbins - 1
    if not 0. <= mean <= dof:
        raise ValueError(f"mean must be ≥ 0 and ≤ {dof}, but got {mean}")
    fmean = mean / dof
    _validate_fraction(fvar, "fvar")
    if fvar < 1.:
        # Calculate the variance.
        var = fvar * fmean * (1. - fmean)
        if var > 0.:
            # If the variance is > 0, then use a beta distribution
            # (which is continuous) and then integrate over each window
            # between consecutive integers to calculate the probability
            # of each bin.
            from scipy.stats import beta
            pbins = np.diff(beta.cdf(np.linspace(0., 1., nbins + 1),
                                     *calc_beta_params(fmean, var)))
        else:
            # If the variance is 0, all reads must have the same length:
            # namely, mean read length (rounded to the nearest integer).
            # Find the bin in which the mean read length will occur.
            pbins = np.zeros(nbins)
            mean_bin = min(int(fmean * nbins), nbins - 1)
            pbins[mean_bin] = 1.
    else:
        # The variance is maximal, so only the lowest and highest bins
        # can have non-zero probability.
        pbins = np.zeros(nbins)
        pbins[0] = 1. - fmean
        pbins[-1] = fmean
    return pbins


def sim_pends(end5: int,
              end3: int,
              center_fmean: float,
              center_fvar: float,
              length_fmean: float,
              length_fvar: float,
              keep_empty_reads: bool = True):
    """ Simulate segment end coordinate probabilities.

    Parameters
    ----------
    end5: int
        5' end of the region (minimum allowed 5' end coordinate).
    end3: int
        3' end of the region (maximum allowed 5' end coordinate).
    center_fmean: float
        Mean read center, as a fraction of the reference length.
    center_fvar: float
        Variance of the read center, as a fraction of its maximum.
    length_fmean: float
        Mean read length, as a fraction of the available length.
    length_fvar: float
        Variance of the read length, as a fraction of its maximum.
    keep_empty_reads: bool
        Whether to keep reads whose lengths are 0.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        5' and 3' coordinates and their probabilities.
    """
    # Length of the region.
    _validate_int(end5, "end5", minimum=1)
    _validate_int(end3, "end3", minimum=0)
    difference = end3 - end5
    region_length = difference + 1
    if region_length < 0:
        raise ValueError("Length of the region must be ≥ 0, but got "
                         f"{region_length} (end5={end5}, end3={end3})")
    # Mean center (average of 5' and 3' ends) among all reads.
    _validate_fraction(center_fmean, "center_fmean")
    mean_read_center = center_fmean * difference + end5
    # Central position of the region (can be a half-integer).
    region_center = (end5 + end3) / 2
    # Maximum possible mean read length given the mean read center.
    max_mean_read_length = region_length - 2 * abs(mean_read_center
                                                   - region_center)
    # Mean length among all reads.
    _validate_fraction(length_fmean, "length_fmean")
    mean_read_length = length_fmean * max_mean_read_length
    # Calculate the probability of each read length.
    _validate_fraction(length_fvar, "length_fvar")
    p_read_length = _calc_p_bins(region_length + 1,
                                 mean_read_length,
                                 length_fvar)
    # For each read length, calculate the probability of each pair of
    # 5' and 3' ends.
    _validate_fraction(center_fvar, "center_fvar")
    end5s = list()
    end3s = list()
    pends = list()
    for read_length in map(int, np.flatnonzero(p_read_length)):
        if keep_empty_reads or read_length >= 1:
            # Number of positions at which the read can start.
            n_position_options = region_length - read_length + 1
            # Mean 5' end of reads with this length.
            end5_mean = mean_read_center - (read_length - 1) / 2
            # Calculate the probability of each 5'/3' end position.
            p_bins = _calc_p_bins(n_position_options,
                                  end5_mean - end5,
                                  center_fvar) * p_read_length[read_length]
            # Use the positions with non-zero probabilities.
            select = np.flatnonzero(p_bins)
            pends.append(p_bins[select])
            end5s.append(np.arange(end5, end5 + n_position_options)[select])
            end3s.append(end5s[-1] + (read_length - 1))
    if pends:
        end5s = np.concatenate(end5s)
        end3s = np.concatenate(end3s)
        pends = np.concatenate(pends)
    else:
        end5s = np.array([], dtype=int)
        end3s = np.array([], dtype=int)
        pends = np.array([], dtype=float)
    return end5s, end3s, pends


def sim_pends_ct(ct_file: Path, *,
                 center_fmean: float,
                 center_fvar: float,
                 length_fmean: float,
                 length_fvar: float,
                 force: bool):
    pends_file = ct_file.with_suffix(path.PARAM_ENDS_EXT)
    if need_write(pends_file, force):
        region = find_ct_region(ct_file)
        uniq_end5s, uniq_end3s, pends = sim_pends(region.end5,
                                                  region.end3,
                                                  center_fmean,
                                                  center_fvar,
                                                  length_fmean,
                                                  length_fvar,
                                                  keep_empty_reads=False)
        uniq_ends = {END5_COORD: uniq_end5s, END3_COORD: uniq_end3s}
        pends = pd.Series(pends,
                          pd.MultiIndex.from_arrays(list(uniq_ends.values()),
                                                    names=list(uniq_ends)),
                          name=PROPORTION)
        pends.to_csv(pends_file)
    return pends_file


def load_pends(pends_file: Path):
    """ Load end coordinate proportions from a file. """
    data = pd.read_csv(pends_file, index_col=list(range(2)))[PROPORTION]
    index = data.index
    num_segments, odd = divmod(index.nlevels, 2)
    if odd:
        raise ValueError("Number of end coordinates must be even, "
                         f"but got {index.nlevels}")
    uniq_end5s = index.get_level_values(END5_COORD).values
    uniq_end3s = index.get_level_values(END3_COORD).values
    pends = data.values
    return uniq_end5s, uniq_end3s, pends


@run_func(COMMAND)
def run(*,
        ct_file: tuple[str, ...],
        center_fmean: float,
        center_fvar: float,
        length_fmean: float,
        length_fvar: float,
        force: bool,
        max_procs: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(sim_pends_ct,
                    max_procs=max_procs,
                    pass_n_procs=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(center_fmean=center_fmean,
                                center_fvar=center_fvar,
                                length_fmean=length_fmean,
                                length_fvar=length_fvar,
                                force=force))


params = [
    opt_ct_file,
    opt_center_fmean,
    opt_center_fvar,
    opt_length_fmean,
    opt_length_fvar,
    opt_force,
    opt_max_procs
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the proportions of 5' and 3' end coordinates. """
    run(*args, **kwargs)

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
