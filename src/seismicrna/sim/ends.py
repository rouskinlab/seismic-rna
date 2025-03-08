import os
from pathlib import Path
from typing import Iterable

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
                        opt_num_cpus)
from ..core.batch import END5_COORD, END3_COORD
from ..core.rna import find_ct_region
from ..core.run import run_func
from ..core.stats import calc_beta_params
from ..core.task import as_list_of_tuples, dispatch
from ..core.validate import require_atleast, require_between, require_fraction
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]

rng = np.random.default_rng()

PROPORTION = "Proportion"


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
    require_atleast("nbins", nbins, 1, classes=int)
    if nbins == 1:
        # There is only one bin, which must have probability 1.
        return np.ones(1)
    # Calculate the mean as a fraction of the degrees of freedom.
    dof = nbins - 1
    require_between("mean", mean, 0., dof, maximum_name="nbins - 1")
    fmean = mean / dof
    require_fraction("fvar", fvar)
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
    require_atleast("end5", end5, 1, classes=int)
    require_atleast("end3", end3, 0, classes=int)
    difference = end3 - end5
    region_length = difference + 1
    require_atleast("region_length", region_length, 0)
    # Mean center (average of 5' and 3' ends) among all reads.
    require_fraction("center_fmean", center_fmean)
    mean_read_center = center_fmean * difference + end5
    # Central position of the region (can be a half-integer).
    region_center = (end5 + end3) / 2
    # Maximum possible mean read length given the mean read center.
    max_mean_read_length = region_length - 2 * abs(mean_read_center
                                                   - region_center)
    # Mean length among all reads.
    require_fraction("length_fmean", length_fmean)
    mean_read_length = length_fmean * max_mean_read_length
    # Calculate the probability of each read length.
    require_fraction("length_fvar", length_fvar)
    p_read_length = _calc_p_bins(region_length + 1,
                                 mean_read_length,
                                 length_fvar)
    # For each read length, calculate the probability of each pair of
    # 5' and 3' ends.
    require_fraction("center_fvar", center_fvar)
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
        ct_file: Iterable[str | Path],
        center_fmean: float,
        center_fvar: float,
        length_fmean: float,
        length_fvar: float,
        force: bool,
        num_cpus: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(sim_pends_ct,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
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
    opt_num_cpus
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the proportions of 5' and 3' end coordinates. """
    run(*args, **kwargs)
