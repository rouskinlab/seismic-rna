import os
from pathlib import Path

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (docdef,
                        opt_ct_file,
                        opt_end3_fmean,
                        opt_read_fmean,
                        opt_ends_var,
                        opt_force,
                        opt_parallel,
                        opt_max_procs)
from ..core.array import stochastic_round
from ..core.batch import END5_COORD, END3_COORD, match_reads_segments
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rna import find_ct_section
from ..core.stats import calc_dirichlet_params
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]

rng = np.random.default_rng()

PROPORTION = "Proportion"


def _sim_ends(end5: int,
              end3: int,
              end3_mean: float,
              read_mean: float,
              variance: float,
              num_reads: int):
    """ Simulate segment end coordinates using a Dirichlet distribution.

    Parameters
    ----------
    end5: int
        5' end of the section (minimum allowed 5' end coordinate).
    end3: int
        3' end of the section (maximum allowed 5' end coordinate).
    end3_mean: float
        Mean 3' end coordinate; must be ≥ `end5` and ≤ `end3`.
    read_mean: float
        Mean read length.
    variance: float
        Variance as a fraction of its supremum; must be ≥ 0 and < 1.
    num_reads: int
        Number of reads to simulate.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        5' and 3' end coordinates of each read.
    """
    if not 1 <= end5 <= end3_mean <= end3:
        raise ValueError("Must have 1 ≤ end5 ≤ end3_mean ≤ end3, but got "
                         f"end5={end5}, end3_mean={end3_mean}, end3={end3}")
    gap3_mean = end3 - end3_mean
    if not 0 <= read_mean <= end3_mean - (end5 - 1):
        raise ValueError("Must have 1 ≤ read_mean ≤ (end3_mean - (end5 - 1)), "
                         f"but got read_mean={read_mean}")
    gap5_mean = end3_mean - (end5 - 1) - read_mean
    interval = end3 - (end5 - 1)
    means = np.array([gap5_mean, read_mean, gap3_mean]) / interval
    variances = variance * (means * (1. - means))
    is_nonzero = variances != 0.
    if is_nonzero.any():
        alpha = calc_dirichlet_params(means[is_nonzero],
                                      variances[is_nonzero])
        sample = rng.dirichlet(alpha, size=num_reads).T * interval
        if is_nonzero.all():
            gap5s, _, gap3s = sample
        elif not is_nonzero[0]:
            diffs, gap3s = sample
            gap5s = np.zeros(num_reads)
        elif not is_nonzero[1]:
            gap5s, gap3s = sample
        elif not is_nonzero[2]:
            gap5s, diffs = sample
            gap3s = np.full(num_reads, gap3_mean)
        else:
            raise
        end5s = stochastic_round(gap5s) + end5
        end3s = end3 - stochastic_round(gap3s)
    else:
        end5s = np.full(num_reads, end3_mean - (read_mean - 1))
        end3s = np.full(num_reads, end3_mean)
    return end5s[:, np.newaxis], end3s[:, np.newaxis]


def sim_pends(end5: int,
              end3: int,
              end3_mean: float,
              read_mean: float,
              variance: float,
              num_reads: int | None = None):
    """ Simulate segment end coordinate probabilities.

    Parameters
    ----------
    end5: int
        5' end of the section (minimum allowed 5' end coordinate).
    end3: int
        3' end of the section (maximum allowed 5' end coordinate).
    end3_mean: float
        Mean 3' end coordinate; must be ≥ `end5` and ≤ `end3`.
    read_mean: float
        Mean read length.
    variance: float
        Variance as a fraction of its supremum; must be ≥ 0 and < 1.
    num_reads: int | None = None
        Number of reads to use for simulation; if omitted, will

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        5' and 3' coordinates and their probabilities.
    """
    if num_reads is None:
        num_reads = 1_000_000
    elif num_reads < 1:
        raise ValueError(f"num_reads must be ≥ 1, but got {num_reads}")
    end5s, end3s = _sim_ends(end5,
                             end3,
                             end3_mean,
                             read_mean,
                             variance,
                             num_reads)
    _, num_segs = match_reads_segments(end5s, end3s)
    ends = np.hstack([end5s, end3s])
    uniq_ends, counts = np.unique(ends, return_counts=True, axis=0)
    uniq_end5s = uniq_ends[:, :num_segs]
    uniq_end3s = uniq_ends[:, num_segs:]
    pends = counts / num_reads
    return uniq_end5s, uniq_end3s, pends


def _format_coord_name(end: str, segment: int):
    return f"{end} {segment + 1}"


def sim_pends_ct(ct_file: Path, *,
                 end3_fmean: float,
                 read_fmean: float,
                 ends_var: float,
                 force: bool):
    pends_file = ct_file.with_suffix(path.PARAM_ENDS_EXT)
    if need_write(pends_file, force):
        section = find_ct_section(ct_file)
        end3_mean = end3_fmean * (section.length - 1) + section.end5
        read_mean = read_fmean * section.length
        uniq_end5s, uniq_end3s, pends = sim_pends(section.end5,
                                                  section.end3,
                                                  end3_mean,
                                                  read_mean,
                                                  ends_var)
        _, num_segs = match_reads_segments(uniq_end5s, uniq_end3s)
        uniq_ends = {_format_coord_name(end, i): coords[:, i]
                     for end, coords in [(END5_COORD, uniq_end5s),
                                         (END3_COORD, uniq_end3s)]
                     for i in range(num_segs)}
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
    uniq_end5s = np.stack(
        [index.get_level_values(_format_coord_name(END5_COORD, i)).values
         for i in range(num_segments)],
        axis=1
    )
    uniq_end3s = np.stack(
        [index.get_level_values(_format_coord_name(END3_COORD, i)).values
         for i in range(num_segments)],
        axis=1
    )
    pends = data.values
    return uniq_end5s, uniq_end3s, pends


@docdef.auto()
def run(ct_file: tuple[str, ...],
        end3_fmean: float,
        read_fmean: float,
        ends_var: float,
        force: bool,
        parallel: bool,
        max_procs: int):
    """ Simulate the rate of each kind of mutation at each position. """
    return dispatch(sim_pends_ct,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=as_list_of_tuples(map(Path, ct_file)),
                    kwargs=dict(end3_fmean=end3_fmean,
                                read_fmean=read_fmean,
                                ends_var=ends_var,
                                force=force))


params = [
    opt_ct_file,
    opt_end3_fmean,
    opt_read_fmean,
    opt_ends_var,
    opt_force,
    opt_parallel,
    opt_max_procs
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate the proportions of 5' and 3' end coordinates. """
    run(*args, **kwargs)
