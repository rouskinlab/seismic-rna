import os
from logging import getLogger
from shutil import rmtree
from typing import Any

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (opt_struct_file,
                        opt_p_paired,
                        opt_pmut_paired,
                        opt_pmut_unpaired)
from ..core.array import stochastic_round
from ..core.batch import match_reads_segments
from ..core.header import format_clust_name, index_order_clusts
from ..core.rel import (MATCH,
                        NOCOV,
                        DELET,
                        SUB_A,
                        SUB_C,
                        SUB_G,
                        SUB_T,
                        ANY_B,
                        ANY_D,
                        ANY_H,
                        ANY_V,
                        ANY_N,
                        REL_TYPE)
from ..core.rna import RNAProfile, from_ct
from ..core.seq import (BASE_NAME,
                        BASEA,
                        BASEC,
                        BASEG,
                        BASET,
                        BASEN,
                        DNA,
                        Section)
from ..core.stats import calc_beta_params, calc_dirichlet_params
from ..fold.rnastructure import fold

COMMAND = __name__.split(os.path.extsep)[-1]

rng = np.random.default_rng()


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
        num_reads = 100_000
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
