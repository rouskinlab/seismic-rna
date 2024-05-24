import numpy as np
import pandas as pd

from ..core.array import stochastic_round
from ..core.batch import match_reads_segments
from ..core.rel import NOCOV, REL_TYPE
from ..core.seq import Section, index_to_pos
from ..core.stats import calc_dirichlet_params
from ..core.types import fit_uint_type

rng = np.random.default_rng()


def sim_ends(end5: int,
             end3: int,
             end3_mean: float,
             read_mean: float,
             variance: float,
             num_reads: int):
    """ Simulate end coordinates using a Dirichlet distribution.

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


def sim_muts(pmut: pd.DataFrame,
             seg_end5s: np.ndarray,
             seg_end3s: np.ndarray):
    """ Simulate mutation data.

    Parameters
    ----------
    pmut: Section
        Rate of each type of mutation at each position.
    seg_end5s:
        5' end coordinate of each segment.
    seg_end3s:
        3' end coordinate of each segment.
    """
    num_reads, _ = match_reads_segments(seg_end5s, seg_end3s)
    read_nums = np.arange(num_reads, dtype=fit_uint_type(num_reads))
    rels = np.asarray(pmut.columns, dtype=REL_TYPE)
    if NOCOV in rels:
        raise ValueError(f"Mutation types include no coverage: {rels}")
    muts = dict()
    for pos in index_to_pos(pmut.index):
        muts[pos] = dict()
        # Find the reads that cover this position.
        usable_reads = read_nums[np.any(np.logical_and(seg_end5s <= pos,
                                                       pos <= seg_end3s),
                                        axis=1)]
        # Choose a number of reads for each type of relationship.
        num_reads_pos_rels = pd.Series(rng.multinomial(usable_reads.size,
                                                       pmut.loc[pos])[0],
                                       index=rels)
        for rel in rels:
            num_reads_pos_rel = num_reads_pos_rels[rel]
            if num_reads_pos_rel > 0:
                # Randomly select reads with this relationship.
                reads_pos_rel = rng.choice(usable_reads,
                                           num_reads_pos_rel,
                                           replace=False,
                                           shuffle=False)
                muts[pos][rel] = reads_pos_rel
                # Prevent those reads from being chosen for another
                # relationship.
                usable_reads = np.setdiff1d(usable_reads,
                                            reads_pos_rel,
                                            assume_unique=True)
    return muts
