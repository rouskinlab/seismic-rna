""" Illustrate the steps of the workflow. """

import os
import unittest as ut
from functools import reduce
from operator import or_
from typing import Iterable

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from seismicrna.core.array import find_dims
from seismicrna.core.random import stochastic_round
from seismicrna.core.rel import NOCOV, MATCH, DELET, SUB_A, SUB_C, SUB_G, SUB_T
from seismicrna.core.rna import RNAStructure, parse_db_structure
from seismicrna.core.seq import DNA, Region, BASEA, BASEC, BASEG, BASET
from seismicrna.graph.color import get_cmap, RelColorMap, SeqColorMap
from seismicrna.core.table import (COVER_REL,
                                   MATCH_REL,
                                   MUTAT_REL,
                                   DELET_REL,
                                   SUB_A_REL,
                                   SUB_C_REL,
                                   SUB_G_REL,
                                   SUB_T_REL)

rng = np.random.default_rng(42)

ALIGN_CLIP = 1
READ = "Read"
POSITION = "Position"
CLUSTER = "Cluster"
DELET_RATIO = 0.25  # number of deletions relative to each substitution
ALIGN_GAP = "-"
VGRID_MM = 4.  # vertical grid spacing (millimeters)
HGRID_MM = 3.  # horizontal grid spacing (millimeters)
GRID_GAP_MM = 0.5  # gap between grid elements (millimeters)
COVER_TICK_INC = 10
MUTAT_TICK_INC = 0.1
TICK_COLOR = "#E0E0E0"
TICK_WIDTH = 0.5  # millimeters
MASK_COLOR = "#D0D0D0"
FILE_FORMAT = "svg"
DRAW_DIR = "draw"
USE_MUTS = [SUB_A, SUB_C, SUB_G, SUB_T]
MUTAT = reduce(or_, USE_MUTS)
UNINF = MUTAT | DELET
UNINF_COLOR = "#D0D0D0"


def mm_to_inch(mm: float):
    """ Convert millimeters to inches. """
    return mm / 25.4


def inch_to_point(inch: float):
    """ Convert inches to points. """
    return inch * 72.


def mm_to_point(mm: float):
    """ Convert millimeters to points. """
    return inch_to_point(mm_to_inch(mm))


def calc_p_clust_given_read(end5s: np.ndarray,
                            end3s: np.ndarray,
                            muts: np.ndarray,
                            mu: np.ndarray,
                            pi: np.ndarray):
    """ Calculate the probability that each read came from each cluster.

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end of each read.
    end3s: np.ndarray
        1D (reads) array of the 3' end of each read.
    muts: np.ndarray
        2D (reads x positions) boolean array of whether each position is
        mutated in each read.
    mu: np.ndarray
        2D (positions x clusters) array of the probability that each
        position is mutated in each cluster.
    pi: np.ndarray
        1D (cluster) array of the proportion of each cluster.

    Returns
    -------
    np.ndarray
        2D (reads x clusters) array of the probability that each read
        came from each cluster.
    """
    dims = find_dims([(READ,),
                      (READ,),
                      (READ, POSITION),
                      (POSITION, CLUSTER),
                      (CLUSTER,)],
                     [end5s, end3s, muts, mu, pi],
                     ["end5s", "end3s", "muts", "mu", "pi"])
    n_reads = dims[READ]
    n_clust = dims[CLUSTER]
    # Calculate the joint probability that each read would have its
    # mutations and come from each cluster.
    log_mu = np.log(mu)
    log_nu = np.log(1. - mu)
    log_pi = np.log(pi)
    logp_clust_and_read = np.empty((n_reads, n_clust))
    for i in range(n_reads):
        end5 = end5s[i]
        end3_plus_1 = end3s[i] + 1
        for k in range(n_clust):
            logp_clust_and_read[i, k] = log_pi[k] + np.sum(
                np.where(muts[i, end5: end3_plus_1],
                         log_mu[end5: end3_plus_1, k],
                         log_nu[end5: end3_plus_1, k])
            )
    # Calculate the conditional probability that each read would come
    # from each cluster.
    logp_read = np.logaddexp.reduce(logp_clust_and_read, axis=1)
    logp_clust_given_read = logp_clust_and_read - logp_read[:, np.newaxis]
    return np.exp(logp_clust_given_read)


class TestCalcPClustGivenRead(ut.TestCase):

    def test_1_cluster(self):
        pi = np.ones(1)
        for n_pos in range(1, 5):
            mu = rng.random((n_pos, 1))
            for n_reads in range(5):
                end5s = rng.integers(n_pos, size=n_reads)
                end3s = rng.integers(end5s, n_pos, size=n_reads)
                muts = rng.integers(2, size=(n_reads, n_pos)).astype(bool)
                result = calc_p_clust_given_read(end5s, end3s, muts, mu, pi)
                expect = np.ones((n_reads, 1))
                self.assertTrue(np.array_equal(result, expect))

    def test_2_clusters(self):
        end5s = np.array(
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2]
        )
        end3s = np.array(
            [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        )
        muts = np.array([
            # 0 - 0
            [False, False, False],
            [True, False, False],
            # 0 - 1
            [False, False, False],
            [False, True, False],
            [True, False, False],
            [True, True, False],
            # 0 - 2
            [False, False, False],
            [False, False, True],
            [False, True, False],
            [False, True, True],
            [True, False, False],
            [True, False, True],
            [True, True, False],
            [True, True, True],
            # 1 - 2
            [False, False, False],
            [False, False, True],
            [False, True, False],
            [False, True, True],
            # 2 - 2
            [False, False, False],
            [False, False, True]
        ])
        mu = np.array([[0.45, 0.20],
                       [0.25, 0.15],
                       [0.10, 0.30]])
        pi = np.array([0.6, 0.4])
        result = calc_p_clust_given_read(end5s, end3s, muts, mu, pi)
        expect_joint = np.array([
            [0.33, 0.32],
            [0.27, 0.08],
            [0.2475, 0.272],
            [0.0825, 0.048],
            [0.2025, 0.068],
            [0.0675, 0.012],
            [0.22275, 0.19040],
            [0.02475, 0.08160],
            [0.07425, 0.03360],
            [0.00825, 0.01440],
            [0.18225, 0.04760],
            [0.02025, 0.02040],
            [0.06075, 0.00840],
            [0.00675, 0.00360],
            [0.4050, 0.238],
            [0.0450, 0.102],
            [0.1350, 0.042],
            [0.0150, 0.018],
            [0.54, 0.28],
            [0.06, 0.12]
        ])
        expect = expect_joint / expect_joint.sum(axis=1)[:, np.newaxis]
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


def simulate_mu(is_paired: np.ndarray,
                pa: float = 2.,
                pb: float = 64.,
                ua: float = 8.,
                ub: float = 16.):
    """ Simulate mutation rates. """
    n_pos, n_clust = is_paired.shape
    mu = np.zeros((n_pos, n_clust))
    for k in range(n_clust):
        n_paired = np.count_nonzero(is_paired[:, k])
        n_unpaired = n_pos - n_paired
        mu[is_paired[:, k], k] = rng.beta(pa, pb, size=n_paired)
        mu[~is_paired[:, k], k] = rng.beta(ua, ub, size=n_unpaired)
    return mu


def simulate_pi(n_clust: int, alpha: float = 1.):
    """ Simulate cluster proportions. """
    return rng.dirichlet(np.repeat(alpha, n_clust))


def simulate_ends(n_reads: int,
                  n_pos: int,
                  p5: float = 1 / 3,
                  p3: float = 1 / 3):
    """ Simulate 5' and 3' end coordinates. """
    partition = rng.multinomial(n_pos, [p5, 1. - (p5 + p3), p3], n_reads)
    end5s = partition[:, 0]
    end3s = (n_pos - 1) - partition[:, -1]
    return end5s, end3s


def simulate_muts(mu: np.ndarray,
                  pi: np.ndarray,
                  end5s: np.ndarray,
                  end3s: np.ndarray):
    """ Simulate mutations as a boolean array. """
    dims = find_dims([(POSITION, CLUSTER), (CLUSTER,), (READ,), (READ,)],
                     [mu, pi, end5s, end3s],
                     ["mu", "pi", "end5s", "end3s"])
    n_reads = dims[READ]
    n_pos = dims[POSITION]
    n_clust = dims[CLUSTER]
    # Assign each read to a cluster.
    n_reads_per_clust = stochastic_round(pi * n_reads, preserve_sum=True)
    ks = rng.permutation(np.repeat(np.arange(n_clust), n_reads_per_clust))
    # Choose which positions are mutated in each read.
    muts = list()
    for end5, end3, k in zip(end5s, end3s, ks, strict=True):
        muts.append(rng.random(n_pos) < mu[:, k])
        # Remove any mutations outside the 5'/3' ends.
        muts[-1][:end5] = False
        muts[-1][end3 + 1:] = False
    return np.vstack(muts)


def simulate_rels(seq: DNA,
                  muts: np.ndarray,
                  end5s: np.ndarray,
                  end3s: np.ndarray):
    """ Simulate the relationship for each mutation. """
    dims = find_dims([(READ, POSITION), (READ,), (READ,)],
                     [muts, end5s, end3s],
                     ["muts", "end5s", "end3s"])
    n_reads = dims[READ]
    n_pos = dims[POSITION]
    if DELET_RATIO > 1.:
        n_del = round(DELET_RATIO)
        n_sub = 1
    elif DELET_RATIO > 0.:
        n_del = 1
        n_sub = round(1 / DELET_RATIO)
    else:
        n_del = 0
        n_sub = 1
    rels = np.full((n_reads, n_pos), NOCOV)
    for i in range(n_reads):
        rels[i, end5s[i]: end3s[i] + 1] = MATCH
        for j in np.flatnonzero(muts[i]):
            base = seq[j]
            if base == BASEA:
                sub_options = [SUB_C, SUB_G, SUB_T]
            elif base == BASEC:
                sub_options = [SUB_A, SUB_G, SUB_T]
            elif base == BASEG:
                sub_options = [SUB_A, SUB_C, SUB_T]
            elif base == BASET:
                sub_options = [SUB_A, SUB_C, SUB_G]
            else:
                raise ValueError(base)
            mut_options = sub_options * n_sub + [DELET] * n_del
            rels[i, j] = rng.choice(mut_options)
    return rels


def remove_edge_muts(rels: np.ndarray,
                     end5s: np.ndarray,
                     end3s: np.ndarray,
                     min_match: int = ALIGN_CLIP):
    """ Remove mutations that are too close to the 5'/3' ends. """
    n_reads, n_pos = rels.shape
    for i in range(n_reads):
        while (end5s[i] <= end3s[i] and not np.all(
                rels[i, end5s[i]: end5s[i] + min_match] == MATCH
        )):
            end5s[i] += 1
        while (end3s[i] >= end5s[i] and not np.all(
                rels[i, max(end3s[i] - min_match + 1, 0): end3s[i] + 1] == MATCH
        )):
            end3s[i] -= 1
        rels[i, :end5s[i]] = NOCOV
        rels[i, end3s[i] + 1:] = NOCOV


def clip_rels(rels: np.ndarray,
              end5s: np.ndarray,
              end3s: np.ndarray,
              n_clip: int = ALIGN_CLIP):
    """ Clip 5'/3' bases from the reads. """
    n_reads, n_pos = rels.shape
    for i in range(n_reads):
        end5s[i] = min(end5s[i] + n_clip, end3s[i] + 1)
        end3s[i] = max(end3s[i] - n_clip, end5s[i] - 1)
        rels[i, :end5s[i]] = NOCOV
        rels[i, end3s[i] + 1:] = NOCOV


def calc_is_paired(seq: DNA, ss_dbs: Iterable[str]):
    """ Determine whether each base is paired or unpaired in the given
    dot-bracket structures. """
    region = Region("", seq)
    return np.stack(
        [RNAStructure(region=region,
                      title=str(k),
                      pairs=parse_db_structure(ss)).is_paired.values
         for k, ss in enumerate(ss_dbs)],
        axis=1
    )


def generate_reads(seq: DNA, rels: np.ndarray, gaps: bool = False):
    """ Generate read sequences. """
    n_reads, n_pos = rels.shape
    reads = list()
    for i in range(n_reads):
        bases = list()
        for j in range(n_pos):
            rel = rels[i, j]
            if rel == NOCOV:
                continue
            elif rel == MATCH:
                bases.append(seq[j])
            elif rel == SUB_A:
                bases.append(BASEA)
            elif rel == SUB_C:
                bases.append(BASEC)
            elif rel == SUB_G:
                bases.append(BASEG)
            elif rel == SUB_T:
                bases.append(BASET)
            elif rel == DELET:
                if gaps:
                    bases.append(ALIGN_GAP)
            else:
                raise ValueError(rel)
        reads.append("".join(bases))
    return reads


def calc_coverage(rels: np.ndarray):
    """ Calculate the coverage at each position. """
    return np.count_nonzero(rels != NOCOV, axis=0)


def calc_mus_avg(mu: np.ndarray, pi: np.ndarray):
    """ Calculate the mutation rates of the population average. """
    return mu @ pi


def no_margins():
    """ Eliminate margins from the graph. """
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


def draw_read(ax: plt.Axes,
              i: int,
              read: str,
              end5: int,
              fontsize_pt: float,
              fontfamily="monospace"):
    """ Draw one read on the axis. """
    for j, base in enumerate(read, start=end5):
        ax.text(j, -i, base,
                fontfamily=fontfamily,
                fontsize=fontsize_pt,
                color=get_cmap(SeqColorMap).get(base))


def draw_reads(filename: str,
               seq: DNA,
               rels: np.ndarray,
               end5s: np.ndarray | None = None):
    """ Draw reads on an axis and save to a file. """
    # Set parameters.
    fontsize_mm = 4.  # font size (millimeters)
    # Generate reads.
    reads = generate_reads(seq, rels, end5s is not None)
    n_reads = len(reads)
    if end5s is None:
        end5s = np.repeat(0, n_reads)
    end3s = end5s + np.array(list(map(len, reads))) - 1
    # Graph reads.
    fig, ax = plt.subplots()
    fontsize_pt = mm_to_point(fontsize_mm)
    x_max = end3s.max() + 1
    y_min = 1 - n_reads
    for i, (read, end5) in enumerate(zip(reads, end5s, strict=True)):
        draw_read(ax, i, read, end5, fontsize_pt)
    # Adjust dimensions and spacing.
    ax.set_xlim((0, x_max))
    ax.set_ylim((y_min, 1))
    ax.set_aspect(VGRID_MM / HGRID_MM)
    ax.axis("off")
    fig.set_size_inches(x_max * mm_to_inch(HGRID_MM),
                        n_reads * mm_to_inch(VGRID_MM))
    no_margins()
    # Save the figure.
    plt.savefig(os.path.join(DRAW_DIR, filename))
    plt.close()


def draw_rels(filename: str,
              rels: np.ndarray,
              alpha: np.ndarray | None = None,
              mask_pos: np.ndarray | None = None,
              mask_reads: np.ndarray | None = None):
    """ Draw reads on an axis and save to a file. """
    n_reads, n_pos = rels.shape
    # Graph reads.
    fig, ax = plt.subplots()
    x_max = n_pos
    y_min = -n_reads
    reads = np.arange(n_reads)
    if mask_reads is not None:
        reads = reads[~mask_reads]
    positions = np.arange(n_pos)
    if mask_pos is not None:
        positions = positions[~mask_pos]
    hgap = GRID_GAP_MM / HGRID_MM
    vgap = GRID_GAP_MM / VGRID_MM
    width = 1. - hgap
    height = 1. - vgap
    for i in reads:
        for j in positions:
            rel = rels[i, j]
            if rel != NOCOV:
                if rel == UNINF:
                    rel_color = UNINF_COLOR
                else:
                    if rel == MATCH:
                        rel_name = MATCH_REL
                    elif rel == MUTAT:
                        rel_name = MUTAT_REL
                    elif rel == DELET:
                        rel_name = DELET_REL
                    elif rel == SUB_A:
                        rel_name = SUB_A_REL
                    elif rel == SUB_C:
                        rel_name = SUB_C_REL
                    elif rel == SUB_G:
                        rel_name = SUB_G_REL
                    elif rel == SUB_T:
                        rel_name = SUB_T_REL
                    else:
                        raise ValueError(rel)
                    rel_color = get_cmap(RelColorMap)[rel_name]
                xy = j + hgap / 2, vgap / 2 - (i + 1)
                ax.add_patch(
                    Rectangle(xy=xy,
                              width=width,
                              height=height,
                              facecolor=rel_color,
                              alpha=(alpha[i] if alpha is not None else None))
                )
    # Strike out positions and reads.
    if mask_reads is not None:
        for i in -np.flatnonzero(mask_reads) - 0.5:
            ax.plot([0, x_max],
                    [i, i],
                    color=MASK_COLOR,
                    linewidth=mm_to_point(VGRID_MM - GRID_GAP_MM))
    if mask_pos is not None:
        for j in np.flatnonzero(mask_pos) + 0.5:
            ax.plot([j, j],
                    [y_min, 0],
                    color=MASK_COLOR,
                    linewidth=mm_to_point(HGRID_MM - GRID_GAP_MM))
    # Adjust dimensions and spacing.
    ax.set_xlim((0, x_max))
    ax.set_ylim((y_min, 0))
    ax.set_aspect(VGRID_MM / HGRID_MM)
    ax.axis("off")
    fig.set_size_inches(x_max * mm_to_inch(HGRID_MM),
                        n_reads * mm_to_inch(VGRID_MM))
    no_margins()
    # Save the figure.
    plt.savefig(os.path.join(DRAW_DIR, filename))
    plt.close()


def calc_mu_y_max(mu: np.ndarray):
    """ Calculate the maximum y tick. """
    return np.ceil(mu.max(initial=0) / MUTAT_TICK_INC) * MUTAT_TICK_INC


def graph_profile(filename: str,
                  data: np.ndarray,
                  color: str,
                  y_tick_inc: float,
                  y_max: float | None = None,
                  start: int = 1):
    """ Graph a profile. """
    n_pos, = data.shape
    end = n_pos + start
    if y_max is None:
        y_max = data.max(initial=0)
    x_min = start - 0.5
    x_max = end - 0.5
    fig, ax = plt.subplots()
    # Graph the profile.
    ax.bar(np.arange(start, end), data, facecolor=color)
    # Add tick marks.
    for y in np.arange(y_tick_inc, y_max + y_tick_inc, y_tick_inc):
        ax.plot([x_min, x_max],
                [y, y],
                color=TICK_COLOR,
                linewidth=mm_to_point(TICK_WIDTH),
                zorder=0.)
    # Adjust dimensions and spacing.
    ax.set_xlim((x_min, x_max))
    ax.set_ylim((0, y_max))
    side_length = (x_max - x_min) * mm_to_inch(HGRID_MM)
    fig.set_size_inches(side_length, side_length)
    no_margins()
    # Save the figure.
    plt.savefig(os.path.join(DRAW_DIR, filename))
    plt.close()


def graph_cov(*args, **kwargs):
    """ Graph a coverage profile. """
    graph_profile(*args, **kwargs,
                  color=get_cmap(RelColorMap)[COVER_REL],
                  y_tick_inc=COVER_TICK_INC)


def graph_mus(filename: str,
              data: np.ndarray,
              *args,
              y_max: float | None = None,
              **kwargs):
    """ Graph a mutation rate profile. """
    graph_profile(filename, data, *args, **kwargs,
                  y_max=(y_max if y_max is not None else calc_mu_y_max(data)),
                  color=get_cmap(RelColorMap)[MUTAT_REL],
                  y_tick_inc=MUTAT_TICK_INC)


def main():
    try:
        os.mkdir(DRAW_DIR)
    except FileExistsError:
        pass
    # Parameters.
    n_reads = 30
    n_clust = 2
    seq = DNA("GACCGAGTCACCTACGGA")
    ss_dbs = ["..(((((....)).))).",
              "(((...))).((...))."]
    region_end5 = 1
    region_end3 = 14
    min_ncov_read = 6
    min_ninfo_pos = 10
    min_mut_gap = 1
    # Derive properties from the sequence and structure.
    n_pos = len(seq)
    is_paired = calc_is_paired(seq, ss_dbs)
    # Give Gs and Us low mutation rates.
    for j, base in enumerate(seq):
        if base != BASEA and base != BASEC:
            is_paired[j] = True
    # Simulate data.
    mu = simulate_mu(is_paired)
    pi = simulate_pi(n_clust, alpha=6.)
    end5s, end3s = simulate_ends(n_reads, n_pos, p5=0.1, p3=0.1)
    rels = simulate_rels(seq, simulate_muts(mu, pi, end5s, end3s), end5s, end3s)
    # Illustrate align.
    draw_reads(f"reads.{FILE_FORMAT}", seq, rels)
    remove_edge_muts(rels, end5s, end3s)
    draw_reads(f"align.{FILE_FORMAT}", seq, rels, end5s)
    # Illustrate relate.
    draw_rels(f"relate-0.{FILE_FORMAT}", rels)
    clip_rels(rels, end5s, end3s)
    draw_rels(f"relate-1.{FILE_FORMAT}", rels)
    graph_cov(f"relate-1-cov.{FILE_FORMAT}", calc_coverage(rels))
    graph_mus(f"relate-1-mus.{FILE_FORMAT}", calc_mus_avg(mu, pi))
    # Illustrate mask: define region.
    end5s = np.maximum(end5s, region_end5)
    end3s = np.minimum(end3s, region_end3)
    rels = rels[:, region_end5: region_end3 + 1]
    mu = mu[region_end5: region_end3 + 1]
    seq = seq[region_end5: region_end3 + 1]
    draw_rels(f"mask-0.{FILE_FORMAT}", rels)
    # Illustrate mask: exclude positions.
    mask_pos = np.array([base not in "AC" for base in seq])
    draw_rels(f"mask-1.{FILE_FORMAT}", rels,
              mask_pos=mask_pos)
    # Illustrate mask: define mutations.
    muts = np.isin(rels, USE_MUTS)
    matches = rels == MATCH
    rels = np.full_like(rels, UNINF)
    rels[np.nonzero(matches)] = MATCH
    rels[np.nonzero(muts)] = MUTAT
    rels[:, mask_pos] = NOCOV
    for i, (end5, end3) in enumerate(zip(end5s, end3s, strict=True)):
        rels[i, :end5 - region_end5] = NOCOV
        rels[i, end3 - region_end5 + 1:] = NOCOV
    draw_rels(f"mask-2.{FILE_FORMAT}", rels,
              mask_pos=mask_pos)
    # Illustrate mask: mask reads.
    mask_reads = np.count_nonzero(rels != NOCOV, axis=1) < min_ncov_read
    cumsum = np.cumsum(rels == MUTAT, axis=1)
    mask_reads |= np.max(
        cumsum[:, min_mut_gap + 1:] - cumsum[:, :-(min_mut_gap + 1)],
        axis=1
    ) > 1
    rels[mask_reads] = NOCOV
    draw_rels(f"mask-3.{FILE_FORMAT}", rels,
              mask_pos=mask_pos, mask_reads=mask_reads)
    graph_cov(f"mask-3-cov.{FILE_FORMAT}",
              calc_coverage(rels),
              start=region_end5)
    # Illustrate mask: mask positions.
    mask_pos |= np.count_nonzero(rels != NOCOV, axis=0) < min_ninfo_pos
    rels[:, mask_pos] = NOCOV
    mu[mask_pos] = 0.
    draw_rels(f"mask-4.{FILE_FORMAT}", rels,
              mask_pos=mask_pos, mask_reads=mask_reads)
    graph_cov(f"mask-4-cov.{FILE_FORMAT}",
              calc_coverage(rels),
              start=region_end5)
    graph_mus(f"mask-4-mus.{FILE_FORMAT}",
              calc_mus_avg(mu, pi),
              start=region_end5)
    # Illustrate cluster:
    rels[np.nonzero(rels == UNINF)] = MATCH
    p_clust = calc_p_clust_given_read(
        end5s, end3s, np.equal(rels, MUTAT), mu, pi
    )
    for k in range(n_clust):
        draw_rels(f"cluster-0-{k + 1}.{FILE_FORMAT}", rels,
                  alpha=p_clust[:, k], mask_pos=mask_pos, mask_reads=mask_reads)
        graph_mus(f"cluster-0-{k + 1}-mus.{FILE_FORMAT}", mu[:, k],
                  start=region_end5, y_max=calc_mu_y_max(mu))
        print(f"Mutation rates for cluster {k + 1}:")
        print(np.round(mu[:, k], 3))


if __name__ == "__main__":
    main()
    ut.main(verbosity=2)
