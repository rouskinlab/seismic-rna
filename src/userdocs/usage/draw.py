""" Illustrate the steps of the workflow. """

import unittest as ut
from typing import Iterable

import numpy as np
from matplotlib import pyplot as plt

from seismicrna.core.array import find_dims
from seismicrna.core.random import stochastic_round
from seismicrna.core.rel import NOCOV, MATCH, DELET, SUB_A, SUB_C, SUB_G, SUB_T
from seismicrna.core.rna import RNAStructure, parse_db_structure
from seismicrna.core.seq import DNA, Section, BASEA, BASEC, BASEG, BASET
from seismicrna.graph.color import tetra

rng = np.random.default_rng(0)

READS = "reads"
POSITIONS = "positions"
CLUSTERS = "clusters"
GAP = "-"


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
    dims = find_dims([(READS,),
                      (READS,),
                      (READS, POSITIONS),
                      (POSITIONS, CLUSTERS),
                      (CLUSTERS,)],
                     [end5s, end3s, muts, mu, pi],
                     ["end5s", "end3s", "muts", "mu", "pi"])
    n_reads = dims[READS]
    n_clust = dims[CLUSTERS]
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
                pa: float = 1.,
                pb: float = 99.,
                ua: float = 2.,
                ub: float = 6.):
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


def simulate_ends(n_reads: int, n_pos: int, p5: float = 0.25, p3: float = 0.25):
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
    dims = find_dims([(POSITIONS, CLUSTERS), (CLUSTERS,), (READS,), (READS,)],
                     [mu, pi, end5s, end3s],
                     ["mu", "pi", "end5s", "end3s"])
    n_reads = dims[READS]
    n_pos = dims[POSITIONS]
    n_clust = dims[CLUSTERS]
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
    dims = find_dims([(READS, POSITIONS), (READS,), (READS,)],
                     [muts, end5s, end3s],
                     ["muts", "end5s", "end3s"])
    n_reads = dims[READS]
    n_pos = dims[POSITIONS]
    rels = np.full((n_reads, n_pos), NOCOV)
    for i in range(n_reads):
        end5 = end5s[i]
        end3 = end3s[i]
        rels[i, end5: end3 + 1] = MATCH
        for j in np.flatnonzero(muts[i]):
            base = seq[j]
            if base == BASEA:
                mut_options = [SUB_C, SUB_G, SUB_T]
            elif base == BASEC:
                mut_options = [SUB_A, SUB_G, SUB_T]
            elif base == BASEG:
                mut_options = [SUB_A, SUB_C, SUB_T]
            elif base == BASET:
                mut_options = [SUB_A, SUB_C, SUB_G]
            else:
                raise ValueError(base)
            if end5 < j < end3:
                # Deletions are allowed except at the ends of the read.
                mut_options.append(DELET)
            rels[i, j] = rng.choice(mut_options)
    return rels


def calc_is_paired(seq: DNA, ss_dbs: Iterable[str]):
    """ Determine whether each base is paired or unpaired in the given
    dot-bracket structures. """
    section = Section("", seq)
    return np.stack(
        [RNAStructure(section=section,
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
                    bases.append(GAP)
            else:
                raise ValueError(rel)
        reads.append("".join(bases))
    return reads


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
                color=tetra.get(base))


def draw_reads(filename: str,
               seq: DNA,
               rels: np.ndarray,
               end5s: np.ndarray | None = None):
    """ Draw reads on an axis and save to a file. """
    # Set parameters.
    vspacing_mm = 6.  # vertical spacing (millimeters)
    hspacing_mm = 5.  # horizontal spacing (millimeters)
    fontsize_mm = 4.  # font size (millimeters)
    # Generate reads.
    is_alignment = end5s is not None
    reads = generate_reads(seq, rels, is_alignment)
    n_reads = len(reads)
    # Graph reads.
    fig, ax = plt.subplots()
    fontsize_pt = mm_to_point(fontsize_mm)
    if is_alignment:
        x_max = len(seq)
        y_min = -n_reads
        y_line = y_min + (1. + fontsize_mm / vspacing_mm) / 2.
        ax.plot([0, x_max],
                [y_line, y_line],
                color="#3f3f3f",
                linewidth=mm_to_point(0.5))
        draw_read(ax, n_reads, str(seq), 0, fontsize_pt)
    else:
        x_max = max(map(len, reads))
        y_min = 1 - n_reads
        end5s = [0] * n_reads
    for i, (read, end5) in enumerate(zip(reads, end5s, strict=True)):
        draw_read(ax, i, read, end5, fontsize_pt)
    # Adjust dimensions and spacing.
    ax.set_xlim((0, x_max))
    ax.set_ylim((y_min, 1))
    ax.set_aspect(vspacing_mm / hspacing_mm)
    ax.axis("off")
    fig.set_size_inches(x_max * mm_to_inch(hspacing_mm),
                        n_reads * mm_to_inch(vspacing_mm))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    # Save the figure.
    plt.savefig(filename)
    plt.close()


def main():
    # Parameters.
    n_reads = 40
    n_clust = 2
    seq = DNA("GACCGAGTCACCTACGGA")
    ss_dbs = ["..(((((....)).))).",
              "(((...))).((...))."]
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
    print(pi)
    end5s, end3s = simulate_ends(n_reads, n_pos, p5=0.1, p3=0.1)
    muts = simulate_muts(mu, pi, end5s, end3s)
    rels = simulate_rels(seq, muts, end5s, end3s)
    # Draw data for each step.
    draw_reads("reads.pdf", seq, rels)
    draw_reads("align.pdf", seq, rels, end5s)


if __name__ == "__main__":
    main()
    ut.main()
