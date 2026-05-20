from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .batch import format_read_name
from .io import ReadNamesBatchIO, RelateBatchIO, RefseqIO
from .report import RelateReport
from ..core import path
from ..core.arg import MUT_COLLISIONS_DROP
from ..core.header import parse_header
from ..core.random import stochastic_round, get_random_integer_generator
from ..core.rel import MATCH, RelPattern
from ..core.seq import BASE_NAME, DNA, index_to_pos
from ..core.tmp import release_to_out
from ..core.unbias import (calc_p_noclose_given_clust,
                           calc_p_noclose_given_ends_auto,
                           calc_p_clust_given_noclose)
from ..core.write import need_write


def make_p_ends_2d(pends: np.ndarray,
                   uniq_end5s: np.ndarray,
                   uniq_end3s: np.ndarray,
                   pmut_index) -> np.ndarray:
    """ Convert a 1-D end-pair probability array to a 2-D (positions ×
    positions) matrix required by calc_p_noclose_given_clust. """
    positions = index_to_pos(pmut_index)
    min_pos = int(positions.min())
    num_pos = len(positions)
    p_ends_2d = np.zeros((num_pos, num_pos))
    np.add.at(p_ends_2d, (uniq_end5s - min_pos, uniq_end3s - min_pos), pends)
    return p_ends_2d


def _calc_p_noclose_given_clust(pmut_given_clust: np.ndarray,
                     p_ends: np.ndarray,
                     min_mut_gap: int) -> float:
    """ Probability that a read from a cluster has no two mutations too close. """
    p_noclose_given_ends = calc_p_noclose_given_ends_auto(pmut_given_clust, min_mut_gap)
    return calc_p_noclose_given_clust(p_ends, p_noclose_given_ends)


def simulate_batch(sample: str,
                   branches: dict[str, str],
                   ref: str,
                   batch: int,
                   write_read_names: bool,
                   formatter: Callable[[int, int], str] = format_read_name,
                   **kwargs):
    """ Simulate a pair of RelateBatchIO and ReadNamesBatchIO. """
    relate_batch = RelateBatchIO.simulate(sample=sample,
                                          branches=branches,
                                          ref=ref,
                                          batch=batch,
                                          **kwargs)
    if write_read_names:
        name_batch = ReadNamesBatchIO.simulate(sample=sample,
                                               branches=branches,
                                               ref=ref,
                                               batch=batch,
                                               num_reads=relate_batch.num_reads,
                                               formatter=formatter)
    else:
        name_batch = None
    return relate_batch, name_batch


def simulate_cluster(first_batch: int,
                     batch_size: int,
                     num_reads: int,
                     seed: int | None,
                     **kwargs):
    """ Simulate all batches for one cluster. """
    seeds = get_random_integer_generator(seed)
    # Determine the numbers of batches and reads per batch.
    num_full_batches, last_batch_size = divmod(int(num_reads), int(batch_size))
    last_batch = first_batch + num_full_batches
    # Simulate every full-size batch.
    for batch in range(first_batch, last_batch):
        yield simulate_batch(batch=batch,
                             num_reads=batch_size,
                             seed=next(seeds),
                             **kwargs)
    # Simulate the last batch, which may have fewer reads.
    if last_batch_size > 0:
        yield simulate_batch(batch=last_batch,
                             num_reads=last_batch_size,
                             seed=next(seeds),
                             **kwargs)


def calc_pmut_pattern(pmut: pd.DataFrame, pattern: RelPattern,
                      normalize: bool = True):
    """ Calculate the rate of a given type of mutation.

    Parameters
    ----------
    normalize: bool
        If True (default), divide by (MATCH + fmut) to get the conditional
        mutation rate among informative positions.  If False, return the
        unconditional probability sum fmut directly — use this when computing
        p_noclose for simulation, where non-informative (AMB) positions are
        treated as non-mutations by reads_noclose_muts.
    """
    header = parse_header(pmut.columns)
    rels = list(map(int, header.get_rel_header().index))
    # Accumulate the frequencies of all selected mutations.
    fmut = pd.DataFrame(0., pmut.index, header.get_clust_header().index)
    for base in DNA.alph():
        positions = fmut.index[fmut.index.get_level_values(BASE_NAME) == base]
        for rel in rels:
            if all(pattern.fits(base, rel)):
                fmut.loc[positions] += pmut.loc[positions, rel]
    if normalize:
        # Conditional rate: mutations / (matches + mutations), i.e. among
        # informative positions only.
        return fmut / (pmut.loc[:, MATCH] + fmut)
    return fmut


def simulate_batches(batch_size: int,
                     pmut: pd.DataFrame,
                     pclust: pd.Series,
                     pends: np.ndarray,
                     uniq_end5s: np.ndarray,
                     uniq_end3s: np.ndarray,
                     num_reads: int,
                     min_mut_gap: int,
                     min_mut_gap_weights: dict[int, float],
                     mut_collisions: str,
                     seed: int | None,
                     **kwargs):
    """
    Simulate batches of reads for all clusters.

    Parameters
    ----------
    batch_size: int
        Number of reads per batch.
    pmut: pd.DataFrame
        Mutation rate DataFrame with columns indexed by (rel, k, clust).
    pclust: pd.Series
        Proportion of reads belonging to each cluster; index is (k, clust).
    num_reads: int
        Total number of reads to simulate across all clusters.
    seed: int | None
        Random seed for reproducibility; None for no fixed seed.
    min_mut_gap_weights: dict[int, float]
        Mapping of min_mut_gap value to weight. When non-empty, reads for
        each cluster are split across the gap values proportionally;
        overrides the single min_mut_gap passed via kwargs. An empty dict
        uses the single min_mut_gap from kwargs.
    mut_collisions: str
        How to handle reads with mutations closer than min_mut_gap.
    **kwargs
        Additional keyword arguments forwarded to `simulate_cluster`.

    Yields
    ------
    tuple[RelateBatchIO, ReadNamesBatchIO | None]
        Pairs of relate and (optionally) read-name batch objects.
    """
    seeds = get_random_integer_generator(seed)
    rng = np.random.default_rng(next(seeds))
    # Simulate the number of reads per cluster, adjusting pclust to
    # pclust_given_noclose when reads with close mutations are dropped.
    first_batch = 0
    if min_mut_gap_weights:
        gap_keys = list(min_mut_gap_weights.keys())
        p_gap = np.array(list(min_mut_gap_weights.values()))
        if mut_collisions == MUT_COLLISIONS_DROP:
            pmut_given_clust = calc_pmut_pattern(pmut, RelPattern.muts(),
                                                 normalize=False).values
            p_ends_2d = make_p_ends_2d(
                pends, uniq_end5s, uniq_end3s, pmut.index
            )
            # The weights are the probability that a read came from the
            # distibution with the minimum gap given that it has no two
            # mutations too close.
            if not np.isclose(p_gap.sum(), 1.0):
                raise ValueError("min_mut_gap_weights must sum to 1, "
                                 f"but got {p_gap.sum()}")
            # Calculate the probability that a read has no mutations too
            # close for each cluster and mutation gap.
            # Shape: (num_gaps, num_clusters)
            p_noclose_given_gap_clust = np.stack([
                _calc_p_noclose_given_clust(pmut_given_clust, p_ends_2d, gap)
                for gap in gap_keys
            ], axis=0)
            # Calculate the unconditional proportion of each gap.
            # P(gap | noclose) = P(gap) / P(noclose) * P(noclose | gap)
            # P(gap) = p(noclose) * P(gap | noclose) / P(noclose | gap)
            p_noclose_given_gap = p_noclose_given_gap_clust @ pclust.values
            p_gap_uncond = p_gap / p_noclose_given_gap
            p_gap_uncond /= p_gap_uncond.sum()  # normalize to sum = 1
            # Calculate the proportion of each cluster after dropping
            # reads with two mutations too close.
            p_noclose_given_clust = (
                p_gap_uncond[np.newaxis, :] @ p_noclose_given_gap_clust
            ).reshape(-1)
            p_clust_given_noclose = calc_p_clust_given_noclose(
                pclust.values, p_noclose_given_clust
            )
            # Calculate the number of reads needed after dropping those
            # with close mutations for each combination of number of
            # clusters and minimum gap between mutations.
            # Shape: (num_gaps, num_clusters)
            num_reads_final = (
                num_reads
                * p_gap[:, np.newaxis]
                * p_clust_given_noclose[np.newaxis, :]
            )
            # Calculate the number of reads that must be simulated to
            # obtain approximately num_reads_final after dropping.
            num_reads_sim = num_reads_final / p_noclose_given_gap_clust
            n_sim = int(stochastic_round(num_reads_sim.sum(), seed=next(seeds)))
            flat = rng.multinomial(n_sim,
                                   (num_reads_sim / num_reads_sim.sum()).ravel())
            shape = num_reads_sim.shape
        else:
            joint = p_gap[:, np.newaxis] * pclust.values[np.newaxis, :]
            flat = rng.multinomial(num_reads, joint.ravel())
            shape = joint.shape
        num_reads_clusters = pd.DataFrame(
            flat.reshape(shape),
            index=gap_keys,
            columns=pclust.index
        )
        # Simulate batches for each cluster and gap.
        for k, clust in pclust.index:
            for gap in min_mut_gap_weights:
                pmut_cluster = pmut.loc[:, (slice(None), k, clust)]
                n_gap = int(num_reads_clusters.at[gap, (k, clust)])
                if n_gap > 0:
                    for batch in simulate_cluster(first_batch,
                                                  batch_size,
                                                  n_gap,
                                                  pmut=pmut_cluster,
                                                  pends=pends,
                                                  seed=next(seeds),
                                                  min_mut_gap=gap,
                                                  mut_collisions=mut_collisions,
                                                  uniq_end5s=uniq_end5s,
                                                  uniq_end3s=uniq_end3s,
                                                  **kwargs):
                        yield batch
                        first_batch += 1
    else:
        if mut_collisions == MUT_COLLISIONS_DROP:
            pmut_given_clust = calc_pmut_pattern(pmut, RelPattern.muts(),
                                                 normalize=False).values
            p_ends_2d = make_p_ends_2d(
                pends, uniq_end5s, uniq_end3s, pmut.index
            )
            p_noclose_given_clust = _calc_p_noclose_given_clust(
                pmut_given_clust, p_ends_2d, min_mut_gap
            )
            p_clust_given_noclose = calc_p_clust_given_noclose(pclust.values,
                                                               p_noclose_given_clust)
            num_reads_final = num_reads * p_clust_given_noclose
            num_reads_sim = num_reads_final / p_noclose_given_clust
            n_sim = int(stochastic_round(num_reads_sim.sum(), seed=next(seeds)))
            num_reads_clusters = pd.Series(
                rng.multinomial(n_sim, num_reads_sim / num_reads_sim.sum()),
                pclust.index
            )
        else:
            num_reads_clusters = pd.Series(rng.multinomial(num_reads, pclust),
                                           pclust.index)
        # Simulate batches for each cluster.
        for k, clust in pclust.index:
            pmut_cluster = pmut.loc[:, (slice(None), k, clust)]
            num_reads_cluster = num_reads_clusters.at[(k, clust)]
            for batch in simulate_cluster(first_batch,
                                          batch_size,
                                          num_reads_cluster,
                                          pmut=pmut_cluster,
                                          pends=pends,
                                          seed=next(seeds),
                                          min_mut_gap=min_mut_gap,
                                          mut_collisions=mut_collisions,
                                          uniq_end5s=uniq_end5s,
                                          uniq_end3s=uniq_end3s,
                                          **kwargs):
                yield batch
                first_batch += 1


def simulate_relate(*,
                    out_dir: Path,
                    tmp_dir: Path,
                    branch: str,
                    sample: str,
                    ref: str,
                    refseq: DNA,
                    write_read_names: bool,
                    brotli_level: int,
                    force: bool,
                    **kwargs):
    """ Simulate an entire relate step. """
    branches = path.add_branch(path.RELATE_STEP, branch, dict())
    report_file = RelateReport.build_path({path.TOP: out_dir,
                                           path.SAMPLE: sample,
                                           path.BRANCHES: branches,
                                           path.REF: ref})
    if need_write(report_file, force):
        began = datetime.now()
        # Write the reference sequence to a file.
        refseq_file = RefseqIO(sample=sample,
                               branches=branches,
                               ref=ref,
                               refseq=refseq)
        _, refseq_checksum = refseq_file.save(tmp_dir,
                                              brotli_level=brotli_level,
                                              force=True)
        # Simulate and write the batches.
        checksums = {RelateBatchIO.btype(): list()}
        if write_read_names:
            checksums[ReadNamesBatchIO.btype()] = list()
        read_count = 0
        for relate_batch, name_batch in simulate_batches(
            sample=sample,
            branches=branches,
            ref=ref,
            write_read_names=write_read_names,
            **kwargs
        ):
            _, relate_checksum = relate_batch.save(tmp_dir,
                                                   brotli_level=brotli_level,
                                                   force=True)
            checksums[RelateBatchIO.btype()].append(relate_checksum)
            if write_read_names:
                assert isinstance(name_batch, ReadNamesBatchIO)
                _, name_checksum = name_batch.save(tmp_dir,
                                                   brotli_level=brotli_level,
                                                   force=True)
                checksums[ReadNamesBatchIO.btype()].append(name_checksum)
            else:
                assert name_batch is None
            read_count += relate_batch.num_reads
        if write_read_names:
            assert ReadNamesBatchIO.btype() in checksums
            assert (len(checksums[ReadNamesBatchIO.btype()])
                    == len(checksums[RelateBatchIO.btype()]))
        else:
            assert ReadNamesBatchIO.btype() not in checksums
        ended = datetime.now()
        # Write the report.
        report = RelateReport(sample=sample,
                              branches=branches,
                              ref=ref,
                              min_mapq=0,
                              phred_enc=0,
                              min_phred=0,
                              insert3=False,
                              ambindel=False,
                              clip_end5=0,
                              clip_end3=0,
                              min_reads=0,
                              n_reads_xam=0,
                              n_reads_rel=read_count,
                              n_batches=len(checksums[RelateBatchIO.btype()]),
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=ended,
                              overhangs=True)
        report_saved = report.save(tmp_dir, force=True)
        release_to_out(out_dir, tmp_dir, report_saved.parent)
    return report_file
