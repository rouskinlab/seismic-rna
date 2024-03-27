from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from .io import QnamesBatchIO, RelateBatchIO
from .report import RelateReport
from .write import get_reads_per_batch, mib_to_bytes
from ..core.batch import END_COORDS, END5_COORD, END3_COORD, get_length
from ..core.dims import triangular
from ..core.io import RefseqIO
from ..core.rel import SUB_A, SUB_C, SUB_G, SUB_T
from ..core.seq import (BASEA,
                        BASEC,
                        BASEG,
                        BASET,
                        DNA,
                        index_to_pos,
                        index_to_seq,
                        seq_pos_to_index)
from ..core.write import need_write

rng = np.random.default_rng()


def _list_naturals(n: int):
    """ List natural numbers. """
    return np.arange(1, n + 1)


def _check_naturals(values: np.ndarray, what: str = "values"):
    """ Raise ValueError if the values are not monotonically increasing
    natural numbers. """
    length = get_length(values, what)
    if not np.array_equal(values, np.arange(1, length + 1)):
        raise ValueError(f"{what} must be numbered 1 to {length}, "
                         f"but got {values}")
    return np.asarray(values, dtype=int)


def sub_options(base: str):
    """ Valid substitutions for a given reference base. """
    if base == BASEA:
        return [SUB_C, SUB_G, SUB_T]
    if base == BASEC:
        return [SUB_A, SUB_G, SUB_T]
    if base == BASEG:
        return [SUB_A, SUB_C, SUB_T]
    if base == BASET:
        return [SUB_A, SUB_C, SUB_G]
    raise ValueError(f"Invalid base: {repr(base)}")


def index_to_refseq_pos(index: pd.MultiIndex):
    # Determine the reference sequence.
    refseq = index_to_seq(index)
    # Validate the sequence positions.
    pos = _check_naturals(index_to_pos(index), "positions")
    return refseq, pos


def choose_clusters(p_clust: pd.Series, n_reads: int):
    """ Choose a cluster for each read. """
    # Validate the cluster proportions.
    if not isinstance(p_clust, pd.Series):
        raise TypeError("p_clust must be Series, "
                        f"but got {type(p_clust).__name__}")
    clusters = _check_naturals(p_clust.index, "clusters")
    # Choose a cluster for each read.
    return rng.choice(clusters, n_reads, p=p_clust.values)


def simulate_p_mut(refseq: DNA, ncls: int):
    """ Simulate mutation rates. """
    return pd.DataFrame(
        rng.random((len(refseq), ncls)),
        seq_pos_to_index(refseq, _list_naturals(len(refseq)), 1),
        _list_naturals(ncls)
    )


def simulate_p_ends(npos: int, min_size: int = 1, max_size: int | None = None):
    """ Simulate proportions of end coordinates. """
    p_ends = pd.Series(
        1. - rng.random(triangular(npos)),
        pd.MultiIndex.from_arrays(
            [a + 1 for a in np.triu_indices(npos)],
            names=END_COORDS
        )
    )
    # Validate the read size limits.
    if not 1 <= min_size <= npos:
        raise ValueError(f"min_size must be ≥ 1 and ≤ {npos}, "
                         f"but got {min_size}")
    if max_size is None:
        max_size = npos
    elif not min_size <= max_size <= npos:
        raise ValueError(f"max_size must be ≥ {min_size} and ≤ {npos}, "
                         f"but got {min_size}")
    # Remove coordinates with improper sizes.
    if min_size > 1 or max_size < npos:
        end5s = p_ends.index.get_level_values(END5_COORD).values
        end3s = p_ends.index.get_level_values(END3_COORD).values
        sizes = end3s - end5s + 1
        p_ends = p_ends.loc[np.logical_and(min_size <= sizes,
                                           sizes <= max_size)]
    return p_ends / p_ends.sum()


def simulate_p_clust(ncls: int):
    p_clust = pd.Series(1. - rng.random(ncls), _list_naturals(ncls))
    return p_clust / p_clust.sum()


def simulate_read_name(batch_num: int, read_num: int):
    return f"Batch:{batch_num};Read:{read_num}"


def simulate_qnames_batch(sample: str,
                          ref: str,
                          batch: int,
                          n_reads: int):
    return QnamesBatchIO(sample=sample,
                         ref=ref,
                         batch=batch,
                         names=[simulate_read_name(batch, read)
                                for read in range(n_reads)])


def simulate_relate_batch(sample: str,
                          ref: str,
                          batch: int,
                          n_reads: int,
                          p_mut: pd.DataFrame,
                          p_ends: pd.Series,
                          cluster_choices: np.ndarray):
    """ Simulate a relate batch. """
    # Validate the mutation rates.
    if not isinstance(p_mut, pd.DataFrame):
        raise TypeError("p_mut must be DataFrame, "
                        f"but got {type(p_mut).__name__}")
    refseq, positions = index_to_refseq_pos(p_mut.index)
    clusters = _check_naturals(p_mut.columns, "clusters")
    # Validate the end coordinates.
    if not isinstance(p_ends, pd.Series):
        raise TypeError("p_ends must be Series, "
                        f"but got {type(p_ends).__name__}")
    if tuple(p_ends.index.names) != END_COORDS:
        raise ValueError(f"p_ends index must have names {END_COORDS}, "
                         f"but got {p_ends.index.names}")
    reads = np.arange(n_reads)
    # Choose 5' and 3' end coordinates for each read.
    coords = rng.choice(p_ends.size, n_reads, p=p_ends.values)
    end5_choices = p_ends.index.get_level_values(END5_COORD).values[coords]
    end3_choices = p_ends.index.get_level_values(END3_COORD).values[coords]
    # Validate the cluster choices.
    if not isinstance(cluster_choices, np.ndarray):
        raise TypeError("cluster_choices must be ndarray, "
                        f"but got {type(cluster_choices).__name__}")
    # Determine which reads come from each cluster.
    in_cluster = {cluster: reads[cluster_choices == cluster]
                  for cluster in clusters}
    # Simulate mutations for each position.
    muts = dict()
    for pos, base in zip(positions, refseq, strict=True):
        pos_sub_options = sub_options(base)
        muts[pos] = {sub: np.array([], dtype=int) for sub in pos_sub_options}
        # Determine for which reads the position is within the ends.
        in_ends = reads[np.logical_and(end5_choices <= pos,
                                       pos <= end3_choices)]
        # Simulate mutations for each cluster.
        for cluster in clusters:
            # Determine which reads can be mutated.
            read_options = np.intersect1d(in_ends,
                                          in_cluster[cluster],
                                          assume_unique=True)
            # Simulate the number of reads with mutations.
            n_muts = rng.binomial(get_length(read_options),
                                  p_mut.at[(pos, base), cluster])
            # Choose reads with mutations at the position.
            read_choices = rng.choice(read_options, n_muts, replace=False)
            # Choose a type of subsitution for each mutated read.
            sub_choices = rng.choice(pos_sub_options, read_choices.size)
            # Record the reads with each type of substitution.
            for sub in pos_sub_options:
                muts[pos][sub] = np.hstack([muts[pos][sub],
                                            read_choices[sub_choices == sub]])
    # Assemble a simulated RelateBatchIO.
    return RelateBatchIO(sample=sample,
                         ref=ref,
                         reflen=len(refseq),
                         batch=batch,
                         end5s=end5_choices,
                         mid5s=end5_choices,
                         mid3s=end3_choices,
                         end3s=end3_choices,
                         muts=muts)


def simulate_batch(sample: str,
                   ref: str,
                   batch: int,
                   n_reads: int,
                   p_mut: pd.DataFrame,
                   p_ends: pd.Series,
                   cluster_choices: np.ndarray):
    qnames_batch = simulate_qnames_batch(sample, ref, batch, n_reads)
    relate_batch = simulate_relate_batch(sample, ref, batch, n_reads,
                                         p_mut, p_ends, cluster_choices)
    return qnames_batch, relate_batch


def generate_batch(out_dir: Path, brotli_level: int, *args, **kwargs):
    qnames_batch, relate_batch = simulate_batch(*args, **kwargs)
    _, relate_check = relate_batch.save(out_dir, brotli_level, force=True)
    _, qnames_check = qnames_batch.save(out_dir, brotli_level, force=True)
    return qnames_check, relate_check


def generate_batches(refseq: DNA,
                     n_reads: int,
                     batch_size: float,
                     cluster_choices: np.ndarray,
                     *args,
                     **kwargs):
    qnames_checks = list()
    relate_checks = list()
    # Determine the numbers of batches and reads per batch.
    n_reads_per_batch = get_reads_per_batch(mib_to_bytes(batch_size),
                                            len(refseq))
    n_batches_full, n_reads_extra = divmod(n_reads, n_reads_per_batch)
    # Generate every full-size batch.
    for batch in range(n_batches_full):
        first = n_reads_per_batch * batch
        last = n_reads_per_batch + first
        qnames_check, relate_check = generate_batch(
            *args,
            batch=batch,
            n_reads=n_reads_per_batch,
            cluster_choices=cluster_choices[first: last],
            **kwargs
        )
        qnames_checks.append(qnames_check)
        relate_checks.append(relate_check)
    # Generate the last batch, which may have fewer reads.
    if n_reads_extra:
        first = n_reads_per_batch * n_batches_full
        last = n_reads
        qnames_check, relate_check = generate_batch(
            *args,
            batch=n_batches_full,
            n_reads=n_reads_extra,
            cluster_choices=cluster_choices[first: last],
            **kwargs
        )
        qnames_checks.append(qnames_check)
        relate_checks.append(relate_check)
    # Collect the checksums into one dictionary.
    checksums = {RelateBatchIO.btype(): relate_checks,
                 QnamesBatchIO.btype(): qnames_checks}
    return checksums


def simulate_relate(out_dir: Path,
                    sample: str,
                    ref: str,
                    refseq: DNA,
                    n_reads: int,
                    brotli_level: int,
                    force: bool,
                    *args,
                    **kwargs):
    """ Simulate an entire relate step. """
    report_file = RelateReport.build_path(top=out_dir,
                                          sample=sample,
                                          ref=ref)
    if need_write(report_file, force):
        began = datetime.now()
        # Write the reference sequence to a file.
        refseq_file = RefseqIO(sample=sample, ref=ref, refseq=refseq)
        _, refseq_checksum = refseq_file.save(out_dir, brotli_level, force=True)
        # Simulate and write the batches.
        checksums = generate_batches(out_dir=out_dir,
                                     sample=sample,
                                     ref=ref,
                                     refseq=refseq,
                                     n_reads=n_reads,
                                     brotli_level=brotli_level,
                                     *args,
                                     **kwargs)
        ns_batches = sorted(set(map(len, checksums.values())))
        if len(ns_batches) != 1:
            raise ValueError("Expected exactly 1 number of batches, "
                             f"but got {ns_batches}")
        n_batches = ns_batches[0]
        ended = datetime.now()
        # Write the report.
        report = RelateReport(sample=sample,
                              ref=ref,
                              min_mapq=-1,
                              min_reads=n_reads,
                              n_reads_rel=n_reads,
                              n_batches=n_batches,
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=ended)
        report.save(out_dir, force=True)
    return report_file
