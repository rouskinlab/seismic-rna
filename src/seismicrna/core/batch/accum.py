from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from .confusion import init_confusion_matrix
from .ends import END_COORDS
from .muts import RegionMutsBatch
from ..header import make_header
from ..logs import logger
from ..rel import RelPattern
from ..seq import DNA, seq_pos_to_index
from ..task import as_list_of_tuples, dispatch
from ..validate import require_isinstance, require_index_equals


def accumulate_counts(
    batch_counts: Iterable[tuple[Any, Any, Any, Any]],
    refseq: DNA,
    pos_nums: np.ndarray,
    patterns: dict[str, RelPattern],
    ks: Iterable[int] | None = None,
    *,
    count_ends: bool = True,
    count_pos: bool = True,
    count_read: bool = True,
    validate: bool = True,
):
    """Accumulate counts from batches into total counts.

    Parameters
    ----------
    batch_counts: Iterable[tuple]
        Iterable of (num_reads, end_counts, count_per_pos, count_per_read)
        tuples, one per batch.
    refseq: DNA
        Reference sequence.
    pos_nums: np.ndarray
        Position numbers to include in the per-position counts.
    patterns: dict[str, RelPattern]
        Mapping from relationship name to relationship pattern.
    ks: Iterable[int] | None
        Numbers of clusters; None for unclustered data.
    count_ends: bool = True
        Whether to accumulate end coordinate counts.
    count_pos: bool = True
        Whether to accumulate per-position counts.
    count_read: bool = True
        Whether to accumulate per-read counts.
    validate: bool = True
        Whether to validate the index and column labels of each batch.

    Returns
    -------
    tuple
        Total (num_reads, end_counts, count_per_pos, count_per_read).
    """
    rels = list(patterns)
    with logger.debug.single_context("accumulating counts of patterns {}", rels):
        header = make_header(rels=rels, ks=ks)
        end_counts_index = pd.MultiIndex.from_arrays(
            [np.array([], dtype=int) for _ in END_COORDS], names=END_COORDS
        )
        # Initialize the total read counts and end coordinate counts.
        if header.get_is_clustered():
            dtype = float
            rel_header = header.get_rel_header()
            clust_index = header.get_clust_header().index
            num_reads = pd.Series(0, index=clust_index)
            end_counts = (
                pd.DataFrame(index=end_counts_index, columns=clust_index, dtype=dtype)
                if count_ends
                else None
            )
        else:
            dtype = int
            rel_header = header
            num_reads = 0
            end_counts = (
                pd.Series(index=end_counts_index, dtype=dtype) if count_ends else None
            )
        zero = dtype(0)
        # Initialize the counts per position.
        count_per_pos = (
            pd.DataFrame(zero, seq_pos_to_index(refseq, pos_nums, 1), header.index)
            if count_pos
            else None
        )
        # Initialize the counts per read.
        count_per_batch_read = list() if count_read else None
        # Accumulate the counts from the batches.
        for (
            num_reads_i,
            end_counts_i,
            count_per_pos_i,
            count_per_read_i,
        ) in batch_counts:
            require_isinstance("num_reads_i", num_reads_i, type(num_reads))
            if validate and isinstance(num_reads, pd.Series):
                require_index_equals(
                    "num_reads_i.index",
                    num_reads_i.index,
                    num_reads.index,
                    "num_reads.index",
                )
            num_reads += num_reads_i
            if end_counts is not None:
                require_isinstance("end_counts_i", end_counts_i, type(end_counts))
                if validate and isinstance(end_counts, pd.DataFrame):
                    require_index_equals(
                        "end_counts_i.columns",
                        end_counts_i.columns,
                        end_counts.columns,
                        "end_counts.columns",
                    )
                end_counts = end_counts.add(end_counts_i, fill_value=zero).astype(
                    dtype, copy=False
                )
            if count_per_pos is not None:
                require_isinstance("count_per_pos_i", count_per_pos_i, pd.DataFrame)
                if validate:
                    require_index_equals(
                        "count_per_pos_i.index",
                        count_per_pos_i.index,
                        count_per_pos.index,
                        "count_per_pos.index",
                    )
                    require_index_equals(
                        "count_per_pos_i.columns",
                        count_per_pos_i.columns,
                        count_per_pos.columns,
                        "count_per_pos.columns",
                    )
                count_per_pos += count_per_pos_i
            if count_per_batch_read is not None:
                require_isinstance("count_per_read_i", count_per_read_i, pd.DataFrame)
                if validate:
                    require_index_equals(
                        "count_per_read_i.columns",
                        count_per_read_i.columns,
                        rel_header.index,
                        "rel_header.index",
                    )
                count_per_batch_read.append(count_per_read_i)
        # Concatenate the per-read counts for the batches.
        if count_per_batch_read:
            count_per_read = pd.concat(count_per_batch_read, axis=0)
        elif count_per_batch_read is not None:
            count_per_read = pd.DataFrame(columns=rel_header.index, dtype=dtype)
        else:
            count_per_read = None
        logger.trace("num_reads={}", num_reads)
        if count_per_pos is not None:
            logger.trace(
                "count_per_pos={}({}, {})",
                type(count_per_pos).__name__,
                count_per_pos.shape,
                count_per_pos.dtypes,
            )
        if count_per_read is not None:
            logger.trace(
                "count_per_read={}({}, {})",
                type(count_per_read).__name__,
                count_per_read.shape,
                count_per_read.dtypes,
            )
    return num_reads, end_counts, count_per_pos, count_per_read


def accumulate_batches(
    get_batch_count_all: Callable[[int], tuple[Any, Any, Any, Any]],
    num_batches: int,
    refseq: DNA,
    pos_nums: np.ndarray,
    patterns: dict[str, RelPattern],
    ks: Iterable[int] | None = None,
    *,
    count_ends: bool = True,
    count_pos: bool = True,
    count_read: bool = True,
    validate: bool = True,
    num_cpus: int = 1,
):
    """Compute and accumulate counts from all batches, optionally in
    parallel.

    Parameters
    ----------
    get_batch_count_all: Callable[[int], tuple]
        Callable that takes a batch number and returns
        (num_reads, end_counts, count_per_pos, count_per_read).
    num_batches: int
        Total number of batches to process.
    refseq: DNA
        Reference sequence.
    pos_nums: np.ndarray
        Position numbers to include in the per-position counts.
    patterns: dict[str, RelPattern]
        Mapping from relationship name to relationship pattern.
    ks: Iterable[int] | None
        Numbers of clusters; None for unclustered data.
    count_ends: bool = True
        Whether to accumulate end coordinate counts.
    count_pos: bool = True
        Whether to accumulate per-position counts.
    count_read: bool = True
        Whether to accumulate per-read counts.
    validate: bool = True
        Whether to validate the index and column labels of each batch.
    num_cpus: int = 1
        Number of CPUs to use for parallel processing.

    Returns
    -------
    tuple
        Total (num_reads, end_counts, count_per_pos, count_per_read).
    """
    with logger.debug.single_context("accumulating counts of {} batches", num_batches):
        # Generate the counts for the batches in parallel.
        batch_counts = dispatch(
            get_batch_count_all,
            num_cpus=num_cpus,
            pass_num_cpus=False,
            as_list=False,
            ordered=False,
            raise_on_error=True,
            args=as_list_of_tuples(range(num_batches)),
            kwargs=dict(
                patterns=patterns,
                ks=ks,
                count_ends=count_ends,
                count_pos=count_pos,
                count_read=count_read,
            ),
        )
        # Accumulate the counts for all batches.
        accum_counts = accumulate_counts(
            batch_counts,
            refseq,
            pos_nums,
            patterns,
            ks,
            count_ends=count_ends,
            count_pos=count_pos,
            count_read=count_read,
            validate=validate,
        )
    return accum_counts


def _calc_batch_confusion_matrix(
    batch_num: int,
    *,
    get_batch: Callable[[int], RegionMutsBatch],
    pattern: RelPattern,
    min_gap: int,
):
    batch = get_batch(batch_num)
    return batch.calc_confusion_matrix(pattern, min_gap)


def accumulate_confusion_matrices(
    get_batch: Callable[[int], RegionMutsBatch],
    num_batches: int,
    pattern: RelPattern,
    pos_index: pd.Index,
    clusters: pd.Index | None,
    min_gap: int = 0,
    num_cpus: int = 1,
):
    """Accumulate confusion matrices from all batches.

    Parameters
    ----------
    get_batch: Callable[[int], RegionMutsBatch]
        Callable that returns the batch for a given batch number.
    num_batches: int
        Total number of batches to process.
    pattern: RelPattern
        Relationship pattern defining which reads count as mutated.
    pos_index: pd.Index
        Index of positions to include in the confusion matrix.
    clusters: pd.Index | None
        Cluster index for the confusion matrix columns; None if not
        clustered.
    min_gap: int = 0
        Minimum gap between positions to include in position pairs.
    num_cpus: int = 1
        Number of CPUs to use for parallel processing.

    Returns
    -------
    tuple
        Total (n, a, b, ab) confusion matrix components.
    """
    with logger.debug.single_context(
        "Accumulating confusion matrices of {} batches", num_batches
    ):
        n, a, b, ab = init_confusion_matrix(pos_index, clusters, min_gap)
        for n_batch, a_batch, b_batch, ab_batch in dispatch(
            _calc_batch_confusion_matrix,
            num_cpus=num_cpus,
            pass_num_cpus=False,
            as_list=False,
            ordered=False,
            raise_on_error=True,
            args=as_list_of_tuples(range(num_batches)),
            kwargs=dict(get_batch=get_batch, pattern=pattern, min_gap=min_gap),
        ):
            n += n_batch
            a += a_batch
            b += b_batch
            ab += ab_batch
    return n, a, b, ab
