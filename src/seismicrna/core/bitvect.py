"""
Core -- Bit Vector Module
========================================================================
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from itertools import chain
from logging import getLogger
from typing import Callable, Generator, Iterable

import numpy as np
import pandas as pd

from .mu import calc_f_obs_df
from .sect import Section

logger = getLogger(__name__)

BIT_VECTOR_NAME = "Bit Vector"


class FeedClosedBitAccumError(Exception):
    """ Feed another batch to a closed BitAccum. """


class InconsistentSectionError(Exception):
    """ Sections do not match. """


class DuplicateIndexError(Exception):
    """ A value in an index is duplicated. """


class UniqMutBits(object):
    """ Collection of unique bit vectors indicating only mutations. """

    def __init__(self, muts: np.ndarray):
        """ Find the unique bit vectors. """
        if muts.ndim != 2:
            raise ValueError(f"muts needs 2 dimensions, but got {muts.ndim}")
        # Count each unique bit vector.
        uniq, self._inv, self._count = np.unique(muts, axis=0,
                                                 return_inverse=True,
                                                 return_counts=True)
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the entire sparse matrix (uniq) because
        # mutations are relatively rare.
        self._uniq_idxs = tuple(map(np.flatnonzero, uniq.T))

    @property
    def indexes(self):
        """ For each position, the indexes of all unique bit vectors
        that have a mutation at that position. """
        return self._uniq_idxs

    @property
    def counts(self) -> np.ndarray:
        """ Number of times each unique bit vector occurs in the given
        set of possibly redundant bit vectors. """
        return self._count

    @property
    def inverse(self) -> np.ndarray:
        """ Indexes to map the unique bit vectors back to the given set
        of possibly redundant bit vectors. """
        return self._inv

    @property
    def n_uniq(self):
        """ Number of unique bit vectors. """
        return self.counts.size

    @property
    def n_pos(self):
        """ Number of positions in each bit vector. """
        return len(self.indexes)

    def get_full(self):
        """ Full boolean matrix of the unique bit vectors. """
        # Initialize an all-False matrix with one row for each unique
        # bit vector and one column for each position.
        full = np.zeros((self.n_uniq, self.n_pos), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in enumerate(self.indexes):
            full[indexes, j] = True
        return full

    def get_uniq_names(self):
        """ Return the unique bit vectors as byte strings. """
        # Get the full boolean matrix of the unique bit vectors and cast
        # the data from boolean to unsigned 8-bit integer type.
        chars = self.get_full().astype(np.uint8, copy=False)
        if chars.size > 0:
            # Add ord('0') to transform every 0 into b'0' and every 1
            # into # b'1', and convert each row (bit vector) into a
            # bytes object of b'0' and b'1' characters.
            names = np.apply_along_axis(np.ndarray.tobytes, 1, chars + ord('0'))
        else:
            # If there are no unique bit vectors, then apply_along_axis
            # will fail, so set names to an empty list.
            names = list()
        return pd.Index(names, name=BIT_VECTOR_NAME)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (self.n_pos == other.n_pos
                and np.array_equal(self.counts, other.counts)
                and np.array_equal(self.inverse, other.inverse)
                and all(np.array_equal(self_idxs, other_idxs)
                        for self_idxs, other_idxs in zip(self.indexes,
                                                         other.indexes,
                                                         strict=True)))


class BitVectorBase(ABC):
    """ Base class for bit vectors. """

    def __init__(self, section: Section):
        self._section = section

    @property
    def section(self):
        """ Section from which the bit vectors came. """
        return self._section

    @property
    @abstractmethod
    def nmasked(self) -> Counter[str]:
        """ Number of reads removed for each mask. """
        return Counter()

    @property
    @abstractmethod
    def reads(self):
        """ Names of the reads (after masking). """
        return pd.Index()

    @property
    def nreads(self):
        """ Number of reads (after masking). """
        return self.reads.size

    @property
    def nreads_given(self):
        """ Number of reads given (before masking). """
        return self.nreads + self.nmasked.total()

    @property
    @abstractmethod
    def n_info_per_pos(self) -> pd.Series:
        """ Number of informative bits at each position. """
        return pd.Series(0, index=self.section.unmasked)

    @property
    @abstractmethod
    def n_info_per_read(self) -> pd.Series:
        """ Number of informative bits in each read. """
        return pd.Series([], dtype=int)

    @property
    @abstractmethod
    def n_affi_per_pos(self) -> pd.Series:
        """ Number of affirmative bits at each position. """
        return pd.Series(0, index=self.section.unmasked)

    @property
    @abstractmethod
    def n_affi_per_read(self) -> pd.Series:
        """ Number of affirmative bits in each read. """
        return pd.Series([], dtype=int)

    @property
    def f_info_per_pos(self):
        """ Fraction of informative bits at each position. """
        return self.n_info_per_pos / self.nreads

    @property
    def f_info_per_read(self):
        """ Fraction of informative bits in each read. """
        return self.n_info_per_read / self.section.size

    @property
    def f_affi_per_pos(self):
        """ Fraction of affirmative bits at each position. """
        return self.n_affi_per_pos / self.n_info_per_pos

    @property
    def f_affi_per_read(self):
        """ Fraction of affirmative bits for each read. """
        return self.n_affi_per_read / self.n_info_per_read

    @abstractmethod
    def _drop_duplicate_reads(self, *args):
        """ Drop any read with a duplicate name. """


class BitMatrix(BitVectorBase, ABC):
    """ Bit vectors represented with two explicit boolean matrices of
    informative and affirmative bits. """

    @property
    @abstractmethod
    def info(self):
        """ Boolean DataFrame indicating informative bits. """
        return pd.DataFrame()

    @property
    @abstractmethod
    def affi(self):
        """ Boolean DataFrame indicating affirmative bits. """
        return pd.DataFrame()

    @property
    def reads(self):
        return self.info.index

    @property
    def n_info_per_pos(self):
        """ Number of informative bits for each position. """
        return pd.Series(np.count_nonzero(self.info, axis=0),
                         index=self.section.unmasked)

    @property
    def n_info_per_read(self):
        """ Number of informative bits for each read. """
        return pd.Series(np.count_nonzero(self.info, axis=1),
                         index=self.reads)

    @property
    def n_affi_per_pos(self):
        """ Number of affirmative bits for each position. """
        return pd.Series(np.count_nonzero(self.affi, axis=0),
                         index=self.section.unmasked)

    @property
    def n_affi_per_read(self):
        """ Number of affirmative bits for each read. """
        return pd.Series(np.count_nonzero(self.affi, axis=1),
                         index=self.reads)

    def _drop(self, drop: pd.Index):
        """ Drop the reads in `drop`; return the number dropped. """
        if drop.size > 0:
            logger.debug(f"Dropping {drop} from {self}")
            self.info.drop(index=drop, inplace=True)
            self.affi.drop(index=drop, inplace=True)
        return drop.size

    def _drop_duplicate_reads(self):
        dups = self.reads[self.reads.duplicated(keep="first")]
        if self._drop(dups):
            logger.warning(f"{self} got duplicate read names: {dups.to_list()}")
        return dups


class BitBatch(BitMatrix):
    """ One batch of bit vectors. """

    def __init__(self, section: Section,
                 info: pd.DataFrame, affi: pd.DataFrame,
                 mask: dict[str, Callable[[BitBatch], pd.Index]] | None = None):
        super().__init__(section)
        # Initialize the reads.
        self._info = info
        self._affi = affi
        # Validate the reads.
        self._check_indexes()
        self._drop_duplicate_reads()
        # Mask the reads.
        self._masked = Counter(self._mask(mask or dict()))

    def _check_indexes(self):
        """ Verify that the read names and positions in info and muts
        are consistent with each other and with the section. """
        logger.debug(f"Checking indexes of {self}")
        if not self.info.index.equals(self.affi.index):
            raise InconsistentSectionError(
                f"Read names for informative\n{self.info.index}\nand "
                f"affirmative\n{self.affi.index}\nbits did not match")
        if not self.info.columns.equals(self.affi.columns):
            raise InconsistentSectionError(
                f"Positions for informative\n{self.info.columns}\nand "
                f"affirmative\n{self.affi.columns}\nbits did not match")
        if not self.info.columns.equals(self.section.unmasked):
            raise InconsistentSectionError(
                f"Positions for bits\n{self.info.columns}\nand section\n"
                f"{self.section.unmasked}\ndid not match")

    @property
    def info(self):
        """ DataFrame of informative bits. """
        return self._info

    @property
    def affi(self):
        """ DataFrame of affirmative bits. """
        return self._affi

    def drop_reads(self, drop: pd.Index):
        """ Drop the reads in `drop`; return the number dropped. """
        return self._drop(drop)

    @property
    def nmasked(self):
        return self._masked

    def _mask(self, masks: dict[str, Callable[[BitBatch], np.ndarray]]):
        """ Drop reads selected by any of the masks, which should be
        boolean NumPy arrays. """
        logger.debug(f"Masking {self} with {masks}")
        return {name: self._drop(self.reads[mask(self)])
                for name, mask in masks.items()}


class ClusterBitBatch(BitBatch):
    """ One batch of bit vectors with cluster membership weights. """

    def __init__(self, section: Section,
                 info: pd.DataFrame, affi: pd.DataFrame,
                 resps: pd.DataFrame):
        super().__init__(section, info, affi)
        self._resps = resps

    @property
    def clusters(self):
        """ Index of the clusters. """
        return self._resps.columns

    @property
    def n_info_per_pos(self):
        # Informative bits per position (row) and cluster (column).
        return self.info.T @ self._resps

    @property
    def n_affi_per_pos(self):
        # Affirmative bits per position (row) and cluster (column).
        return self.affi.T @ self._resps

    @property
    def n_info_per_read(self):
        # Informative bits per read (row) and cluster (column).
        return self._resps.mul(super().n_info_per_read, axis=0)

    @property
    def n_affi_per_read(self):
        # Affirmative bits per read (row) and cluster (column).
        return self._resps.mul(super().n_affi_per_read, axis=0)


class BitAccum(BitVectorBase, ABC):
    """ Accumulates batches of bit vectors. """

    def __init__(self, section: Section, batches: Iterable[BitBatch]):
        super().__init__(section)
        # Initialize the number of batches given.
        self._nbatches = 0
        # Initialize the numbers of reads masked.
        self._masked: Counter[str] = Counter()
        for batch in batches:
            logger.debug(f"Add batch {self.nbatches + 1} ({batch}) to {self}")
            # Confirm that the sections from which the batch derives
            # matches the section from which the accumulator derives.
            if batch.section != self.section:
                raise InconsistentSectionError(
                    f"Sections of the batch ({batch.section}) and accumulator "
                    f"({self.section}) do not match")
            # Drop any reads whose names match reads that already exist.
            self._drop_duplicate_reads(batch)
            # Update the counts of the numbers of reads masked.
            self._masked += batch.nmasked
            logger.debug(f"Current masked counts for {self}: {self.nmasked}")
            # Add the counts of informative and affirmative bits.
            self._count_info_affi(batch)
            # Increment the number of batches given.
            self._nbatches += 1

    @abstractmethod
    def _empty_accum(self):
        """ Empty BitAccum. """

    @abstractmethod
    def _count_info_affi(self, batch: BitBatch):
        """ Count the informative and affirmative bits and add them to
        the accumulator. """

    @property
    def nmasked(self):
        return self._masked

    @property
    def nbatches(self):
        """ Number of batches given to the accumulator. """
        return self._nbatches

    @property
    @abstractmethod
    def _accum_info(self):
        """ Return the accumulated informative bits. """

    def _get_accum_read_names(self):
        """ Names of the reads accumulated so far. """
        if isinstance(self._accum_info, pd.DataFrame):
            # The read names are the index of the DataFrame.
            return self._accum_info.index.to_list()
        if self.nbatches > 0:
            # At least one batch of reads has been given.
            return list(chain.from_iterable(batch.index for batch
                                            in self._accum_info))
        # No reads have been given.
        return list()

    def _drop_duplicate_reads(self, batch: BitBatch):
        dups = batch.reads.intersection(self._get_accum_read_names())
        if batch.drop_reads(dups):
            logger.warning(f"{self} got read(s) in {batch} with name(s) seen "
                           f"in a previous batch: {dups.to_list()}")
        return dups


class BitMonolith(BitAccum, BitMatrix):
    """ Accumulates batches of bit vectors into one monolithic unit. """

    def __init__(self, section: Section, batches: Iterable[BitBatch]):
        # Initialize lists of each batch's total and affirmative bits.
        self._info: list[pd.DataFrame] | pd.DataFrame = list()
        self._affi: list[pd.DataFrame] | pd.DataFrame = list()
        super().__init__(section, batches)
        if self.nbatches > 0:
            # Merge the informative and affirmative bits to DataFrames.
            self._info = pd.concat(self._info, axis=0)
            self._affi = pd.concat(self._affi, axis=0)
        else:
            # No batches were given.
            self._info = self._empty_accum()
            self._affi = self._empty_accum()

    def _empty_accum(self):
        """ Empty DataFrame whose columns are the section indexes. """
        return pd.DataFrame(index=[], columns=self.section.unmasked,
                            dtype=int)

    def _count_info_affi(self, batch: BitBatch):
        # Add the informative and affirmative bits from this batch to
        # the totals among all batches.
        self._info.append(batch.info)
        self._affi.append(batch.affi)

    @property
    def _accum_info(self):
        return self._info

    @property
    def info(self) -> pd.DataFrame:
        return self._info

    @property
    def affi(self) -> pd.DataFrame:
        return self._affi

    def get_uniq_muts(self):
        return UniqMutBits(self.affi.values)


class BitCounter(BitAccum):
    """ Accumulates batches of bit vectors into counts of informative
    and affirmative bits per position and per read. """

    def __init__(self, section: Section, batches: Iterable[BitBatch]):
        # Initialize the counts of informative and affirmative bits.
        self._info_per_pos = self._init_per_pos(section)
        self._affi_per_pos = self._init_per_pos(section)
        self._info_per_read: list[pd.Series] = list()
        self._affi_per_read: list[pd.Series] = list()
        super().__init__(section, batches)

    def _init_per_pos(self, section: Section):
        """ Initialize a count of 0 for each position. """
        return pd.Series(0, index=section.unmasked)

    def _empty_accum(self):
        return pd.Series([], dtype=int)

    def _count_info_affi(self, batch: BitBatch):
        # Add the counts for this batch to the totals.
        self._info_per_pos += batch.n_info_per_pos
        self._affi_per_pos += batch.n_affi_per_pos
        self._info_per_read.append(batch.n_info_per_read)
        self._affi_per_read.append(batch.n_affi_per_read)
        logger.debug(f"Added batch {len(self._info_per_read)} to {self}")
        logger.debug(f"Counts:\n{self._info_per_pos}\n{self._affi_per_pos}")

    @property
    def _accum_info(self):
        return self._info_per_read

    @property
    def n_info_per_pos(self):
        return self._info_per_pos

    @cached_property
    def n_info_per_read(self):
        if self.nbatches == 0:
            # No batches were given.
            return self._empty_accum()
        return pd.concat(self._info_per_read, axis=0)

    @property
    def n_affi_per_pos(self):
        return self._affi_per_pos

    @cached_property
    def n_affi_per_read(self) -> pd.Series:
        if self.nbatches == 0:
            # No batches were given.
            return self._empty_accum()
        return pd.concat(self._affi_per_read, axis=0)

    @property
    def reads(self):
        """ Read names. """
        return self.n_info_per_read.index

    @property
    def nreads(self):
        return sum(read_batch.size for read_batch in self._info_per_read)

    @property
    def read_batches(self):
        """ Read names in each batch. """
        for read_batch in self._info_per_read:
            yield read_batch.index


class ClustBitCounter(BitCounter):
    """ Accumulates batches of bit vectors and cluster memberships into
    counts of informative and affirmative bits per position and per read
    for each cluster. """

    def __init__(self, section: Section, clusters: pd.Index,
                 batches: Iterable[ClusterBitBatch]):
        self._clusters = clusters
        super().__init__(section, batches)

    @property
    def clusters(self):
        """ Cluster orders and numbers. """
        return self._clusters

    def _init_per_pos(self, section: Section):
        """ Initialize a count of 0 for each position and cluster. """
        return pd.DataFrame(0, index=section.unmasked, columns=self.clusters)

    def _empty_accum(self):
        return pd.DataFrame(columns=self.clusters, dtype=int)


def iter_all_bit_vectors(mu: pd.Series, section: Section, min_mut_gap: int):
    """ Yield every bit vector in decreasing order of likelihood. """
    # Determine the number of positions.
    if mu.size == 0:
        # If there are no mutation rates, then return an empty iterator.
        return iter(list())
    # Compute the log probability that a random bit vector would have no
    # pair of mutations too close together.
    logf_obs = float(np.log(calc_f_obs_df(mu.to_frame(),
                                          section,
                                          min_mut_gap).squeeze(axis=0)))
    with np.errstate(divide="ignore"):
        # For each position, compute the log probability that the base
        # is mutated (mu) or is not mutated (nu). Ignore warnings about
        # taking the log of 0 because any positions for which mu or nu
        # equals 0 will be skipped later by checking for np.isfinite.
        log_mu = np.log(mu.values)
        log_nu = np.log(1. - mu.values)
    # Determine the base bit for each position: 1 if mu > 0.5, else 0.
    base_bits = np.greater(log_mu, log_nu)
    if min_mut_gap > 0 and np.any(base_bits):
        raise ValueError("Mutation rates > 0.5 are not currently supported "
                         "when min_mut_gap > 0")
    # Compute the log probability that a random bit vector will be the
    # base bit vector.
    base_logp = float(np.where(base_bits, log_mu, log_nu).sum()) - logf_obs
    # Compute the absolute difference in log probability between the two
    # states of each bit (i.e. 0 and 1).
    diff_logp = np.abs(log_mu - log_nu)

    class BitVectorSeries(object):

        def __init__(self, bit: int, logp: float,
                     sub: Callable[[], Iterable[tuple[np.ndarray, float]]]):
            if not 0 <= bit < mu.size:
                raise ValueError(
                    f"bit must be in [0, {mu.size}), but got {bit}")
            self._bit = bit
            self._onehot = np.zeros(mu.size, dtype=bool)
            self._onehot[bit] = 1
            self._logp = logp
            self._sub = sub
            self._bvecs: list[np.ndarray] = list()
            self._logps: list[float] = list()

        def _cache(self, bvec: np.ndarray, logp: float):
            """ Cache a bit vector and its log probability. """
            self._bvecs.append(bvec)
            self._logps.append(logp)

        def iter(self) -> Generator[tuple[np.ndarray, float], None, None]:
            """ Iterate over the values. """
            # Check if any bit vectors have already been computed.
            if self._bvecs:
                # If so, then yield them.
                yield from zip(self._bvecs, self._logps, strict=True)
            else:
                # Otherwise, make iterators for the subordinate series
                # and the current series.
                sub = iter(self._sub())
                if min_mut_gap > 0:
                    # Exclude any bit vectors with two mutations closer
                    # than min_mut_gap positions.
                    new = iter((np.logical_xor(bvec, self._onehot),
                                logp - self._logp)
                               for bvec, logp in self._sub()
                               if np.all(np.abs(np.flatnonzero(bvec)
                                                - self._bit)
                                         > min_mut_gap))
                else:
                    new = iter((np.logical_xor(bvec, self._onehot),
                                logp - self._logp)
                               for bvec, logp in self._sub())
                # Get the first bit vector from each iterator.
                sub_bvec, sub_logp = next(sub, (None, None))
                new_bvec, new_logp = next(new, (None, None))
                # Yield the most likely of the new and subordinate bit
                # vectors until either iterator is exhausted.
                while new_bvec is not None and sub_bvec is not None:
                    if new_logp > sub_logp:
                        # Yield the new bit vector and iterate.
                        self._cache(new_bvec, new_logp)
                        yield new_bvec, new_logp
                        new_bvec, new_logp = next(new, (None, None))
                    else:
                        # Yield the subordinate bit vector and iterate.
                        self._cache(sub_bvec, sub_logp)
                        yield sub_bvec, sub_logp
                        sub_bvec, sub_logp = next(sub, (None, None))
                # Yield the remaining bit vectors from the iterator that
                # has not been exhausted.
                while sub_bvec is not None:
                    # Yield the subordinate bit vector and iterate.
                    self._cache(sub_bvec, sub_logp)
                    yield sub_bvec, sub_logp
                    sub_bvec, sub_logp = next(sub, (None, None))
                while new_bvec is not None:
                    # Yield the new bit vector and iterate.
                    self._cache(new_bvec, new_logp)
                    yield new_bvec, new_logp
                    new_bvec, new_logp = next(new, (None, None))

    bvs = None
    # Iterate through the positions in increasing order of the absolute
    # difference of the log probabilities of the two states of the bit.
    # In this order, the first positions encountered have probabilities
    # close to 0.5 (such that flipping the bit has the smallest effect
    # on the probability of obtaining the bit vector), while the latest
    # positions encountered have probabilities close to 0 or 1.
    for index in np.argsort(diff_logp):
        # Determine the absolute difference in log probability between
        # the two states for the bit at the current index.
        index_logp = diff_logp[index]
        if np.isfinite(index_logp):
            # Create a bit vector series for this index iff the absolute
            # difference in log probability is not infinite, i.e. the
            # probability is neither 0 nor 1.
            bvs = BitVectorSeries(bit=index, logp=index_logp,
                                  sub=((lambda: [(base_bits, base_logp)])
                                       if bvs is None else bvs.iter))
    return iter(bvs.iter() if bvs is not None else list())
