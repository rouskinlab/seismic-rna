from abc import ABC
from collections import defaultdict
from functools import cached_property
from typing import Any, Iterable

import numpy as np

from .batch import ReadNamesBatch, RelateMutsBatch, RelateRegionMutsBatch
from ..core import path
from ..core.array import calc_inverse
from ..core.io import (MutsBatchIO,
                       ReadBatchIO,
                       RefBrickleIO,
                       RefFileIO)
from ..core.logs import logger
from ..core.seq import DNA, Region
from ..core.types import fit_uint_type
from ..core.validate import require_isinstance


class RelateFile(path.HasRefFilePath, ABC):

    @classmethod
    def get_step(cls):
        return path.RELATE_STEP


class RelateIO(RelateFile, RefFileIO, ABC):
    pass


class RefseqIO(RefBrickleIO, RelateIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.RefseqFileSeg

    def __init__(self, *args, refseq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()


class ReadNamesBatchIO(ReadNamesBatch, ReadBatchIO, RefBrickleIO, RelateIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.ReadNamesBatSeg

    def __getstate__(self):
        state = super().__getstate__()
        # Compress the names from unicode to bytes before pickling.
        state["names"] = np.char.encode(self.names)
        return state

    def __setstate__(self, state: dict[str, Any]):
        super().__setstate__(state)
        # Decompress the names from bytes to unicode when unpickling.
        self.names = np.char.decode(state["names"])


class RelateBatchIO(RelateMutsBatch, MutsBatchIO, RefBrickleIO, RelateIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.RelateBatSeg

    @classmethod
    def from_region_batch(cls,
                          batch: RelateRegionMutsBatch, *,
                          sample: str,
                          branches: dict[str, str]):
        """ Create an instance from a RelateRegionMutsBatch. """
        require_isinstance("batch", batch, RelateRegionMutsBatch)
        return cls(batch=batch.batch,
                   region=batch.region,
                   seg_end5s=batch.seg_end5s,
                   seg_end3s=batch.seg_end3s,
                   muts=batch.muts,
                   sample=sample,
                   branches=branches,
                   sanitize=False)

    @classmethod
    def simulate(cls, *args, sample: str, branches: dict[str, str], **kwargs):
        # This class does not have all methods needed to simulate data,
        # so steal them from RelateRegionMutsBatch.
        return cls.from_region_batch(
            RelateRegionMutsBatch.simulate(*args, **kwargs),
            sample=sample,
            branches=branches
        )

    def to_region_batch(self, region: Region):
        """ Create a RelateRegionMutsBatch from this instance. """
        require_isinstance("region", region, Region)
        return RelateRegionMutsBatch(batch=self.batch,
                                     seg_end5s=self.seg_end5s,
                                     seg_end3s=self.seg_end3s,
                                     muts=self.muts,
                                     region=region,
                                     sanitize=False)


def from_reads(reads: Iterable[tuple[str,
                                     tuple[tuple[list[int], list[int]],
                                           dict[int, int]]]], *,
               sample: str,
               branches: dict[str, str],
               ref: str,
               refseq: DNA,
               batch: int,
               write_read_names: bool,
               drop_empty_reads: bool = True):
    """ Gather reads into a batch of relationships. """
    logger.routine(f"Began gathering reads into batch {batch}")
    max_segs = -1
    diff_segs = False
    names = list()
    seg_end5s = list()
    seg_end3s = list()
    muts = {pos: defaultdict(list) for pos in range(1, len(refseq) + 1)}
    pos_dtype = fit_uint_type(len(refseq))
    # Accumulate the data for the remaining reads.
    for read, (name, ((end5s, end3s), poss)) in enumerate(reads):
        # Find the number of segments in this read.
        num_segs = len(end5s)
        if len(end3s) != num_segs:
            raise ValueError(f"Read {repr(name)} has {num_segs} 5' ends "
                             f"and {len(end3s)} 3' ends")
        # It is unlikely that a read has a different number of segments
        # than other reads, so for efficiency, check for equality first.
        if num_segs != max_segs:
            # If the number of segments is not equal, then check whether
            # this read has more or fewer segments than other reads.
            if num_segs > max_segs:
                # If this read is the first read, then max_segs == -1,
                # so consider this read to have a different number of
                # segments from other reads only if max_segs >= 0.
                if max_segs >= 0:
                    logger.detail(f"Read {repr(name)} has {num_segs} segments, "
                                  f"more than {max_segs} in previous reads")
                    diff_segs = True
                max_segs = num_segs
            else:
                logger.detail(f"Read {repr(name)} has {num_segs} segments, "
                              f"fewer than {max_segs} in previous reads")
                diff_segs = True
        names.append(name)
        seg_end5s.append(end5s)
        seg_end3s.append(end3s)
        for pos, rel in poss.items():
            muts[pos][rel].append(read)
    if max_segs < 0:
        # If there were no reads, then set max_segs to 0 by default.
        max_segs = 0
    # Ensure all reads have the same number of segments, padding with
    # extra no-coverage segments (end5=1, end3=0) as necessary.
    logger.detail(f"The most segments in one read is {max_segs}")
    if diff_segs:
        for end5s, end3s in zip(seg_end5s, seg_end3s, strict=True):
            num_segs = len(end5s)
            assert len(end3s) == num_segs
            if num_segs < max_segs:
                padding = max_segs - num_segs
                end5s.extend([1] * padding)
                end3s.extend([0] * padding)
        logger.detail(f"Padded all reads to {max_segs} segments")
    # Convert seg_end5s and seg_end3s to NumPy arrays.
    num_reads = len(names)
    assert len(seg_end5s) == len(seg_end3s) == num_reads
    if num_reads > 0:
        seg_end5s = np.array(seg_end5s, dtype=pos_dtype)
        seg_end3s = np.array(seg_end3s, dtype=pos_dtype)
    else:
        seg_end5s = np.zeros((num_reads, max_segs), dtype=pos_dtype)
        seg_end3s = np.zeros((num_reads, max_segs), dtype=pos_dtype)
    assert seg_end5s.ndim == seg_end3s.ndim == 2
    assert seg_end5s.shape == seg_end3s.shape == (num_reads, max_segs)
    if drop_empty_reads:
        # Keep only reads that have at least 1 segment with length ≥ 1.
        keep = np.any(seg_end3s >= seg_end5s, axis=1)
        if not keep.all():
            # Calculate the numbers of the reads to keep.
            keep_read_nums = np.flatnonzero(keep)
            # Select the names and 5'/3' ends of the reads that are kept.
            names = [names[i] for i in keep_read_nums]
            seg_end5s = seg_end5s[keep_read_nums]
            seg_end3s = seg_end3s[keep_read_nums]
            # Map the read numbers in muts to the new read numbers; drop
            # the numbers of reads that were not kept.
            keep_read_inverse = calc_inverse(keep_read_nums, verify=False)
            for rels in muts.values():
                # Iterate over list(rels.items()) to avoid iterating
                # over rels.items() directly as rels is being modified.
                for rel, curr_read_nums in list(rels.items()):
                    # Find the new index for each read that is kept.
                    # This algorithm works the same as locate_elements()
                    # in core.array; it is reimplemented here to avoid
                    # needing to calculate keep_read_inverse repeatedly.
                    keep_curr_read_nums = np.intersect1d(curr_read_nums,
                                                         keep_read_nums,
                                                         assume_unique=True)
                    rels[rel] = keep_read_inverse[keep_curr_read_nums]
            logger.detail(
                f"Dropped {num_reads - keep_read_nums.size} empty reads"
            )
    # Assemble and return the batches.
    relate_batch = RelateBatchIO(sample=sample,
                                 branches=branches,
                                 batch=batch,
                                 region=Region(ref, refseq),
                                 seg_end5s=seg_end5s,
                                 seg_end3s=seg_end3s,
                                 muts=muts)
    if write_read_names:
        name_batch = ReadNamesBatchIO(sample=sample,
                                      branches=branches,
                                      ref=ref,
                                      batch=batch,
                                      names=names)
    else:
        name_batch = None
    logger.routine(f"Ended gathering reads into batch {batch}")
    return relate_batch, name_batch
