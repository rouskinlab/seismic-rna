from abc import ABC
from collections import defaultdict
from functools import cached_property
from typing import Any, Iterable

import numpy as np

from .batch import ReadNamesBatch, RelateBatch
from ..core import path
from ..core.array import calc_inverse
from ..core.io import (MutsBatchIO,
                       ReadBatchIO,
                       RefBrickleIO,
                       RefFileIO)
from ..core.logs import logger
from ..core.seq import DNA, Region
from ..core.types import fit_uint_type


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


class RelateBatchIO(RelateBatch, MutsBatchIO, RefBrickleIO, RelateIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.RelateBatSeg


def from_reads(reads: Iterable[tuple[str,
                                     tuple[tuple[list[int], list[int]],
                                           dict[int, int]]]], *,
               sample: str,
               branches: dict[str, str],
               ref: str,
               refseq: DNA,
               batch: int,
               write_read_names: bool):
    """ Accumulate reads into relation vectors. """
    reads = iter(reads)
    muts = {pos: defaultdict(list) for pos in range(1, len(refseq) + 1)}
    try:
        name, ((end5s, end3s), poss) = next(reads)
    except StopIteration:
        # There are no reads.
        min_segs = 0
        max_segs = 0
        names = list()
        seg_end5s = list()
        seg_end3s = list()
        logger.detail("There are no reads")
    else:
        # Initialize with the first read.
        num_segs = len(end5s)
        if len(end3s) != num_segs:
            raise ValueError(f"Read {repr(name)} has {num_segs} 5' ends "
                             f"and {len(end3s)} 3' ends")
        min_segs = num_segs
        max_segs = num_segs
        names = [name]
        seg_end5s = [end5s]
        seg_end3s = [end3s]
        for pos, rel in poss.items():
            muts[pos][rel].append(0)
        logger.detail(f"Read {repr(name)} has {num_segs} segments")
    # Accumulate the data for the remaining reads.
    for (name, ((end5s, end3s), poss)) in reads:
        # Find the number of segments and update the min/max numbers.
        num_segs = len(end5s)
        if len(end3s) != num_segs:
            raise ValueError(f"Read {repr(name)} has {num_segs} 5' ends "
                             f"and {len(end3s)} 3' ends")
        if num_segs < min_segs:
            logger.detail(
                f"Read {repr(name)} has {num_segs} segments (new minimum)"
            )
            min_segs = num_segs
        if num_segs > max_segs:
            logger.detail(
                f"Read {repr(name)} has {num_segs} segments (new maximum)"
            )
            max_segs = num_segs
        read = len(names)
        names.append(name)
        seg_end5s.append(end5s)
        seg_end3s.append(end3s)
        for pos, rel in poss.items():
            muts[pos][rel].append(read)
    # Ensure all reads have the same number of segments, padding with
    # extra no-coverage segments (end5=1, end3=0) as necessary.
    logger.detail(f"All reads have ≥ {min_segs} and ≤ {max_segs} segments")
    assert min_segs <= max_segs
    if min_segs < max_segs:
        for end5s, end3s in zip(seg_end5s, seg_end3s, strict=True):
            num_segs = len(end5s)
            assert len(end3s) == num_segs
            if num_segs < max_segs:
                padding = max_segs - num_segs
                end5s.extend([1] * padding)
                end3s.extend([0] * padding)
        logger.detail(f"Padded reads with < {max_segs} segments")
    # Make sure seg_end5s and seg_end3s have two dimensions and at least
    # 1 column each (if there are 0 reads, then seg_end5s and seg_end3s
    # will initially be 1-dimensional when converted to arrays).
    pos_dtype = fit_uint_type(max(muts))
    seg_end5s = np.array(seg_end5s, dtype=pos_dtype)
    if seg_end5s.ndim < 2:
        seg_end5s = seg_end5s[:, np.newaxis]
    assert seg_end5s.ndim == 2
    seg_end3s = np.array(seg_end3s, dtype=pos_dtype)
    if seg_end3s.ndim < 2:
        seg_end3s = seg_end3s[:, np.newaxis]
    assert seg_end3s.ndim == 2
    assert seg_end5s.shape == seg_end3s.shape == (len(names), max_segs)
    # Keep only reads that have at least 1 segment with length ≥ 1.
    keep = np.any(seg_end3s >= seg_end5s, axis=1)
    if not keep.all():
        # Calculate the numbers of the reads to keep.
        keep_read_nums = np.flatnonzero(keep)
        # Select only the 5' and 3' ends of the reads that are kept.
        seg_end5s = seg_end5s[keep_read_nums]
        seg_end3s = seg_end3s[keep_read_nums]
        # Map the read numbers in muts to the new read numbers; drop the
        # numbers of reads that were not kept.
        keep_read_inverse = calc_inverse(keep_read_nums, verify=False)
        for rels in muts.values():
            # Iterate over list(rels.items()) to avoid iterating over
            # rels.items() directly as rels is being modified.
            for rel, curr_read_nums in list(rels.items()):
                # Find the new index for each current read that is kept.
                # This algorithm works the same as locate_elements() in
                # core.array; it is reimplemented here to avoid needing
                # to calculate keep_read_inverse repeatedly.
                keep_curr_read_nums = np.intersect1d(curr_read_nums,
                                                     keep_read_nums,
                                                     assume_unique=True)
                rels[rel] = keep_read_inverse[keep_curr_read_nums]
        logger.detail(f"Removed {len(names) - keep_read_nums.size} reads "
                      f"containing 0 segments with positive lengths")
        # Select the names last so that before this, len(names) is still
        # the original number of reads.
        names = [names[i] for i in keep_read_nums]
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
    return relate_batch, name_batch
