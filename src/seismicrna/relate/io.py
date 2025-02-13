from abc import ABC
from collections import defaultdict
from typing import Any, Iterable

import numpy as np

from .batch import ReadNamesBatch, RelateBatch
from ..core import path
from ..core.io import MutsBatchIO, ReadBatchIO, RefIO
from ..core.logs import logger
from ..core.seq import DNA, Region
from ..core.types import fit_uint_type


class RelateIO(RefIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.RELATE_STEP}


class ReadNamesBatchIO(ReadBatchIO, RelateIO, ReadNamesBatch):

    @classmethod
    def file_seg_type(cls):
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


class RelateBatchIO(MutsBatchIO, RelateIO, RelateBatch):

    @classmethod
    def file_seg_type(cls):
        return path.RelateBatSeg


def from_reads(reads: Iterable[tuple[str,
                                     tuple[tuple[list[int], list[int]],
                                           dict[int, int]]]],
               sample: str,
               ref: str,
               refseq: DNA,
               batch: int,
               write_read_names: bool):
    """ Accumulate reads into relation vectors. """
    # Initialize empty data.
    names = list()
    seg_end5s = list()
    seg_end3s = list()
    muts = {pos: defaultdict(list) for pos in range(1, len(refseq) + 1)}
    # Collect the mutation data from the reads.
    for (name, ((end5s, end3s), poss)) in reads:
        if all(end5 > end3 for end5, end3 in zip(end5s, end3s, strict=True)):
            # Skip a read if no segment has any coverage.
            logger.warning(f"Skipped read {repr(name)} with 5' end(s) {end5s} "
                           f"> 3' end(s) {end3s}")
        else:
            # Otherwise, add the data for the read to the batch.
            read = len(names)
            names.append(name)
            seg_end5s.append(end5s)
            seg_end3s.append(end3s)
            for pos, rel in poss.items():
                muts[pos][rel].append(read)
    # TODO: Check segment number is consistent between batches
    max_segs = max(max(len(segs) for segs in seg_end5s), max(len(segs) for segs in seg_end3s))
    seg_end5s = [segs + [segs[0]] * (max_segs - len(segs)) for segs in seg_end5s]
    seg_end3s = [segs + [segs[0]] * (max_segs - len(segs)) for segs in seg_end3s]

    # Make sure seg_end5s and seg_end3s have two dimensions and at least
    # one column each.
    pos_dtype = fit_uint_type(max(muts))
    seg_end5s = np.array(seg_end5s, dtype=pos_dtype)
    if seg_end5s.ndim < 2:
        seg_end5s = seg_end5s[:, np.newaxis]
    seg_end3s = np.array(seg_end3s, dtype=pos_dtype)
    if seg_end3s.ndim < 2:
        seg_end3s = seg_end3s[:, np.newaxis]
    # Assemble and return the batches.
    relate_batch = RelateBatchIO(sample=sample,
                                 batch=batch,
                                 region=Region(ref, refseq),
                                 seg_end5s=seg_end5s,
                                 seg_end3s=seg_end3s,
                                 muts=muts)
    if write_read_names:
        name_batch = ReadNamesBatchIO(sample=sample,
                                      ref=ref,
                                      batch=batch,
                                      names=names)
    else:
        name_batch = None
    return relate_batch, name_batch
