from abc import ABC
from collections import defaultdict
from typing import Any, Iterable

import numpy as np

from .batch import QnamesBatch, RelateBatch
from ..core import path
from ..core.batch import POS_INDEX
from ..core.io import MutsBatchIO, ReadBatchIO, RefIO
from ..core.seq import DNA


class RelateIO(RefIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_REL_DIR}


class QnamesBatchIO(RelateIO, ReadBatchIO, QnamesBatch):

    @classmethod
    def file_seg_type(cls):
        return path.QnamesBatSeg

    def __getstate__(self):
        state = super().__getstate__()
        # Compress the names from unicode to bytes before pickling.
        state["names"] = np.char.encode(self.names)
        return state

    def __setstate__(self, state: dict[str, Any]):
        super().__setstate__(state)
        # Decompress the names from bytes to unicode when unpickling.
        self.names = np.char.decode(state["names"])


class RelateBatchIO(RelateIO, MutsBatchIO, RelateBatch):

    @classmethod
    def file_seg_type(cls):
        return path.RelateBatSeg


def from_reads(reads: Iterable[tuple[str, int, int, int, int, dict[int, int]]],
               sample: str,
               ref: str,
               refseq: DNA,
               batch: int):
    """ Accumulate reads into relation vectors. """
    # Initialize empty data.
    names = list()
    end5s = list()
    mid5s = list()
    mid3s = list()
    end3s = list()
    muts = defaultdict(lambda: defaultdict(list))
    # Collect the mutation data from the reads.
    for read, (name, end5, mid5, mid3, end3, poss) in enumerate(reads):
        names.append(name)
        end5s.append(end5)
        mid5s.append(mid5)
        mid3s.append(mid3)
        end3s.append(end3)
        for pos, rel in poss.items():
            muts[pos][rel].append(read)
    # Validate the positions.
    if min(muts) < POS_INDEX:
        raise ValueError(f"All positions must be ≥ {POS_INDEX}, but got "
                         f"{[x for x in sorted(muts) if x < POS_INDEX]}")
    if max(muts) > len(refseq):
        raise ValueError(f"All positions must be ≤ {len(refseq)}, but got "
                         f"{[x for x in sorted(muts) if x > len(refseq)]}")
    # Assemble and return the batches.
    name_batch = QnamesBatchIO(sample=sample,
                               ref=ref,
                               batch=batch,
                               names=names)
    rel_batch = RelateBatchIO(sample=sample,
                              ref=ref,
                              batch=batch,
                              end5s=end5s,
                              mid5s=mid5s,
                              mid3s=mid3s,
                              end3s=end3s,
                              muts=muts,
                              seqlen=len(refseq))
    return name_batch, rel_batch
