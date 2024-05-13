from abc import ABC
from collections import defaultdict
from typing import Any, Iterable

import numpy as np

from .batch import QnamesBatch, RelateBatch
from ..core import path
from ..core.io import MutsBatchIO, ReadBatchIO, RefIO
from ..core.seq import DNA, Section
from ..core.types import fit_uint_type


class RelateIO(RefIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_REL_DIR}


class QnamesBatchIO(ReadBatchIO, RelateIO, QnamesBatch):

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


class RelateBatchIO(MutsBatchIO, RelateIO, RelateBatch):

    @classmethod
    def file_seg_type(cls):
        return path.RelateBatSeg


def from_reads(reads: Iterable[tuple[str, tuple[list[int], [list[int]]], dict[int, int]]],
               sample: str,
               ref: str,
               refseq: DNA,
               batch: int):
    """ Accumulate reads into relation vectors. """
    # Initialize empty data.
    names = list()
    seg_end5s = list()
    seg_end3s = list()
    muts = {pos: defaultdict(list) for pos in range(1, len(refseq) + 1)}
    # Collect the mutation data from the reads.
    for read, (name, (end5s, end3s), poss) in enumerate(reads):
        names.append(name)
        seg_end5s.append(end5s)
        seg_end3s.append(end3s)
        for pos, rel in poss.items():
            muts[pos][rel].append(read)
    # Assemble and return the batches.
    pos_dtype = fit_uint_type(max(muts))
    name_batch = QnamesBatchIO(sample=sample,
                               ref=ref,
                               batch=batch,
                               names=names)
    rel_batch = RelateBatchIO(sample=sample,
                              batch=batch,
                              section=Section(ref, refseq),
                              seg_end5s=np.array(seg_end5s, dtype=pos_dtype),
                              seg_end3s=np.array(seg_end3s, dtype=pos_dtype),
                              muts=muts)
    return name_batch, rel_batch

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
