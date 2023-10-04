from abc import ABC
from collections import defaultdict
from functools import cached_property
from typing import Any, Iterable

import numpy as np
import pandas as pd

from ..core import path
from ..core.batch import ReadBatch, MutsBatch, POS_INDEX, READ_INDEX
from ..core.cmd import CMD_REL
from ..core.files import RefOutput
from ..core.rel import MATCH, NOCOV
from ..core.seq import DNA


class RelateOutput(RefOutput, ABC):

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_REL}


def get_num_pos(max_pos: int):
    return max_pos - POS_INDEX + 1


def get_max_read(num_reads: int):
    return num_reads + READ_INDEX - 1


def get_read_nums(num_reads: int, dtype: type):
    return np.arange(READ_INDEX, num_reads + READ_INDEX, dtype=dtype)


class QnamesBatch(ReadBatch, RelateOutput):

    @classmethod
    def file_seg_type(cls):
        return path.QnamesBatSeg

    def __init__(self, *args, names: list[str] | np.ndarray, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = np.asarray(names, dtype=str)

    @property
    def num_reads(self):
        return self.names.size

    @property
    def max_read(self):
        return get_max_read(self.num_reads)

    @cached_property
    def read_nums(self):
        return get_read_nums(self.num_reads, self.read_dtype)

    @cached_property
    def array(self):
        return self.names

    @cached_property
    def table(self):
        return pd.Series(self.array, self.read_nums, copy=False)

    def __getstate__(self):
        state = super().__getstate__()
        # Compress the names from unicode to bytes before pickling.
        state["names"] = np.char.encode(self.names)
        return state

    def __setstate__(self, state: dict[str, Any]):
        super().__setstate__(state)
        # Decompress the names from bytes to unicode when unpickling.
        self.names = np.char.decode(state["names"])


class RelateBatch(MutsBatch, RelateOutput):

    @classmethod
    def file_seg_type(cls):
        return path.RelateBatSeg

    @property
    def num_pos(self):
        return get_num_pos(self.max_pos)

    @property
    def num_reads(self):
        nums_reads = list(set(x.size for x in (self.end5s,
                                               self.mid5s,
                                               self.mid3s,
                                               self.end3s)))
        if len(nums_reads) != 1:
            raise ValueError(
                f"Need exactly one number of reads, but got {nums_reads}")
        return nums_reads[0]

    @property
    def max_read(self):
        return get_max_read(self.num_reads)

    @cached_property
    def pos(self):
        return np.arange(POS_INDEX,
                         self.num_pos + POS_INDEX,
                         dtype=self.pos_dtype)

    @cached_property
    def read_nums(self):
        return get_read_nums(self.num_reads, self.read_dtype)

    @cached_property
    def array(self):
        """ Matrix of vectors in use. """
        array = np.full((self.num_reads, self.num_pos), NOCOV)
        # Mark the covered positions as matches.
        for read, (end5i, mid5i, mid3p, end3p) in enumerate(zip(self.end5s - 1,
                                                                self.mid5s - 1,
                                                                self.mid3s,
                                                                self.end3s,
                                                                strict=True)):
            array[read, end5i: mid3p] = MATCH
            array[read, mid5i: end3p] = MATCH
        # Mark the mutated positions.
        for pos, rels in self.muts.items():
            for rel, reads in rels.items():
                array[reads, pos - 1] = rel
        return array


def from_reads(reads: Iterable[tuple[str, int, int, int, int, dict[int, int]]],
               sample: str,
               ref: str,
               refseq: DNA,
               batch: int):
    """ Accumulate reads into relation vectors. """
    names = list()
    end5s = list()
    mid5s = list()
    mid3s = list()
    end3s = list()
    muts = {pos: defaultdict(list)
            for pos in range(POS_INDEX, len(refseq) + POS_INDEX)}
    for read, (name, end5, mid5, mid3, end3, poss) in enumerate(reads,
                                                                READ_INDEX):
        names.append(name)
        end5s.append(end5)
        mid5s.append(mid5)
        mid3s.append(mid3)
        end3s.append(end3)
        for pos, mut in poss.items():
            muts[pos][mut].append(read)
    name_batch = QnamesBatch(batch, sample, ref, names=names)
    rel_batch = RelateBatch(batch,
                            sample,
                            ref,
                            refseq=refseq,
                            end5s=end5s,
                            mid5s=mid5s,
                            mid3s=mid3s,
                            end3s=end3s,
                            muts=muts)
    return name_batch, rel_batch
