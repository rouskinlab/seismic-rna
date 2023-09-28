from collections import defaultdict
from functools import cached_property
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .output import RelateOutput
from ..core import path
from ..core.batch import Batch, POS_INDEX, READ_INDEX
from ..core.rel import MATCH, NOCOV, REL_TYPE
from ..core.seq import DNA


def get_max_read(num_reads: int):
    return num_reads + READ_INDEX - 1


def get_read_nums(num_reads: int, dtype: type):
    return np.arange(READ_INDEX, num_reads + READ_INDEX, dtype=dtype)


def get_pos(num_pos: int, dtype: type):
    return np.arange(POS_INDEX, num_pos + POS_INDEX, dtype=dtype)


class NameBatch(Batch, RelateOutput):

    @classmethod
    def file_seg_type(cls):
        return path.QnamesBatSeg

    def __init__(self, *args, names: list[str] | np.ndarray, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = np.asarray(names, dtype=str)

    @property
    def num_pos(self):
        return len(self.refseq)

    @property
    def num_reads(self):
        return self.names.size

    @property
    def max_read(self):
        return get_max_read(self.num_reads)

    @cached_property
    def pos(self):
        return get_pos(self.num_pos, self.pos_dtype)

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
        # Compress the names from unicode to bytes.
        state["names"] = np.char.encode(self.names)
        return state

    def __setstate__(self, state: dict[str, Any]):
        super().__setstate__(state)
        # Decompress the names from bytes to unicode.
        self.names = np.char.decode(state["names"])


class RelateBatch(Batch, RelateOutput):

    @classmethod
    def file_seg_type(cls):
        return path.RelateBatSeg

    def __init__(self,
                 *args,
                 end5s: list[int] | np.ndarray,
                 mid5s: list[int] | np.ndarray,
                 mid3s: list[int] | np.ndarray,
                 end3s: list[int] | np.ndarray,
                 muts: dict[int, dict[int, list[int] | np.ndarray]],
                 **kwargs):
        super().__init__(*args, **kwargs)
        # 5' ends
        self.end5s = np.asarray(end5s, dtype=self.pos_dtype)
        if self.end5s.ndim != 1:
            raise ValueError(
                f"end5s must have 1 dimension, but got {self.end5s.ndim}")
        if np.any(self.end5s < POS_INDEX):
            errors = self.end5s[self.end5s < POS_INDEX]
            raise ValueError(f"Got end5 < {POS_INDEX} for reads {errors}")
        # 3' ends
        self.end3s = np.asarray(end3s, dtype=self.pos_dtype)
        if self.end3s.ndim != 1:
            raise ValueError(
                f"end3s must have 1 dimension, but got {self.end3s.ndim}")
        if np.any(self.end3s < self.end5s):
            errors = self.end3s[self.end3s < self.end5s]
            raise ValueError(f"Got end3 < end5 for reads {errors}")
        # 5' mids
        self.mid5s = np.asarray(mid5s, dtype=self.pos_dtype)
        if self.mid5s.ndim != 1:
            raise ValueError(
                f"mid5s must have 1 dimension, but got {self.mid5s.ndim}")
        if np.any(self.mid5s < self.end5s):
            errors = self.mid5s[self.mid5s < self.end5s]
            raise ValueError(f"Got mid5 < end5 for reads {errors}")
        # 3' mids
        self.mid3s = np.asarray(mid3s, dtype=self.pos_dtype)
        if self.mid3s.ndim != 1:
            raise ValueError(
                f"mid3s must have 1 dimension, but got {self.mid3s.ndim}")
        if np.any(self.mid3s > self.end3s):
            errors = self.mid3s[self.mid3s > self.end3s]
            raise ValueError(f"Got mid3 > end3 for reads {errors}")
        # Mutations
        self.muts = {
            pos: {
                REL_TYPE(rel): np.asarray(reads, self.read_dtype)
                for rel, reads in muts.get(pos, dict()).items()
            }
            for pos in get_pos(self.num_pos, self.pos_dtype)
        }

    @property
    def num_pos(self):
        return len(self.refseq)

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
        return self.num_reads + READ_INDEX - 1

    @cached_property
    def pos(self):
        return get_pos(self.num_pos, self.pos_dtype)

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
    read = READ_INDEX
    for name, end5, mid5, mid3, end3, poss in reads:
        names.append(name)
        end5s.append(end5)
        mid5s.append(mid5)
        mid3s.append(mid3)
        end3s.append(end3)
        for pos, mut in poss.items():
            muts[pos][mut].append(read)
        read += 1
    name_batch = NameBatch(sample, ref, refseq, batch, names=names)
    rel_batch = RelateBatch(sample,
                            ref,
                            refseq,
                            batch,
                            end5s=end5s,
                            mid5s=mid5s,
                            mid3s=mid3s,
                            end3s=end3s,
                            muts=muts)
    return name_batch, rel_batch
