from functools import cached_property
from typing import Iterable

import numpy as np

from ..core.batch import Batch, RelsBatch, get_length, sanitize_pos

from ..core.types import fit_uint_type


class MaskBatch(Batch):

    def __init__(self, *, read_nums: np.ndarray, max_read: int, **kwargs):
        self._max_read = max_read
        super().__init__(**kwargs)
        if read_nums.size > 0 and np.max(read_nums) > max_read:
            raise ValueError(f"All read numbers must be ≤ {max_read}, but got "
                             f"{read_nums[read_nums > max_read]}")
        self._read_nums = np.asarray(read_nums, fit_uint_type(max_read))

    @property
    def read_nums(self):
        return self._read_nums

    @property
    def num_reads(self):
        return get_length(self.read_nums, "read_nums")

    @property
    def max_read(self):
        return self._max_read

    @cached_property
    def read_idx(self):
        read_map = np.zeros(self.max_read + 1, dtype=self.read_dtype)
        read_map[self.read_nums] = np.arange(self.num_reads)
        return read_map


class MaskRelsBatch(MaskBatch, RelsBatch):

    @classmethod
    def from_batch(cls,
                   batch: RelsBatch,
                   positions: Iterable[int] | None = None,
                   reads: Iterable[int] | None = None):
        # Clean and validate the selection.
        if positions is not None:
            positions = sanitize_pos(positions, batch.max_pos)
        # Select mutations at each position.
        muts = dict()
        for pos in positions if positions is not None else batch.pos_nums:
            muts[pos] = dict()
            # Select reads with each type of mutation at this position.
            for mut, pos_mut_reads in batch.muts.get(pos, dict()).items():
                muts[pos][mut] = (np.intersect1d(pos_mut_reads,
                                                 reads,
                                                 assume_unique=True)
                                  if reads is not None
                                  else pos_mut_reads)
        if reads is not None:
            read_nums = np.asarray(reads, dtype=batch.read_dtype)
            read_indexes = batch.read_idx[read_nums]
            end5s = batch.end5s[read_indexes]
            mid5s = batch.mid5s[read_indexes]
            mid3s = batch.mid3s[read_indexes]
            end3s = batch.end3s[read_indexes]
        else:
            read_nums = batch.read_nums
            end5s = batch.end5s
            mid5s = batch.mid5s
            mid3s = batch.mid3s
            end3s = batch.end3s
        return cls(batch=batch.batch,
                   muts=muts,
                   seqlen=batch.seqlen,
                   end5s=end5s,
                   mid5s=mid5s,
                   mid3s=mid3s,
                   end3s=end3s,
                   read_nums=read_nums,
                   max_read=batch.max_read)

    @cached_property
    def pos_nums(self):
        return np.array(list(self.muts))
