from functools import cached_property

import numpy as np

from ..core.batch import AllReadBatch, MutsBatch, get_length, sanitize_pos
from ..core.seq import POS_INDEX


class QnamesBatch(AllReadBatch):

    def __init__(self, *, names: list[str] | np.ndarray, **kwargs):
        super().__init__(**kwargs)
        self.names = np.asarray(names, dtype=str)

    @cached_property
    def num_reads(self):
        return get_length(self.names, "read names")


class RelateBatch(MutsBatch, AllReadBatch):

    def __init__(self, muts: dict[int, dict], seqlen: int, **kwargs):
        super().__init__(muts={pos: muts[pos]
                               for pos in sanitize_pos(muts, seqlen)},
                         seqlen=seqlen,
                         **kwargs)

    @cached_property
    def num_reads(self):
        nums_reads = list(set(map(get_length, (self.end5s,
                                               self.mid5s,
                                               self.mid3s,
                                               self.end3s))))
        if len(nums_reads) != 1:
            raise ValueError(f"Got inconsistent numbers of reads: {nums_reads}")
        return nums_reads[0]

    @cached_property
    def pos_nums(self):
        return np.arange(POS_INDEX,
                         POS_INDEX + self.max_pos,
                         dtype=self.pos_dtype)

    @cached_property
    def pattern(self):
        return None
