from abc import ABC, abstractmethod
from functools import cached_property
from typing import Any, Iterable

import numpy as np

from ..core.array import locate_elements
from ..core.batch import MutsBatch, match_reads_segments
from ..core.data import WideDataset
from ..core.seq import Section

BATCH_NUM = "batch"
READ_NUMS = "read_nums"
SEG_END5S = "seg_end5s"
SEG_END3S = "seg_end3s"
MUTS = "muts"
RESPS = "resps"


def _join_position(muts: dict[int, dict[int, np.ndarray]],
                   add_muts: dict[int, dict[int, np.ndarray]],
                   position: int):
    if pos_muts := muts.get(position):
        if add_pos_muts := add_muts.get(position):
            joined_pos_muts = dict()
            muts = set(pos_muts)
            add_muts = set(add_pos_muts)
            # For types of mutations shared by both pos_muts and add_pos_muts,
            # take the union of read numbers and then remove duplicates.
            for mut in muts & add_muts:
                joined_pos_muts[mut] = np.unique(np.concatenate(
                    [pos_muts[mut], add_pos_muts[mut]]
                ))
            # For types of mutations unique to pos_muts or to add_pos_muts,
            # merely copy the read numbers into joined_pos_muts.
            for mut in muts - add_muts:
                joined_pos_muts[mut] = pos_muts[mut]
            for mut in add_muts - muts:
                joined_pos_muts[mut] = add_pos_muts[mut]
            # Verify that no read has more than one type of mutation.
            reads, counts = np.unique(
                np.concatenate(list(joined_pos_muts.values())),
                return_counts=True
            )
            if counts.max(initial=0) > 1:
                raise ValueError(f"Reads {reads[counts > 1]} have > 1 type of "
                                 f"mutation at position {position}")
            return joined_pos_muts
        return pos_muts
    return add_muts.get(position, dict())


def _join_attrs(attrs: dict[str, Any],
                add_attrs: dict[str, Any],
                section: Section):
    # Verify the batch number matches (but no need to update).
    if attrs[BATCH_NUM] != add_attrs[BATCH_NUM]:
        raise ValueError(f"Inconsistent batch number ({attrs[BATCH_NUM]} "
                         f"≠ {add_attrs[BATCH_NUM]})")
    # Merge the read numbers and find the indexes of the original read
    # numbers in the combined array.
    union_read_nums = np.union1d(attrs[READ_NUMS], add_attrs[READ_NUMS])
    read_indexes, add_read_indexes = locate_elements(union_read_nums,
                                                     attrs[READ_NUMS],
                                                     add_attrs[READ_NUMS],
                                                     what="read_nums",
                                                     verify=False)
    attrs[READ_NUMS] = union_read_nums
    # Merge the end coordinates, setting the 5' and 3' ends of missing
    # segments to one after the 3' end and one before the 5' end of the
    # section, respectively (the 5'/3' coordinates are swapped so that
    # any missing segments will have 5' coordinates greater than their
    # 3' coordinates and thus be masked).
    num_reads, = attrs[READ_NUMS].shape
    _, num_segs = match_reads_segments(attrs[SEG_END5S],
                                       attrs[SEG_END3S])
    _, add_num_segs = match_reads_segments(add_attrs[SEG_END5S],
                                           add_attrs[SEG_END3S])
    seg_end5s = np.full((num_reads, num_segs), section.end3 + 1)
    seg_end3s = np.full((num_reads, num_segs), section.end5 - 1)
    add_seg_end5s = np.full((num_reads, add_num_segs), section.end3 + 1)
    add_seg_end3s = np.full((num_reads, add_num_segs), section.end5 - 1)
    seg_end5s[read_indexes] = attrs[SEG_END5S]
    seg_end3s[read_indexes] = attrs[SEG_END3S]
    add_seg_end5s[add_read_indexes] = add_attrs[SEG_END5S]
    add_seg_end3s[add_read_indexes] = add_attrs[SEG_END3S]
    attrs[SEG_END5S] = np.hstack([seg_end5s, add_seg_end5s])
    attrs[SEG_END3S] = np.hstack([seg_end3s, add_seg_end3s])
    # Merge the mutations.
    muts = attrs[MUTS]
    add_muts = add_attrs[MUTS]
    attrs[MUTS] = {pos: _join_position(muts, add_muts, pos)
                   for pos in section.unmasked_int}


class JoinMutsDataset(WideDataset, ABC):

    @classmethod
    @abstractmethod
    def get_batch_type(cls) -> type[MutsBatch]:
        """ Type of batch. """

    @classmethod
    @abstractmethod
    def name_batch_attrs(cls) -> list[str]:
        """ Name the attributes of each batch. """

    @classmethod
    def check_batch_type(cls, batch: MutsBatch):
        """ Raise TypeError if the batch is the incorrect type. """
        if not isinstance(batch, cls.get_batch_type()):
            raise TypeError(f"Expected {cls.get_batch_type().__name__}, "
                            f"but got {type(batch).__name__}")

    @classmethod
    def _get_first_batch(cls, batches: Iterable[tuple[str, MutsBatch]]):
        """ Get the first batch; raise ValueError if no batches. """
        for sect, batch in batches:
            cls.check_batch_type(batch)
            return batch, sect
        raise ValueError("Cannot get first batch among 0 batches")

    @cached_property
    def min_mut_gap(self):
        return self._get_common_attr("min_mut_gap")

    def _get_batch_attrs(self, batch: MutsBatch, _: str):
        """ Get the values of the attributes from a batch. """
        self.check_batch_type(batch)
        return {attr: getattr(batch, attr) for attr in self.name_batch_attrs()}

    def _join_attrs(self, attrs: dict[str, Any], add_attrs: dict[str, Any]):
        """ Join the attributes from a new batch. """
        _join_attrs(attrs, add_attrs, self.section)

    def _finalize_attrs(self, attrs: dict[str, Any]):
        """ Finalize the attributes. """

    def _join(self, batches: Iterable[tuple[str, MutsBatch]]):
        attrs = self._get_batch_attrs(*self._get_first_batch(batches))
        for sect, batch in batches:
            self._join_attrs(attrs, self._get_batch_attrs(batch, sect))
        self._finalize_attrs(attrs)
        return self.get_batch_type()(section=self.section, **attrs)

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
