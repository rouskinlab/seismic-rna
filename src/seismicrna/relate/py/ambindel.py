from __future__ import annotations

from abc import ABC, abstractmethod

from .encode import encode_relate
from ...core.rel import DELET, INS_5, INS_3, MATCH, SUB_N


def get_ins_rel(insert3: bool):
    return INS_3 if insert3 else INS_5


def calc_lateral5(lateral3: int):
    return lateral3 - 1


def get_lateral(lateral3: int, insert3: int):
    return lateral3 if insert3 else calc_lateral5(lateral3)


def _check_out_of_bounds(next_opposite: int, end5: int, end3: int):
    """ Check if the position to move the indel is out of bounds. """
    return not end5 < next_opposite < end3


def _consistent_rels(rel1: int, rel2: int):
    """ Whether the relationships are consistent. """
    # The relationships are consistent if one of the following is true:
    # - There is at least one type of primary relationship (i.e. match,
    #   substitution to each base) that both of them could both be.
    # - Both could be substitutions to any base.
    return bool((rel1 & rel2) or (rel1 & SUB_N and rel2 & SUB_N))


class Indel(ABC):
    """ Insertion or Deletion. """

    __slots__ = ["opposite", "lateral3", "pod"]

    def __init__(self, opposite: int, lateral3: int, pod: IndelPod):
        self.opposite = opposite
        self.lateral3 = lateral3
        self.pod = pod
        self.pod.add(self)

    @property
    def lateral5(self):
        return calc_lateral5(self.lateral3)

    def get_lateral(self, insert3: bool):
        return get_lateral(self.lateral3, insert3)

    def _calc_positions(self, move5to3: bool):
        # Find the position within the indel's own sequence with which
        # to try to swap the indel.
        if move5to3:
            increment = 1
            swap_lateral = self.lateral3
        else:
            increment = -1
            swap_lateral = self.lateral5
        # If the indel moves, then its position within its own sequence
        # will change.
        next_lateral3 = self.lateral3 + increment
        # Find the position within the opposite sequence to which to try
        # to move the indel.
        next_opposite = self.opposite + increment
        # Keep incrementing until the position does not coincide with
        # any other indel of the same kind.
        while self.pod.get_indel_by_opp(next_opposite):
            next_opposite += increment
        return swap_lateral, next_lateral3, next_opposite

    def move(self, opposite: int, lateral3: int):
        """ Move the indel to a new position. """
        self.opposite = opposite
        self.lateral3 = lateral3

    @abstractmethod
    def _calc_rels(self,
                   ref_seq: str,
                   read_seq: str,
                   read_qual: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int) -> tuple[int, int]:
        """ Calculate the relationships of the swapping base with the
        current opposite and next opposite bases. """

    @abstractmethod
    def try_move(self,
                 rels: dict[int, int],
                 pods: list[IndelPod],
                 insert3: bool,
                 ref_seq: str,
                 read_seq: str,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 move5to3: bool) -> bool:
        """ Try to move the indel a step in one direction. """

    def __bool__(self):
        return True

    def __str__(self):
        return (f"{type(self).__name__} at {self.opposite} "
                f"{self.lateral5, self.lateral3}")


class Deletion(Indel):

    def _calc_rels(self,
                   ref_seq: str,
                   read_seq: str,
                   read_quals: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int):
        # Read base with which the deletion will swap.
        read_base = read_seq[swap_lateral - 1]
        read_qual = read_quals[swap_lateral - 1]
        # Reference base that is currently deleted from the read; will
        # move opposite that read base.
        ref_base_del = ref_seq[self.opposite - 1]
        # Reference base that is currently opposite that read base; will
        # become a deletion after the move.
        ref_base_opp = ref_seq[next_opposite - 1]
        rel_del = encode_relate(ref_base_del,
                                read_base,
                                read_qual,
                                min_qual)
        rel_opp = encode_relate(ref_base_opp,
                                read_base,
                                read_qual,
                                min_qual)
        return rel_del, rel_opp

    def try_move(self,
                 rels: dict[int, int],
                 pods: list[IndelPod],
                 insert3: bool,
                 ref_seq: str,
                 read_seq: str,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 move5to3: bool):
        (swap_lateral,
         next_lateral3,
         next_opposite) = self._calc_positions(move5to3)
        if _check_out_of_bounds(next_opposite, ref_end5, ref_end3):
            return False
        if _check_collisions(pods, self.pod, next_lateral3):
            return False
        rel_del, rel_opp = self._calc_rels(ref_seq,
                                           read_seq,
                                           read_qual,
                                           min_qual,
                                           swap_lateral,
                                           next_opposite)
        # Determine if the relationships before and after moving the
        # deletion are consistent with each other.
        if not _consistent_rels(rel_del, rel_opp):
            return False
        # When the deletion moves, the position from which it moves
        # gains the relationship (rel_indel) between the read base and
        # the reference base that was originally deleted from the read.
        rels[self.opposite] = rels.get(self.opposite, MATCH) | rel_del
        # Move the deletion.
        _move_indels(self, next_opposite, next_lateral3)
        # Mark the position to which the deletion moves.
        rels[self.opposite] = rels.get(self.opposite, MATCH) | DELET
        return True


class Insertion(Indel):

    def _calc_rels(self,
                   ref_seq: str,
                   read_seq: str,
                   read_quals: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int):
        # Reference base to whose position the insertion will move.
        ref_base = ref_seq[swap_lateral - 1]
        # Read base that is currently inserted into the reference; will
        # move opposite that reference base.
        read_base_ins = read_seq[self.opposite - 1]
        read_qual_ins = read_quals[self.opposite - 1]
        # Read base that is currently opposite that reference base; will
        # become an insertion after the move.
        read_base_opp = read_seq[next_opposite - 1]
        read_qual_opp = read_quals[next_opposite - 1]
        # Calculate the relationship between the reference base and each
        # base in the read.
        rel_ins = encode_relate(ref_base,
                                read_base_ins,
                                read_qual_ins,
                                min_qual)
        rel_opp = encode_relate(ref_base,
                                read_base_opp,
                                read_qual_opp,
                                min_qual)
        return rel_ins, rel_opp

    def try_move(self,
                 rels: dict[int, int],
                 pods: list[IndelPod],
                 insert3: bool,
                 ref_seq: str,
                 read_seq: str,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 move5to3: bool):
        (swap_lateral,
         next_lateral3,
         next_opposite) = self._calc_positions(move5to3)
        if _check_out_of_bounds(next_opposite, read_end5, read_end3):
            return False
        if _check_collisions(pods, self.pod, next_lateral3):
            return False
        rel_ins, rel_opp = self._calc_rels(ref_seq,
                                           read_seq,
                                           read_qual,
                                           min_qual,
                                           swap_lateral,
                                           next_opposite)
        # Determine if the relationships before and after moving the
        # insertion are consistent with each other.
        if not _consistent_rels(rel_ins, rel_opp):
            return False
        init_lateral = self.get_lateral(insert3)
        if move5to3 == insert3:
            # The insertion moves to the same side as it is marked.
            # Move the indels.
            _move_indels(self, next_opposite, next_lateral3)
            # Check after the move.
            if rel_ins != MATCH or all(ins.get_lateral(insert3) != init_lateral
                                       for ins in self.pod.indels):
                # If there is another insertion next to the base that just
                # took the place of the insertion that moved, then do not
                # mark the position of that base as a possible match.
                rels[self.get_lateral(not insert3)] = (
                        rels.get(self.get_lateral(not insert3), MATCH)
                        | rel_ins
                )
            rels[self.get_lateral(insert3)] = (
                    rels.get(self.get_lateral(insert3), MATCH)
                    | get_ins_rel(insert3)
            )
        else:
            # The insertion moves to the opposite side as it is marked.
            (swap_lateral_,
             next_lateral3_,
             next_opposite_) = self._calc_positions(insert3)
            rel_ins_, rel_opp_ = self._calc_rels(ref_seq,
                                                 read_seq,
                                                 read_qual,
                                                 min_qual,
                                                 swap_lateral_,
                                                 next_opposite_)
            # Move the indels.
            _move_indels(self, next_opposite, next_lateral3)
            # Check after the move.
            if rel_opp_ != MATCH or all(ins.get_lateral(insert3) != init_lateral
                                        for ins in self.pod.indels):
                # If there is another insertion next to the base that just
                # took the place of the insertion that moved, then do not
                # mark the position of that base as a possible match.
                rels[init_lateral] = rels.get(init_lateral, MATCH) | rel_opp_
            rels[self.get_lateral(insert3)] = (
                    rels.get(self.get_lateral(insert3), MATCH)
                    | get_ins_rel(insert3)
                    | (rel_ins if rel_ins != MATCH else 0)
            )
        return True


class IndelPod(ABC):
    __slots__ = ["n", "indels"]

    @classmethod
    @abstractmethod
    def indel_type(cls) -> type[Indel]:
        """ All indels in the pod must be of this type. """

    def __init__(self, n: int):
        self.n = n
        self.indels: list[Indel] = list()

    def add(self, indel: Indel):
        """ Add an indel to the pod, performing validation. """
        assert isinstance(indel, self.indel_type())
        for other in self.indels:
            assert indel.opposite != other.opposite
        self.indels.append(indel)

    def sort(self, move5to3: bool):
        """ Sort the indels in the pod by their positions. """
        self.indels.sort(key=(lambda indel: indel.opposite),
                         reverse=(not move5to3))

    def get_indel_by_opp(self, opposite: int):
        """ Return the indel that lies opposite the given position
        (opposite), or None if such an indel does not exist. """
        for indel in self.indels:
            if indel.opposite == opposite:
                return indel
        return None

    def __str__(self):
        return f"{type(self).__name__} {self.n}: {list(map(str, self.indels))}"


class DeletionPod(IndelPod):

    @classmethod
    def indel_type(cls):
        return Deletion


class InsertionPod(IndelPod):

    @classmethod
    def indel_type(cls):
        return Insertion


def _check_collisions(pods: list[IndelPod],
                      pod: IndelPod,
                      next_lateral3: int):
    """ Check if moving the indel will make it collide with another
    indel of the opposite kind. """
    pod_index = pods.index(pod)
    next_pod_index = pod_index + 1
    if 0 <= next_pod_index < len(pods):
        next_pod = pods[next_pod_index]
        assert not isinstance(next_pod, type(pod))
        next_indel = next_pod.indels[0]
        return (next_indel.opposite == next_lateral3
                or next_indel.opposite == calc_lateral5(next_lateral3))
    return False


def _move_indels(indel: Indel, opposite: int, lateral3: int):
    """ Move an indel while adjusting the positions of any other indels
    through which the moving indel tunnels. """
    # Move any indels through which this insertion will tunnel.
    for other in indel.pod.indels:
        if (indel.opposite < other.opposite < opposite
                or opposite < other.opposite < indel.opposite):
            other.move(other.opposite, lateral3)
    # Move the indel.
    indel.move(opposite, lateral3)


def _sort_pods(pods: list[IndelPod], move5to3: bool):
    pods.sort(key=(lambda pod: pod.n), reverse=(not move5to3))


def _find_ambindels_recurse(rels: dict[int, int],
                            pods: list[IndelPod],
                            insert3: bool,
                            ref_seq: str,
                            read_seq: str,
                            read_qual: str,
                            min_qual: str,
                            ref_end5: int,
                            ref_end3: int,
                            read_end5: int,
                            read_end3: int,
                            move5to3: bool,
                            pod_index: int,
                            indel_index: int):
    assert 0 <= pod_index < len(pods)
    pod = pods[pod_index]
    assert 0 <= indel_index < len(pod.indels)
    indel = pod.indels[indel_index]
    init_opposite = indel.opposite
    init_lateral3 = indel.lateral3
    if indel.try_move(rels=rels,
                      pods=pods,
                      insert3=insert3,
                      ref_seq=ref_seq,
                      read_seq=read_seq,
                      read_qual=read_qual,
                      min_qual=min_qual,
                      ref_end5=ref_end5,
                      ref_end3=ref_end3,
                      read_end5=read_end5,
                      read_end3=read_end3,
                      move5to3=move5to3):
        # Re-sort the indels in the pod so that they stay in positional
        # order following the move.
        pod.sort(move5to3)
        # Try to move the indel another step.
        _find_ambindels_recurse(rels=rels,
                                pods=pods,
                                insert3=insert3,
                                ref_seq=ref_seq,
                                read_seq=read_seq,
                                read_qual=read_qual,
                                min_qual=min_qual,
                                ref_end5=ref_end5,
                                ref_end3=ref_end3,
                                read_end5=read_end5,
                                read_end3=read_end3,
                                move5to3=move5to3,
                                pod_index=pod_index,
                                indel_index=indel_index)
        if move5to3:
            # Move the indel back to its initial position.
            _move_indels(indel, init_opposite, init_lateral3)
            # Re-sort the indels in the pod so that they stay in
            # positional order following the move.
            pod.sort(move5to3)
    if indel_index + 1 < len(pod.indels):
        # Move the next indel in the pod.
        _find_ambindels_recurse(rels=rels,
                                pods=pods,
                                insert3=insert3,
                                ref_seq=ref_seq,
                                read_seq=read_seq,
                                read_qual=read_qual,
                                min_qual=min_qual,
                                ref_end5=ref_end5,
                                ref_end3=ref_end3,
                                read_end5=read_end5,
                                read_end3=read_end3,
                                move5to3=move5to3,
                                pod_index=pod_index,
                                indel_index=(indel_index + 1))
    if pod_index > 0:
        # Move the indels in the previous pod.
        _find_ambindels_recurse(rels=rels,
                                pods=pods,
                                insert3=insert3,
                                ref_seq=ref_seq,
                                read_seq=read_seq,
                                read_qual=read_qual,
                                min_qual=min_qual,
                                ref_end5=ref_end5,
                                ref_end3=ref_end3,
                                read_end5=read_end5,
                                read_end3=read_end3,
                                move5to3=move5to3,
                                pod_index=(pod_index - 1),
                                indel_index=0)


def find_ambindels(rels: dict[int, int],
                   pods: list[IndelPod],
                   insert3: bool,
                   ref_seq: str,
                   read_seq: str,
                   read_qual: str,
                   min_qual: str,
                   ref_end5: int,
                   ref_end3: int,
                   read_end5: int,
                   read_end3: int):
    if not pods:
        # Nothing to do.
        return
    for move5to3 in [False, True]:
        # Sort the pods.
        _sort_pods(pods, move5to3)
        # Sort the indels in each pod.
        for pod in pods:
            pod.sort(move5to3)
        _find_ambindels_recurse(rels=rels,
                                pods=pods,
                                insert3=insert3,
                                ref_seq=ref_seq,
                                read_seq=read_seq,
                                read_qual=read_qual,
                                min_qual=min_qual,
                                ref_end5=ref_end5,
                                ref_end3=ref_end3,
                                read_end5=read_end5,
                                read_end3=read_end3,
                                move5to3=move5to3,
                                pod_index=(len(pods) - 1),
                                indel_index=0)

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
