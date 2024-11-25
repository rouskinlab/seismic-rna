from __future__ import annotations

from abc import ABC, abstractmethod
from itertools import chain

from .encode import encode_relate
from ...core.rel import DELET, INS_5, INS_3, MATCH, SUB_N
from ...core.seq import DNA


def get_ins_rel(insert3: bool):
    return INS_3 if insert3 else INS_5


def calc_lateral5(lateral3: int):
    return lateral3 - 1


def get_lateral(lateral3: int, insert3: int):
    return lateral3 if insert3 else calc_lateral5(lateral3)


def _check_out_of_bounds(next_opposite: int, end5: int, end3: int):
    """ Check if the position to move the indel is out of bounds. """
    return next_opposite <= end5 or next_opposite >= end3


def _check_collisions(diffs: list[Indel], next_lateral3: int):
    """ Check if moving the indel will make it collide with another
    indel of the opposite kind. """
    return (_get_indel_by_opp(diffs, next_lateral3) or
            _get_indel_by_opp(diffs, calc_lateral5(next_lateral3)))


class Indel(ABC):
    """ Insertion or deletion. """

    __slots__ = ["opposite", "lateral3"]

    def __init__(self, opposite: int, lateral3: int):
        self.opposite = opposite
        self.lateral3 = lateral3

    @property
    def lateral5(self):
        return calc_lateral5(self.lateral3)

    def get_lateral(self, insert3: bool):
        return get_lateral(self.lateral3, insert3)

    def _get_positions(self, sames: list[Indel], from5to3: bool):
        # Find the position within the indel's own sequence with which
        # to try to swap the indel.
        if from5to3:
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
        while _get_indel_by_opp(sames, next_opposite):
            next_opposite += increment
        return swap_lateral, next_lateral3, next_opposite

    def move(self, opposite: int, lateral3: int):
        self.opposite = opposite
        self.lateral3 = lateral3

    def __str__(self):
        return (f"{type(self).__name__} at {self.opposite} "
                f"{self.lateral5, self.lateral3}")

    def __bool__(self):
        return True

    @abstractmethod
    def _calc_rels(self,
                   ref_seq: DNA,
                   read_seq: DNA,
                   read_qual: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int) -> tuple[int, int]:
        """ Calculate the relationships of the swapping base with the
        current opposite and next opposite bases. """

    @abstractmethod
    def try_move(self,
                 rels: dict[int, int],
                 dels: list[Deletion],
                 inns: list[Insertion],
                 insert3: bool,
                 ref_seq: DNA,
                 read_seq: DNA,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 from5to3: bool) -> bool:
        """ Try to move the indel a step in one direction. """


class Deletion(Indel):

    def _calc_rels(self,
                   ref_seq: DNA,
                   read_seq: DNA,
                   read_qual: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int):
        # Read base with which the deletion will swap.
        read_base = read_seq[swap_lateral - 1]
        read_qual = read_qual[swap_lateral - 1]
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
                 dels: list[Deletion],
                 inns: list[Insertion],
                 insert3: bool,
                 ref_seq: DNA,
                 read_seq: DNA,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 from5to3: bool):
        (swap_lateral,
         next_lateral3,
         next_opposite) = self._get_positions(dels, from5to3)
        if _check_out_of_bounds(next_opposite, ref_end5, ref_end3):
            return False
        if _check_collisions(inns, next_lateral3):
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
        _move_indels(self, dels, next_opposite, next_lateral3)
        # Mark the position to which the deletion moves.
        rels[self.opposite] = rels.get(self.opposite, MATCH) | DELET
        return True


class Insertion(Indel):

    def _calc_rels(self,
                   ref_seq: DNA,
                   read_seq: DNA,
                   read_qual: str,
                   min_qual: str,
                   swap_lateral: int,
                   next_opposite: int):
        # Reference base to whose position the insertion will move.
        ref_base = ref_seq[swap_lateral - 1]
        # Read base that is currently inserted into the reference; will
        # move opposite that reference base.
        read_base_ins = read_seq[self.opposite - 1]
        read_qual_ins = read_qual[self.opposite - 1]
        # Read base that is currently opposite that reference base; will
        # become an insertion after the move.
        read_base_opp = read_seq[next_opposite - 1]
        read_qual_opp = read_qual[next_opposite - 1]
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
                 dels: list[Deletion],
                 inns: list[Insertion],
                 insert3: bool,
                 ref_seq: DNA,
                 read_seq: DNA,
                 read_qual: str,
                 min_qual: str,
                 ref_end5: int,
                 ref_end3: int,
                 read_end5: int,
                 read_end3: int,
                 from5to3: bool):
        (swap_lateral,
         next_lateral3,
         next_opposite) = self._get_positions(inns, from5to3)
        if _check_out_of_bounds(next_opposite, read_end5, read_end3):
            return False
        if _check_collisions(dels, next_lateral3):
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
        if from5to3 == insert3:
            # The insertion moves to the same side as it is marked.
            # Move the indels.
            _move_indels(self, inns, next_opposite, next_lateral3)
            # Check after the move.
            if rel_ins != MATCH or all(ins.get_lateral(insert3) != init_lateral
                                       for ins in inns):
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
             next_opposite_) = self._get_positions(inns, insert3)
            rel_ins_, rel_opp_ = self._calc_rels(ref_seq,
                                                 read_seq,
                                                 read_qual,
                                                 min_qual,
                                                 swap_lateral_,
                                                 next_opposite_)
            # Move the indels.
            _move_indels(self, inns, next_opposite, next_lateral3)
            # Check after the move.
            if rel_opp_ != MATCH or all(ins.get_lateral(insert3) != init_lateral
                                        for ins in inns):
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


def _get_indel_by_opp(indels: list[Indel], opposite: int):
    """ Return the indel that lies opposite the given position (opposite),
    or None if such an indel does not exist. """
    matches = [indel for indel in indels if indel.opposite == opposite]
    if matches:
        if len(matches) > 1:
            raise ValueError(f"{len(matches)} indels have opposite == {opposite}: "
                             f"{list(map(str, matches))}")
        return matches[0]
    return None


def _consistent_rels(rel1: int, rel2: int):
    """ Whether the relationships are consistent. """
    # The relationships are consistent if one of the following is true:
    # - There is at least one type of primary relationship (i.e. match,
    #   substitution to each base) that both of them could both be.
    # - Both could be substitutions to any base.
    return bool((rel1 & rel2) or (rel1 & SUB_N and rel2 & SUB_N))


def _move_indels(indel: Indel, sames: list[Indel], opposite: int, lateral3: int):
    """ Move an indel while adjusting the positions of any other indels
    through which the moving indel tunnels. """
    # Move any indels through which this insertion will tunnel.
    for same in sames:
        if (indel.opposite < same.opposite < opposite
                or opposite < same.opposite < indel.opposite):
            same.move(same.opposite, lateral3)
    # Move the indel.
    indel.move(opposite, lateral3)


def _count_muts(ref_seq: DNA,
                read_seq: DNA,
                read_qual: str,
                min_qual: str,
                ref_end5: int,
                ref_end3: int,
                read_end5: int,
                read_end3: int,
                dels: list[Deletion],
                inns: list[Insertion]):
    """ Count the mutations given the positions of the indels. """
    num_muts = 0
    # Positions in reference and read are initially 0-indexed.
    ref_pos = ref_end5 - 1
    read_pos = read_end5 - 1
    while ref_pos < ref_end3 and read_pos < read_end3:
        # Get the sequence information.
        ref_base = ref_seq[ref_pos]
        read_base = read_seq[read_pos]
        qual = read_qual[read_pos]
        # Make the positions 1-indexed.
        ref_pos += 1
        read_pos += 1
        # Check if this position has an indel.
        has_del = _get_indel_by_opp(dels, ref_pos)
        has_ins = _get_indel_by_opp(inns, read_pos)
        if has_del:
            if has_ins:
                raise ValueError(f"Position {ref_pos, read_pos} has both a "
                                 f"deletion and an insertion")
            # Deletion: move back the read position.
            read_pos -= 1
            num_muts += 1
        elif has_ins:
            # Insertion: move back the reference position.
            ref_pos -= 1
            num_muts += 1
        else:
            # Match or substitution: determine which.
            if qual >= min_qual and read_base != ref_base:
                num_muts += 1
    if ref_pos != ref_end3:
        raise ValueError(
            f"Reached only position {ref_pos} of {ref_end3} in the reference"
        )
    if read_pos != read_end3:
        raise ValueError(
            f"Reached only position {read_pos} of {read_end3} in the read"
        )
    return num_muts


def _find_ambindels(rels: dict[int, int],
                    dels: list[Deletion],
                    inns: list[Insertion],
                    insert3: bool,
                    ref_seq: DNA,
                    read_seq: DNA,
                    read_qual: str,
                    min_qual: str,
                    ref_end5: int,
                    ref_end3: int,
                    read_end5: int,
                    read_end3: int,
                    from5to3: bool,
                    restore: bool):
    for indel in chain(dels, inns):
        init_opp = indel.opposite
        init_lateral3 = indel.lateral3
        if indel.try_move(rels=rels,
                          dels=dels,
                          inns=inns,
                          insert3=insert3,
                          ref_seq=ref_seq,
                          read_seq=read_seq,
                          read_qual=read_qual,
                          min_qual=min_qual,
                          ref_end5=ref_end5,
                          ref_end3=ref_end3,
                          read_end5=read_end5,
                          read_end3=read_end3,
                          from5to3=from5to3):
            # Try to move the other indels if moving this one succeeded.
            _find_ambindels(rels=rels,
                            dels=dels,
                            inns=inns,
                            insert3=insert3,
                            ref_seq=ref_seq,
                            read_seq=read_seq,
                            read_qual=read_qual,
                            min_qual=min_qual,
                            ref_end5=ref_end5,
                            ref_end3=ref_end3,
                            read_end5=read_end5,
                            read_end3=read_end3,
                            from5to3=from5to3,
                            restore=restore)
            if restore:
                if isinstance(indel, Deletion):
                    sames = dels
                elif isinstance(indel, Insertion):
                    sames = inns
                else:
                    raise TypeError(indel)
                # Move the indel(s) back to the original position(s).
                _move_indels(indel, sames, init_opp, init_lateral3)


def find_ambindels(rels: dict[int, int],
                   dels: list[Deletion],
                   inns: list[Insertion],
                   insert3: bool,
                   ref_seq: DNA,
                   read_seq: DNA,
                   read_qual: str,
                   min_qual: str,
                   ref_end5: int,
                   ref_end3: int,
                   read_end5: int,
                   read_end3: int):
    for mode in [False, True]:
        _find_ambindels(rels=rels,
                        dels=dels,
                        inns=inns,
                        insert3=insert3,
                        ref_seq=ref_seq,
                        read_seq=read_seq,
                        read_qual=read_qual,
                        min_qual=min_qual,
                        ref_end5=ref_end5,
                        ref_end3=ref_end3,
                        read_end5=read_end5,
                        read_end3=read_end3,
                        from5to3=mode,
                        restore=mode)
