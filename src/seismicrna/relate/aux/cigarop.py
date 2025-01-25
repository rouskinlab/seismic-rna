from typing import Callable

from ..py.cigar import (CIG_ALIGN,
                        CIG_MATCH,
                        CIG_SUBST,
                        CIG_DELET,
                        CIG_INSRT,
                        CIG_SCLIP,
                        parse_cigar,
                        op_consumes_read,
                        op_consumes_ref)


class CigarOp(object):
    """ Represent one operation in a CIGAR string. """

    def __init__(self, op: str):
        if op not in (CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                      CIG_DELET, CIG_INSRT, CIG_SCLIP):
            raise ValueError(f"Invalid CIGAR operation: {repr(op)}")
        self._op = op
        self._len = 1

    @property
    def op(self):
        """ CIGAR operation as a character. """
        return self._op

    def lengthen(self):
        """ Lengthen the operation by 1 base call. """
        self._len += 1

    def __str__(self):
        """ Text that goes into the CIGAR string. """
        return f"{self._len}{self._op}"


def _find_cigar_op_pos(cigar_string: str,
                       find_op: str,
                       init_pos: int,
                       op_consumes_func: Callable[[str], bool]):
    if op_consumes_func(find_op):
        pos = init_pos
        # Check each CIGAR operation, starting at position 1.
        for op, olen in parse_cigar(cigar_string):
            if op == find_op:
                # If the operation matches the code, then yield the position
                # of every base consumed by that operation.
                yield from range(pos, pos := (pos + olen))
            elif op_consumes_func(op):
                # Advance the position by the length of the operation.
                pos += olen


def find_cigar_op_pos_read(cigar_string: str, find_op: str):
    """ Yield the position in the read of every base with a type of
    operation specified by a CIGAR string. """
    yield from _find_cigar_op_pos(cigar_string, find_op, 1, op_consumes_read)


def find_cigar_op_pos_ref(cigar_string: str, find_op: str, end5: int):
    """ Yield the position in the reference of every base with a type of
    operation specified by a CIGAR string. """
    yield from _find_cigar_op_pos(cigar_string, find_op, end5, op_consumes_ref)


def op_is_mutation(op: str):
    """ Whether the CIGAR operation is a mutation. """
    if op == CIG_ALIGN or op == CIG_MATCH or op == CIG_SCLIP:
        return False
    if op == CIG_SUBST or op == CIG_DELET or op == CIG_INSRT:
        return True
    raise ValueError(f"Invalid CIGAR operation: {repr(op)}")


def count_cigar_read(cigar_string: str):
    """ Count the read positions consumed by a CIGAR string. """
    return sum(olen for op, olen in parse_cigar(cigar_string)
               if op_consumes_read(op))


def count_cigar_ref(cigar_string: str):
    """ Count the reference positions consumed by a CIGAR string. """
    return sum(olen for op, olen in parse_cigar(cigar_string)
               if op_consumes_ref(op))


def count_cigar_muts(cigar_string: str):
    """ Count the mutations in a CIGAR string. """
    return sum(olen for op, olen in parse_cigar(cigar_string)
               if op_is_mutation(op))
