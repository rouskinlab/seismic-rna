from .cigar import (CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                    CIG_DELET, CIG_INSRT, CIG_SCLIP,
                    parse_cigar)


class CigarOp(object):
    """ Represent one operation in a CIGAR string. """

    def __init__(self, op: str):
        if op not in (CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                      CIG_DELET, CIG_INSRT, CIG_SCLIP):
            raise ValueError(f"Invalid CIGAR operation: '{op}'")
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


def count_cigar_muts(cigar_string: str):
    """ Return the total number of mutations in a CIGAR string. """
    mutation_types = CIG_SUBST, CIG_DELET, CIG_INSRT
    return sum(olen for op, olen in parse_cigar(cigar_string)
               if op in mutation_types)


def find_cigar_op_pos(cigar_string: str, find_op: str):
    """ Yield the position in the read of every base with a particular
    type of operation specified by a CIGAR string. """
    consume_read = CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_INSRT
    if find_op in consume_read:
        # Only these operations correspond to positions in the read.
        pos = 1
        # Check each CIGAR operation, starting at position 1.
        for op, olen in parse_cigar(cigar_string):
            if op == find_op:
                # If the operation matches the code, then yield the position
                # of every base consumed by that operation.
                yield from range(pos, pos := (pos + olen))
            elif op in consume_read:
                # Advance the position by the length of the operation.
                pos += olen
