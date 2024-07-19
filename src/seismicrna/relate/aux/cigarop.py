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


def find_cigar_op_pos(cigar_string: str, find_op: str):
    """ Yield the position in the read of every base with a particular
    type of operation specified by a CIGAR string. """
    if op_consumes_read(find_op):
        # Only these operations correspond to positions in the read.
        pos = 1
        # Check each CIGAR operation, starting at position 1.
        for op, olen in parse_cigar(cigar_string):
            if op == find_op:
                # If the operation matches the code, then yield the position
                # of every base consumed by that operation.
                yield from range(pos, pos := (pos + olen))
            elif op_consumes_read(op):
                # Advance the position by the length of the operation.
                pos += olen


def op_is_mutation(op: str):
    """ Whether the CIGAR operation is a mutation. """
    if op == CIG_SUBST or op == CIG_DELET or op == CIG_INSRT:
        return True
    if op == CIG_ALIGN or op == CIG_MATCH or op == CIG_SCLIP:
        return False
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
