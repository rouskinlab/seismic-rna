import re

from .error import RelateError

# CIGAR string operation codes
CIG_ALIGN = "M"  # alignment match
CIG_MATCH = "="  # sequence match
CIG_SUBST = "X"  # substitution
CIG_DELET = "D"  # deletion
CIG_INSRT = "I"  # insertion
CIG_SCLIP = "S"  # soft clipping
CIG_INTRN = "N"  # soft clipping

# Regular expression pattern that matches a single CIGAR operation
# (length ≥ 1 and operation code, defined above)
CIG_PATTERN = re.compile("".join([r"(\d+)([",
                                  CIG_ALIGN,
                                  CIG_MATCH,
                                  CIG_SUBST,
                                  CIG_DELET,
                                  CIG_INSRT,
                                  CIG_SCLIP,
                                  CIG_INTRN,
                                  "])"]))


def parse_cigar(cigar_string: str):
    """
    Yield the fields of a CIGAR string as pairs of (operation, length),
    where operation is 1 byte indicating the CIGAR operation and length
    is a positive integer indicating the number of bases from the read
    that the operation consumes. Note that in the CIGAR string itself,
    each length precedes its corresponding operation.

    Parameters
    ----------
    cigar_string: bytes
        CIGAR string from a SAM file. For full documentation, refer to
        https://samtools.github.io/hts-specs/
    Yield
    -----
    bytes (length = 1)
        Current CIGAR operation
    int (≥ 1)
        Length of current CIGAR operation
    """
    # Length-0 CIGAR strings are forbidden.
    if not cigar_string:
        raise RelateError("CIGAR string is empty")
    # If the CIGAR string has any invalid bytes (e.g. an unrecognized
    # operation byte, an operation longer than 1 byte, a length that is
    # not a positive integer, or any extraneous characters), then the
    # regular expression parser will simply skip these invalid bytes.
    # In order to catch such problems, keep track of the number of
    # bytes matched from the CIGAR string. After reading the CIGAR, if
    # the number of bytes matched is smaller than the length of the
    # CIGAR string, then some bytes must have been skipped, indicating
    # that the CIGAR string contained at least one invalid byte.
    num_chars_matched = 0
    # Find every operation in the CIGAR string that matches the regular
    # expression.
    prev_operation = ""
    for match in CIG_PATTERN.finditer(cigar_string):
        length_str, operation = match.groups()
        # Consecutive CIGAR operations cannot be of the same type.
        if prev_operation:
            if operation == prev_operation:
                raise RelateError("Identical consecutive CIGAR operations")
            if {operation, prev_operation} == {CIG_DELET, CIG_INSRT}:
                raise RelateError("Adjacent insertion and deletion")
            if {operation, prev_operation} == {CIG_DELET, CIG_INTRN}:
                raise RelateError("Adjacent intron and deletion")
        # Convert the length field from str to int and verify that it
        # is a positive integer.
        if (length_int := int(length_str)) < 1:
            raise RelateError("CIGAR operation has length 0")
        # Add the total number of characters in the current operation to
        # the total number of characters matched from the CIGAR string.
        num_chars_matched += len(length_str) + len(operation)
        # Note that the fields are yielded as (operation, length), but
        # in the CIGAR string itself, the order is (length, operation).
        yield operation, length_int
        prev_operation = operation
    # Confirm that all bytes in the CIGAR string were matched by the
    # regular expression.
    if num_chars_matched != len(cigar_string):
        raise RelateError("Invalid CIGAR string")


def op_consumes_ref(op: str):
    """ Whether the CIGAR operation consumes the reference. """
    if op == CIG_ALIGN or op == CIG_MATCH or op == CIG_SUBST or op == CIG_DELET:
        return True
    if op == CIG_INSRT or op == CIG_SCLIP:
        return False
    raise RelateError("Unsupported CIGAR operation")


def op_consumes_read(op: str):
    """ Whether the CIGAR operation consumes the read. """
    if op == CIG_DELET:
        return False
    if (op == CIG_ALIGN
            or op == CIG_MATCH
            or op == CIG_SUBST
            or op == CIG_INSRT
            or op == CIG_SCLIP):
        return True
    raise RelateError("Unsupported CIGAR operation")
