"""

Relate Core Module

========================================================================

Convert the relationships between reads and a reference from SAM format
(which encodes relationships implicitly as CIGAR strings) to vectorized
format (which encodes relationships explicitly as elements of arrays).

------------------------------------------------------------------------

"""

import re
from sys import byteorder

import numpy as np
import pandas as pd

from .cli import opt_phred_enc
from .sect import seq_pos_to_index
from .seq import BASEA, BASEC, BASEG, BASET, BASEN, DNA

# Data type for NumPy and NumPy-like arrays of relation vectors.
NP_TYPE = np.uint8

# Minimum number of non-indel bases between an insertion and a deletion.
MIN_INDEL_GAP = 1

# Integer encodings for relation vectors
IRREC = b"\x00"[0]  # 00000000 (000): irreconcilable paired mates
MATCH = b"\x01"[0]  # 00000001 (001): match with reference
DELET = b"\x02"[0]  # 00000010 (002): deletion from reference
INS_5 = b"\x04"[0]  # 00000100 (004): insertion 5' of base in reference
INS_3 = b"\x08"[0]  # 00001000 (008): insertion 3' of base in reference
SUB_A = b"\x10"[0]  # 00010000 (016): substitution to A
SUB_C = b"\x20"[0]  # 00100000 (032): substitution to C
SUB_G = b"\x40"[0]  # 01000000 (064): substitution to G
SUB_T = b"\x80"[0]  # 10000000 (128): substitution to T
NOCOV = b"\xff"[0]  # 11111111 (255): not covered by read
SUB_N = SUB_A | SUB_C | SUB_G | SUB_T
SUB_B = SUB_N ^ SUB_A
SUB_D = SUB_N ^ SUB_C
SUB_H = SUB_N ^ SUB_G
SUB_V = SUB_N ^ SUB_T
ANY_N = SUB_N | MATCH
ANY_B = SUB_B | MATCH
ANY_D = SUB_D | MATCH
ANY_H = SUB_H | MATCH
ANY_V = SUB_V | MATCH
INS_8 = INS_5 | INS_3
INDEL = DELET | INS_8
MINS5 = INS_5 | MATCH
MINS3 = INS_3 | MATCH
ANY_8 = INS_8 | MATCH

# Map bases to integer encodings and vice versa.
BASE_ENCODINGS = {BASEA: SUB_A, BASEC: SUB_C, BASEG: SUB_G, BASET: SUB_T,
                  BASEN: IRREC}
BASE_DECODINGS = {code: base for base, code in BASE_ENCODINGS.items()
                  if base != BASEN}

# Number of unique bytes
N_BYTES = 256  # = 2^8
# Map each common byte in the vector encoding (e.g. MATCH) to a byte of
# human-readable text.
BYTE_SYMBOLS = {NOCOV: b'_',
                MATCH: b'=',
                DELET: b'.',
                INS_5 | MATCH: b'{',
                INS_5 | MATCH | INS_3: b'+',
                INS_3 | MATCH: b'}',
                SUB_A: b'A',
                SUB_C: b'C',
                SUB_G: b'G',
                SUB_T: b'T',
                SUB_A | SUB_C | SUB_G | MATCH: b'N',
                SUB_T | SUB_A | SUB_C | MATCH: b'N',
                SUB_G | SUB_T | SUB_A | MATCH: b'N',
                SUB_C | SUB_G | SUB_T | MATCH: b'N',
                IRREC: b'!'}
# Default for uncommon bytes in the vector encoding.
OTHER_SYMBOL = b'?'
# Create an array that maps each vector byte to its readable character.
map_array = bytearray(OTHER_SYMBOL) * N_BYTES
for byte, symbol in BYTE_SYMBOLS.items():
    map_array[byte] = int.from_bytes(symbol, byteorder)
# Create a translation table from vector to human-readable encodings.
map_table = bytes.maketrans(bytes(range(N_BYTES)), map_array)

# CIGAR string operation codes
CIG_ALIGN = 'M'  # alignment match
CIG_MATCH = '='  # sequence match
CIG_SUBST = 'X'  # substitution
CIG_DELET = 'D'  # deletion
CIG_INSRT = 'I'  # insertion
CIG_SCLIP = 'S'  # soft clipping

# Regular expression pattern that matches a single CIGAR operation
# (length ≥ 1 and operation code, defined above)
CIG_PATTERN = re.compile("".join([r"(\d+)([",
                                  CIG_ALIGN,
                                  CIG_MATCH,
                                  CIG_SUBST,
                                  CIG_DELET,
                                  CIG_INSRT,
                                  CIG_SCLIP,
                                  "])"]))

# Define default values for low-, medium-, and high-quality base calls.
MIN_QUAL = chr(opt_phred_enc.default)
MED_QUAL = chr(opt_phred_enc.default + 20)
MAX_QUAL = chr(opt_phred_enc.default + 40)


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
        raise ValueError("CIGAR string is empty")
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
    for match in CIG_PATTERN.finditer(cigar_string):
        length_str, operation = match.groups()
        # Convert the length field from str to int and verify that it
        # is a positive integer.
        if (length_int := int(length_str)) < 1:
            raise ValueError("length of CIGAR operation must be ≥ 1")
        # Add the total number of characters in the current operation to
        # the total number of characters matched from the CIGAR string.
        num_chars_matched += len(length_str) + len(operation)
        # Note that the fields are yielded as (operation, length), but
        # in the CIGAR string itself, the order is (length, operation).
        yield operation, length_int
    # Confirm that all bytes in the CIGAR string were matched by the
    # regular expression.
    if num_chars_matched != len(cigar_string):
        raise ValueError(f"Invalid CIGAR: '{cigar_string}'")


def translate_relvec(relvec: np.ndarray):
    """ Translate a binary relation vector to human-readable text. """
    return relvec.tobytes().translate(map_table)


def blank_relvec(bases: int | DNA,
                 reads: int | list | np.ndarray | pd.Index | None = None):
    """
    Return blank relation vector(s) of a given length.

    Parameters
    ----------
    bases: int | DNA
        The reference sequence (if DNA) or just its length (if int).
        If the sequence itself is given, then return a Series/DataFrame
        whose main index (index for Series, columns for DataFrame) is a
        MultiIndex of 1-indexed positions and bases in the sequence.
        If only the sequence length is given, then return a NumPy array.
    reads: int | list | np.ndarray | pd.Index | None = None
        Optional names of the relation vectors. If given, then return a
        DataFrame whose indexes are the read names if bases is DNA,
        otherwise a 2D NumPy array with one row for each name in reads.
        If None, then return a Series if bases is DNA, otherwise a 1D
        NumPy array.

    Returns
    -------
    np.ndarray | pd.Series | pd.DataFrame
        The blank relation vector(s).
    """
    # Determine whether to return a Pandas or NumPy object.
    if isinstance(bases, DNA):
        # Make a Series or DataFrame with the sequence as its main axis.
        sequence = seq_pos_to_index(bases, np.arange(1, len(bases) + 1), 1)
        if reads is None:
            # Return a Series representing just one relation vector.
            return pd.Series(NOCOV, index=sequence, dtype=NP_TYPE)
        # Ensure that names is a sequence of read names as str objects.
        if isinstance(reads, int):
            names = [f"Read_{i}" for i in range(1, reads + 1)]
        else:
            names = list(map(str, reads))
        # Return a DataFrame with one row per relation vector.
        return pd.DataFrame(NOCOV, index=names, columns=sequence, dtype=NP_TYPE)
    # Determine the size of the NumPy array.
    size = (bases if reads is None
            else ((reads, bases) if isinstance(reads, int)
                  else (len(reads), bases)))
    return np.full(size, fill_value=NOCOV, dtype=NP_TYPE)


def encode_relate(ref_base: str, read_base: str, read_qual: str, min_qual: str):
    """
    Encode the relation between a base in the read and a base in the
    reference sequence. If the read quality is sufficient, then return
    the match encoding if the read and reference bases match, otherwise
    the encoding of the substitution for the base in the read.
    If the read quality is insufficient, then return the fully ambiguous
    base encoding, that is a match or substitution to any base except
    the reference base, since a "substitution to the reference base"
    would be a match, not a substitution.

    Parameters
    ----------
    ref_base: DNA
        Base in the reference sequence.
    read_base: DNA
        Base in the read sequence.
    read_qual: str
        ASCII encoding for the Phred quality score of the read base.
    min_qual: str
        Minimum value of `read_qual` to not call the relation ambiguous.
    """
    return ((MATCH if ref_base == read_base
             else BASE_ENCODINGS[read_base]) if (read_qual >= min_qual
                                                 and read_base != BASEN
                                                 and ref_base != BASEN)
            else ANY_N ^ BASE_ENCODINGS[ref_base])


def encode_match(read_base: str, read_qual: str, min_qual: str):
    """
    A more efficient version of `encode_compare` given prior knowledge
    from the CIGAR string that the read and reference match at this
    position. Note that there is no analagous version when there is a
    known substitution because substitutions are relatively infrequent,
    so optimizing their processing would speed the program only slightly
    while making the source code more complex and harder to maintain.
    """
    return (MATCH if (read_qual >= min_qual and read_base != BASEN)
            else ANY_N ^ BASE_ENCODINGS[read_base])
