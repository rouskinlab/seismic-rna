from numba import njit, types
from numba.typed import List

from ..core.seq.xna import DNA

# --------------------------
# 2-bit encoding functions
# (for when substitutions to N are NOT allowed)
# --------------------------
def encode_barcode_2bit(barcode: DNA):
    """
    Encode a DNA barcode (string) into an integer using 2-bit encoding per base.
    Mapping: A=00, C=01, G=10, T=11.
    """
    if not isinstance(barcode, DNA):
        raise ValueError(f"Barcode must be of type DNA but got {type(barcode)}")
    barcode = str(barcode)
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    value = 0
    for base in barcode:
        value = (value << 2) | mapping[base]
    return value

def decode_barcode_2bit(barcode: int, length: int):
    """
    Decode the 2-bit encoded integer back into a DNA string.
    """
    mapping = ['A', 'C', 'G', 'T']
    bases = []
    for _ in range(length):
        bases.append(mapping[barcode & 0b11])
        barcode >>= 2
    return ''.join(reversed(bases))

@njit(cache=False)
def rec_neighbors_2bit(orig: int,
                       length: int,
                       max_mismatches: int,
                       pos: int,
                       mismatches: int,
                       current: int,
                       out: List):
    """
    Recursively generate neighbor integers using 2-bit encoding.
    """
    if pos == length:
        out.append(current)
        return
    shift = 2 * (length - pos - 1)
    orig_base = (orig >> shift) & 0b11
    # Option 1: Keep the original base.
    candidate = current | (orig_base << shift)
    rec_neighbors_2bit(orig, length, max_mismatches, pos + 1, mismatches, candidate, out)
    # Option 2: Substitute with any alternative (A, C, G, T).
    if mismatches < max_mismatches:
        for base in range(4):
            if base != orig_base:
                candidate = current | (base << shift)
                rec_neighbors_2bit(orig, length, max_mismatches, pos + 1, mismatches + 1, candidate, out)

@njit(cache=False)
def generate_neighbors_2bit(orig: int, length: int, max_mismatches: int):
    """
    Generate all neighbor integers using 2-bit encoding.
    """
    out = List.empty_list(types.int64)
    rec_neighbors_2bit(orig, length, max_mismatches, 0, 0, 0, out)
    return out

def _get_neighbors_2bit(barcode: DNA, max_mismatches: int):
    """
    Get neighbors (as DNA strings) using 2-bit encoding (no N substitutions).
    """
    length = len(barcode)
    original = encode_barcode_2bit(barcode)
    neighbors_int = generate_neighbors_2bit(original, length, max_mismatches)
    return set([decode_barcode_2bit(n, length) for n in neighbors_int])

# --------------------------
# 3-bit encoding functions
# (for when substitutions to N are allowed)
# --------------------------
def encode_barcode_3bit(barcode: DNA):
    """
    Encode a DNA barcode (string) into an integer using 3-bit encoding per base.
    Mapping: A=000, C=001, G=010, T=011.
    (Note: The barcode itself will never contain 'N'.)
    """
    if not isinstance(barcode, DNA):
        raise ValueError(f"Barcode must be of type DNA but got {type(barcode)}")
    barcode = str(barcode)
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    value = 0
    for base in barcode:
        value = (value << 3) | mapping[base]
    return value

def decode_barcode_3bit(barcode: int, length: int):
    """
    Decode the 3-bit encoded integer back into a DNA string.
    Mapping: 0->A, 1->C, 2->G, 3->T, 4->N.
    """
    mapping = ['A', 'C', 'G', 'T', 'N']
    bases = []
    for _ in range(length):
        bases.append(mapping[barcode & 0b111])
        barcode >>= 3
    return ''.join(reversed(bases))

@njit(cache=False)
def rec_neighbors_3bit(orig: int,
                       length: int,
                       max_mismatches: int,
                       pos: int,
                       mismatches: int,
                       current: int,
                       out: List):
    """
    Recursively generate neighbor integers using 3-bit encoding.
    This version allows substitutions to 'N' (encoded as 4).
    """
    if pos == length:
        out.append(current)
        return
    shift = 3 * (length - pos - 1)
    orig_base = (orig >> shift) & 0b111  # Will be one of 0-3 (barcode has no 'N').
    # Option 1: Keep the original base.
    candidate = current | (orig_base << shift)
    rec_neighbors_3bit(orig, length, max_mismatches, pos + 1, mismatches, candidate, out)
    # Option 2: Substitute with any alternative base (A, C, G, T, or N).
    if mismatches < max_mismatches:
        for base in range(5):
            if base != orig_base:
                candidate = current | (base << shift)
                rec_neighbors_3bit(orig, length, max_mismatches, pos + 1, mismatches + 1, candidate, out)

@njit(cache=False)
def generate_neighbors_3bit(orig: int, length: int, max_mismatches: int):
    """
    Generate all neighbor integers using 3-bit encoding.
    """
    out = List.empty_list(types.int64)
    rec_neighbors_3bit(orig, length, max_mismatches, 0, 0, 0, out)
    return out

def _get_neighbors_3bit(barcode: DNA, max_mismatches: int):
    """
    Get neighbors (as DNA strings) using 3-bit encoding (allowing N substitutions).
    """
    length = len(barcode)
    original = encode_barcode_3bit(barcode)
    neighbors_int = generate_neighbors_3bit(original, length, max_mismatches)
    return set([decode_barcode_3bit(n, length) for n in neighbors_int])

# --------------------------
# Public API: get_neighbors
# --------------------------
def get_neighbors(barcode: DNA, max_mismatches: int, allow_n: bool):
    """
    Get all DNA barcodes within max_mismatches.

    Parameters:
      barcode: DNA barcode (of type DNA). It should not contain 'N'.
      max_mismatches: Maximum allowed mismatches.
      allow_n_substitution: If True, substitutions to 'N' are allowed (using 3-bit encoding).
                            If False, only A, C, G, and T are used (using 2-bit encoding).

    Returns:
      A set of barcode strings representing all neighbors.
    """
    if allow_n:
        return _get_neighbors_3bit(barcode, max_mismatches)
    else:
        return _get_neighbors_2bit(barcode, max_mismatches)
