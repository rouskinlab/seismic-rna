from numba import njit, types
from numba.typed import List

from ..core.seq.xna import DNA

def encode_barcode(barcode: DNA):
    """
    Encode a DNA barcode (string) into an integer using 2-bit encoding per base.
    Mapping: A=00, C=01, G=10, T=11.
    The first base becomes the most significant bits.
    """
    if not isinstance(barcode, DNA):
        raise ValueError(f"Barcode must be of type DNA but got {type(barcode)}")
    barcode = str(barcode)

    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    value = 0
    for base in barcode:
        value = (value << 2) | mapping[base]
    return value

def decode_barcode(barcode: int, length: int):
    """
    Decode the integer x back into a DNA string of given length.
    """
    mapping = ['A', 'C', 'G', 'T']
    bases = []
    for _ in range(length):
        bases.append(mapping[barcode & 0b11])
        barcode >>= 2
    return ''.join(reversed(bases))

@njit(cache=False)
def rec_neighbors(orig: int,
                  length: int,
                  max_mismatches: int,
                  pos: int,
                  mismatches: int,
                  current: int,
                  out: List):
    """
    Recursively build neighbors by working on a binary integer representation.

    original: The original barcode (as an integer).
    length: Number of bases.
    max_mismatches: Maximum allowed mismatches.
    pos: Current position (0-indexed, where 0 is the leftmost base).
    mismatches: Number of mismatches introduced so far.
    current: The neighbor built so far (as an integer).
    out: Numba typed list to store the neighbor integers.
    """
    if pos == length:
        out.append(current)
        return
    # Determine the bit shift for the current position.
    shift = 2 * (length - pos - 1)
    orig_base = (orig >> shift) & 0b11

    # Option 1: Keep the original base.
    candidate = current | (orig_base << shift)
    rec_neighbors(orig, length, max_mismatches, pos + 1, mismatches, candidate, out)

    # Option 2: If allowed, substitute with each alternative base.
    if mismatches < max_mismatches:
        for base in range(4):
            if base != orig_base:
                candidate = current | (base << shift)
                rec_neighbors(orig, length, max_mismatches, pos + 1, mismatches + 1, candidate, out)

@njit(cache=False)
def generate_neighbors(orig: int, length: int, max_mismatches: int):
    """
    Generate all neighbor representations (as integers) using Numba.
    """
    out = List.empty_list(types.int64)  # explicitly specify that we store int64, this caps barcode length at 64.
    rec_neighbors(orig, length, max_mismatches, 0, 0, 0, out)
    return out

def get_neighbors(barcode: DNA, max_mismatches: int):
    """
    Get all DNA barcodes within max_mismatches
    """
    length = len(barcode)
    original = encode_barcode(barcode)
    neighbors_int = generate_neighbors(original, length, max_mismatches)
    return set([decode_barcode(n, length) for n in neighbors_int])
