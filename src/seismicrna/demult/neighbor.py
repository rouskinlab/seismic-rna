from ..core.seq.xna import DNA


def encode_barcode_2bit(barcode: DNA):
    """
    Encode a DNA barcode (string) into an integer using 2-bit encoding per base.
    Mapping: A=00, C=01, G=10, T=11.
    """
    if not isinstance(barcode, DNA):
        raise ValueError(f"Barcode must be of type DNA but got {type(barcode)}")
    barcode = str(barcode)
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    value = 0
    for base in barcode:
        value = (value << 2) | mapping[base]
    return value


def decode_barcode_2bit(barcode: int, length: int):
    """
    Decode the 2-bit encoded integer back into a DNA string.
    """
    mapping = ["A", "C", "G", "T"]
    bases = []
    for _ in range(length):
        bases.append(mapping[barcode & 0b11])
        barcode >>= 2
    return "".join(reversed(bases))


def _get_neighbors_2bit(barcode: DNA, max_mismatches: int):
    """
    Get neighbors (as DNA strings) using 2-bit encoding (no N substitutions).
    """
    # Imported here (not at module level) so importing this module does not
    # import numba via neighbor_jit.
    from .neighbor_jit import generate_neighbors_2bit

    length = len(barcode)
    original = encode_barcode_2bit(barcode)
    neighbors_int = generate_neighbors_2bit(original, length, max_mismatches)
    return set([decode_barcode_2bit(n, length) for n in neighbors_int])


def encode_barcode_3bit(barcode: DNA):
    """
    Encode a DNA barcode (string) into an integer using 3-bit encoding per base.
    Mapping: A=000, C=001, G=010, T=011.
    (Note: The barcode itself will never contain 'N'.)
    """
    if not isinstance(barcode, DNA):
        raise ValueError(f"Barcode must be of type DNA but got {type(barcode)}")
    barcode = str(barcode)
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    value = 0
    for base in barcode:
        value = (value << 3) | mapping[base]
    return value


def decode_barcode_3bit(barcode: int, length: int):
    """
    Decode the 3-bit encoded integer back into a DNA string.
    Mapping: 0->A, 1->C, 2->G, 3->T, 4->N.
    """
    mapping = ["A", "C", "G", "T", "N"]
    bases = []
    for _ in range(length):
        bases.append(mapping[barcode & 0b111])
        barcode >>= 3
    return "".join(reversed(bases))


def _get_neighbors_3bit(barcode: DNA, max_mismatches: int):
    """
    Get neighbors (as DNA strings) using 3-bit encoding (allowing N substitutions).
    """
    # Imported here (not at module level) so importing this module does not
    # import numba via neighbor_jit.
    from .neighbor_jit import generate_neighbors_3bit

    length = len(barcode)
    original = encode_barcode_3bit(barcode)
    neighbors_int = generate_neighbors_3bit(original, length, max_mismatches)
    return set([decode_barcode_3bit(n, length) for n in neighbors_int])


def get_neighbors(barcode: DNA, max_mismatches: int, allow_n: bool):
    """
    Get all DNA barcodes within max_mismatches.

    Parameters:
      barcode: DNA barcode (of type DNA). It should not contain 'N'.
      max_mismatches: Maximum allowed mismatches.
      allow_n: If True, substitutions to 'N' are allowed (using 3-bit encoding).
               If False, only A, C, G, and T are used (using 2-bit encoding).

    Returns:
      A set of barcode strings representing all neighbors.
    """
    if allow_n:
        return _get_neighbors_3bit(barcode, max_mismatches)
    else:
        return _get_neighbors_2bit(barcode, max_mismatches)
