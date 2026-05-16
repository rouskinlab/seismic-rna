import struct
from pathlib import Path
from typing import Iterator

from ..core.seq import DNA

EOF_MARKER = b"[mmeof]"

# Hex nibble digit → nucleotide (RNA Framework encoding: 4=N, 3=T, 2=G, 1=C, 0=A)
_HEX_TO_NUC = str.maketrans("43210", "NTGCA")


def _read_exact(fh, n: int) -> bytes:
    data = fh.read(n)
    if len(data) != n:
        raise EOFError(f"Expected {n} bytes but got {len(data)}")
    return data


def _parse_sequence(packed: bytes, length: int) -> DNA:
    """ Decode 4-bit packed sequence (2 nucleotides per byte, hex nibbles). """
    hex_str = packed.hex().translate(_HEX_TO_NUC)[:length]
    return DNA(hex_str)


def _read_transcript(fh) -> tuple[str, DNA, list[tuple[int, int, list[int]]]]:
    """ Read one transcript block from an open MM file handle.

    Returns (ref_id, refseq, reads) where each read is
    (start, end, mut_positions) with 0-based coordinates.
    """
    # Read ID length (uint16, little-endian).
    id_len = struct.unpack("<H", _read_exact(fh, 2))[0]
    # Read null-terminated ID string.
    ref_id = _read_exact(fh, id_len).rstrip(b"\x00").decode()
    # Read sequence length (uint32).
    seq_len = struct.unpack("<L", _read_exact(fh, 4))[0]
    # Read packed sequence (2 nucleotides per byte, ceil(n/2) bytes).
    packed_len = (seq_len + 1) // 2
    packed = _read_exact(fh, packed_len)
    refseq = _parse_sequence(packed, seq_len)
    # Read number of reads (uint32).
    read_count = struct.unpack("<L", _read_exact(fh, 4))[0]
    # Read all per-read blocks eagerly so the file handle advances fully.
    reads: list[tuple[int, int, list[int]]] = []
    for _ in range(read_count):
        start, end, n_muts = struct.unpack("<LLL", _read_exact(fh, 12))
        if n_muts:
            mut_positions = list(
                struct.unpack(f"<{n_muts}L", _read_exact(fh, 4 * n_muts))
            )
        else:
            mut_positions = []
        reads.append((start, end, mut_positions))
    return ref_id, refseq, reads


def iter_mm_file(
        mm_path: Path
) -> Iterator[tuple[str, DNA, list[tuple[int, int, list[int]]]]]:
    """ Iterate over transcript blocks in an MM file.

    Yields (ref_id, refseq, reads) for each transcript, where reads is a
    list of (start, end, mut_positions) tuples (all positions 0-based).
    """
    with open(mm_path, "rb") as fh:
        while True:
            # Peek at the next bytes to detect the EOF marker or file end.
            peek = fh.read(len(EOF_MARKER))
            if not peek:
                # Reached end of file (no EOF marker — appendable MM files
                # may omit it).
                break
            if peek == EOF_MARKER:
                break
            # Not the EOF marker: rewind and parse as a transcript header.
            fh.seek(-len(peek), 1)
            yield _read_transcript(fh)
