import struct
import tempfile
import unittest as ut
from pathlib import Path

from seismicrna.core.seq import DNA
from seismicrna.importmm.mm import iter_mm_file

# ---------------------------------------------------------------------------
# Shared MM binary fixture
# ---------------------------------------------------------------------------
# Nucleotide encoding used by RNA Framework: A=0, C=1, G=2, T=3, N=4.
# Two nucleotides are packed into one byte as hex nibbles (high nibble first),
# matching Perl's pack("H*", ...) behaviour.
#
# Ref 1: "ACGT" → hex "0123" → b'\x01\x23'
#   3 reads:
#     read 0: start=0, end=3, 1 mutation at pos 1 (0-based; the C base)
#     read 1: start=0, end=3, 2 mutations at pos 0 and 3 (A and T bases)
#     read 2: start=1, end=2, 0 mutations
#
# Ref 2: "TGCA" → hex "3210" → b'\x32\x10'
#   2 reads:
#     read 0: start=0, end=3, 1 mutation at pos 2 (0-based; the C base)
#     read 1: start=0, end=3, 0 mutations

_REF1 = "ref1"
_REF1_SEQ = DNA("ACGT")
_REF2 = "ref2"
_REF2_SEQ = DNA("TGCA")


def _build_mm_bytes() -> bytes:
    """ Build a two-transcript MM binary for use in tests. """

    def _transcript(ref_id: str,
                    seq_packed: bytes,
                    seq_len: int,
                    reads: list[tuple[int, int, list[int]]]) -> bytes:
        id_bytes = ref_id.encode() + b"\x00"
        buf = struct.pack("<H", len(id_bytes))
        buf += id_bytes
        buf += struct.pack("<L", seq_len)
        buf += seq_packed
        buf += struct.pack("<L", len(reads))
        for start, end, muts in reads:
            buf += struct.pack("<LLL", start, end, len(muts))
            if muts:
                buf += struct.pack(f"<{len(muts)}L", *muts)
        return buf

    buf = _transcript(_REF1, bytes.fromhex("0123"), 4, [
        (0, 3, [1]),
        (0, 3, [0, 3]),
        (1, 2, []),
    ])
    buf += _transcript(_REF2, bytes.fromhex("3210"), 4, [
        (0, 3, [2]),
        (0, 3, []),
    ])
    buf += b"[mmeof]"
    return buf


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestIterMMFile(ut.TestCase):

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self._mm_path = Path(self._tmpdir.name) / "test.mm"
        self._mm_path.write_bytes(_build_mm_bytes())
        # Parse once; individual tests pull from this list.
        self._transcripts = list(iter_mm_file(self._mm_path))

    def tearDown(self):
        self._tmpdir.cleanup()

    # --- transcript count and identity ---

    def test_yields_two_transcripts(self):
        self.assertEqual(len(self._transcripts), 2)

    def test_ref1_id(self):
        ref_id, _, _ = self._transcripts[0]
        self.assertEqual(ref_id, _REF1)

    def test_ref2_id(self):
        ref_id, _, _ = self._transcripts[1]
        self.assertEqual(ref_id, _REF2)

    # --- sequence decoding ---

    def test_ref1_sequence(self):
        _, refseq, _ = self._transcripts[0]
        self.assertEqual(refseq, _REF1_SEQ)

    def test_ref2_sequence(self):
        _, refseq, _ = self._transcripts[1]
        self.assertEqual(refseq, _REF2_SEQ)

    # --- read counts ---

    def test_ref1_read_count(self):
        _, _, reads = self._transcripts[0]
        self.assertEqual(len(reads), 3)

    def test_ref2_read_count(self):
        _, _, reads = self._transcripts[1]
        self.assertEqual(len(reads), 2)

    # --- ref1 read coordinates and mutations ---

    def test_ref1_read0_start_end(self):
        _, _, reads = self._transcripts[0]
        start, end, _ = reads[0]
        self.assertEqual(start, 0)
        self.assertEqual(end, 3)

    def test_ref1_read0_single_mutation(self):
        _, _, reads = self._transcripts[0]
        _, _, muts = reads[0]
        self.assertEqual(muts, [1])

    def test_ref1_read1_two_mutations(self):
        _, _, reads = self._transcripts[0]
        _, _, muts = reads[1]
        self.assertEqual(muts, [0, 3])

    def test_ref1_read2_no_mutations(self):
        _, _, reads = self._transcripts[0]
        start, end, muts = reads[2]
        self.assertEqual(start, 1)
        self.assertEqual(end, 2)
        self.assertEqual(muts, [])

    # --- ref2 read coordinates and mutations ---

    def test_ref2_read0_mutation(self):
        _, _, reads = self._transcripts[1]
        start, end, muts = reads[0]
        self.assertEqual(start, 0)
        self.assertEqual(end, 3)
        self.assertEqual(muts, [2])

    def test_ref2_read1_no_mutations(self):
        _, _, reads = self._transcripts[1]
        _, _, muts = reads[1]
        self.assertEqual(muts, [])


if __name__ == "__main__":
    ut.main(verbosity=2)
