import json
import struct
import tempfile
import unittest as ut
from pathlib import Path

from seismicrna.core import path
from seismicrna.core.logs import Level, set_config
from seismicrna.core.rel.code import DELET, INS_3, INS_5, SUB_A, SUB_C, SUB_G, SUB_N, SUB_T
from seismicrna.core.report import NumReadsXamF, NumReadsRelF, RefF, SampleF
from seismicrna.core.seq import DNA
from seismicrna.importmm.write import _build_mut_codes, import_mm
from seismicrna.relate.report import RelateReport

# ---------------------------------------------------------------------------
# Shared MM binary fixture  (same layout as mm_test.py)
# ---------------------------------------------------------------------------
# Ref 1: "ACGT" (4 bp), 3 reads
# Ref 2: "TGCA" (4 bp), 2 reads

_SAMPLE = "test_sample"
_BRANCH = ""
_REF1 = "ref1"
_REF1_SEQ = DNA("ACGT")
_REF2 = "ref2"
_REF2_SEQ = DNA("TGCA")


def _build_mm_bytes() -> bytes:
    """ Build the two-transcript MM binary used in all write tests. """

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


def _report_path(out_dir: Path, ref: str) -> Path:
    """ Return the expected path of a RelateReport JSON file. """
    branches = path.add_branch(path.RELATE_STEP, _BRANCH, {})
    return RelateReport.build_path(
        {path.TOP: out_dir,
         path.SAMPLE: _SAMPLE,
         path.BRANCHES: branches,
         path.REF: ref}
    )


# ---------------------------------------------------------------------------
# Unit tests for _build_mut_codes
# ---------------------------------------------------------------------------

class TestBuildMutCodes(ut.TestCase):
    """ Verify that per-position mutation codes are computed correctly. """

    def _codes(self, seq: str, insert3: bool) -> list[int]:
        """ Return 1-indexed mutation codes for every position. """
        codes = _build_mut_codes(DNA(seq), insert3)
        return [int(codes[i]) for i in range(1, len(seq) + 1)]

    def test_insert3_true_base_a(self):
        # A: all subs except SUB_A, plus DELET and INS_3
        expected = DELET | INS_3 | SUB_C | SUB_G | SUB_T
        self.assertEqual(self._codes("A", insert3=True)[0], expected)

    def test_insert3_true_base_c(self):
        expected = DELET | INS_3 | SUB_A | SUB_G | SUB_T
        self.assertEqual(self._codes("C", insert3=True)[0], expected)

    def test_insert3_true_base_g(self):
        expected = DELET | INS_3 | SUB_A | SUB_C | SUB_T
        self.assertEqual(self._codes("G", insert3=True)[0], expected)

    def test_insert3_true_base_t(self):
        expected = DELET | INS_3 | SUB_A | SUB_C | SUB_G
        self.assertEqual(self._codes("T", insert3=True)[0], expected)

    def test_insert3_true_base_n(self):
        # N: all substitutions allowed (SUB_N = SUB_A|SUB_C|SUB_G|SUB_T)
        expected = DELET | INS_3 | SUB_N
        self.assertEqual(self._codes("N", insert3=True)[0], expected)

    def test_insert3_false_uses_ins5(self):
        # With insert3=False the INS_5 bit should be set, not INS_3.
        codes_true = self._codes("A", insert3=True)
        codes_false = self._codes("A", insert3=False)
        self.assertNotEqual(codes_true[0] & INS_3, 0)   # insert3=True → INS_3 set
        self.assertEqual(codes_false[0] & INS_3, 0)     # insert3=False → INS_3 clear
        self.assertNotEqual(codes_false[0] & INS_5, 0)  # insert3=False → INS_5 set

    def test_acgt_all_positions(self):
        codes = self._codes("ACGT", insert3=True)
        self.assertEqual(codes[0], DELET | INS_3 | SUB_C | SUB_G | SUB_T)  # A
        self.assertEqual(codes[1], DELET | INS_3 | SUB_A | SUB_G | SUB_T)  # C
        self.assertEqual(codes[2], DELET | INS_3 | SUB_A | SUB_C | SUB_T)  # G
        self.assertEqual(codes[3], DELET | INS_3 | SUB_A | SUB_C | SUB_G)  # T


# ---------------------------------------------------------------------------
# Integration tests for import_mm
# ---------------------------------------------------------------------------

class TestImportMM(ut.TestCase):
    """ Integration tests: import a two-reference MM file and verify outputs. """

    def setUp(self):
        set_config(verbosity=Level.FATAL, exit_on_error=True)
        self._tmpdir = tempfile.TemporaryDirectory()
        root = Path(self._tmpdir.name)
        self._mm_path = root / "test.mm"
        self._mm_path.write_bytes(_build_mm_bytes())
        self._out_dir = root / "out"
        self._out_dir.mkdir()
        self._tmp_dir = root / "tmp"
        self._tmp_dir.mkdir()
        self._results = import_mm(
            self._mm_path,
            sample=_SAMPLE,
            out_dir=self._out_dir,
            tmp_dir=self._tmp_dir,
            branch=_BRANCH,
            min_reads=1,
            batch_size=100,
            insert3=True,
            write_read_names=False,
            brotli_level=5,
            force=False,
        )

    def tearDown(self):
        set_config()
        self._tmpdir.cleanup()

    def _load_report(self, ref: str) -> dict:
        return json.loads(_report_path(self._out_dir, ref).read_text())

    # --- top-level return value ---

    def test_returns_two_paths(self):
        self.assertEqual(len(self._results), 2)

    def test_result_paths_are_distinct(self):
        self.assertNotEqual(self._results[0], self._results[1])

    # --- output file existence: ref1 ---

    def test_ref1_report_exists(self):
        self.assertTrue(_report_path(self._out_dir, _REF1).is_file())

    def test_ref1_refseq_exists(self):
        ref1_dir = _report_path(self._out_dir, _REF1).parent
        self.assertTrue((ref1_dir / "refseq.brickle").is_file())

    def test_ref1_batch_exists(self):
        ref1_dir = _report_path(self._out_dir, _REF1).parent
        self.assertTrue((ref1_dir / "relate-batch-0.brickle").is_file())

    # --- output file existence: ref2 ---

    def test_ref2_report_exists(self):
        self.assertTrue(_report_path(self._out_dir, _REF2).is_file())

    def test_ref2_refseq_exists(self):
        ref2_dir = _report_path(self._out_dir, _REF2).parent
        self.assertTrue((ref2_dir / "refseq.brickle").is_file())

    def test_ref2_batch_exists(self):
        ref2_dir = _report_path(self._out_dir, _REF2).parent
        self.assertTrue((ref2_dir / "relate-batch-0.brickle").is_file())

    # --- report field values: ref1 ---

    def test_ref1_report_sample(self):
        self.assertEqual(self._load_report(_REF1)[SampleF.title], _SAMPLE)

    def test_ref1_report_ref(self):
        self.assertEqual(self._load_report(_REF1)[RefF.title], _REF1)

    def test_ref1_report_n_reads_xam(self):
        self.assertEqual(self._load_report(_REF1)[NumReadsXamF.title], 3)

    def test_ref1_report_n_reads_rel(self):
        # All 3 reads have valid coverage (end3 >= end5), so none are dropped.
        self.assertEqual(self._load_report(_REF1)[NumReadsRelF.title], 3)

    # --- report field values: ref2 ---

    def test_ref2_report_sample(self):
        self.assertEqual(self._load_report(_REF2)[SampleF.title], _SAMPLE)

    def test_ref2_report_ref(self):
        self.assertEqual(self._load_report(_REF2)[RefF.title], _REF2)

    def test_ref2_report_n_reads_xam(self):
        self.assertEqual(self._load_report(_REF2)[NumReadsXamF.title], 2)

    def test_ref2_report_n_reads_rel(self):
        self.assertEqual(self._load_report(_REF2)[NumReadsRelF.title], 2)

    # --- min_reads filtering ---

    def test_min_reads_skips_ref2(self):
        """ When min_reads exceeds ref2's read count, only ref1 is imported. """
        import shutil
        out2 = self._out_dir.parent / "out2"
        out2.mkdir()
        tmp2 = self._out_dir.parent / "tmp2"
        tmp2.mkdir()
        try:
            results = import_mm(
                self._mm_path,
                sample=_SAMPLE,
                out_dir=out2,
                tmp_dir=tmp2,
                branch=_BRANCH,
                min_reads=3,   # ref2 has only 2 reads → skipped
                batch_size=100,
                insert3=True,
                write_read_names=False,
                brotli_level=5,
                force=False,
            )
            self.assertEqual(len(results), 1)
            self.assertIn(_REF1, str(results[0]))
        finally:
            shutil.rmtree(out2, ignore_errors=True)
            shutil.rmtree(tmp2, ignore_errors=True)


if __name__ == "__main__":
    ut.main(verbosity=2)
