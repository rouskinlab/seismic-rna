import json
import struct
import tempfile
import unittest as ut
from pathlib import Path

from seismicrna.core import path
from seismicrna.core.logs import Level, logger, set_config
from seismicrna.core.rel.code import DELET, INS_3, INS_5, SUB_A, SUB_C, SUB_G, SUB_N, SUB_T
from seismicrna.core.report import NumReadsXamF, NumReadsRelF, RefF, SampleF
from seismicrna.core.seq import DNA
from seismicrna.importmm.write import _build_mut_codes, import_mm
from seismicrna.relate.report import RelateReport

_SAMPLE = "test_sample"
_BRANCH = ""
_REF1 = "ref1"
_REF2 = "ref2"


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
            relate_pos_table=False,
            relate_read_table=False,
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
                relate_pos_table=False,
                relate_read_table=False,
                brotli_level=5,
                force=False,
            )
            self.assertEqual(len(results), 1)
            self.assertIn(_REF1, str(results[0]))
        finally:
            shutil.rmtree(out2, ignore_errors=True)
            shutil.rmtree(tmp2, ignore_errors=True)


# ---------------------------------------------------------------------------
# rf-count integration test
# ---------------------------------------------------------------------------

class TestRFCountIntegration(ut.TestCase):
    """Verify importmm and relate produce identical position/read-table counts."""

    _SAMPLE = "rfcount_integ"
    _REF = "integ_ref"

    def test_matches_relate(self):
        import shutil
        import subprocess
        import pandas as pd
        from seismicrna.align import run as run_align
        from seismicrna.relate import run as run_relate
        from seismicrna.sim.ref import run as run_sim_ref
        from seismicrna.core.arg.cli import PROBE_DMS
        from seismicrna.sim.fold import run as run_sim_fold
        from seismicrna.sim.params import run as run_sim_params
        from seismicrna.sim.fastq import run as run_sim_fastq

        if shutil.which("rf-count") is None:
            logger.warning("Skipped test of importing Mutation Map files from "
                           "RNAFramework because RNAFramework is not installed")
            return  # rf-count not installed; trivially passes

        set_config(verbosity=Level.FATAL, exit_on_error=True)
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                sim_dir = root / "sim"
                sim_dir.mkdir()
                out_dir = root / "out"
                out_dir.mkdir()
                out_importmm = root / "out_importmm"
                out_importmm.mkdir()
                rf_out = root / "rf_out"
                rf_out.mkdir()
                tmp_pfx = str(root / "tmp")

                # 1. Simulate a 200-bp reference (longer than the 151-bp reads
                #    so most reads align without excessive soft-clipping).
                reflen = 200
                run_sim_ref(refs=self._REF, ref=self._REF, reflen=reflen,
                            sim_dir=sim_dir, seed=42)
                fasta = sim_dir / "refs" / f"{self._REF}.fa"

                # 2. Fold (single MFE structure; requires RNAstructure/ViennaRNA).
                run_sim_fold(fasta, probe=PROBE_DMS, sim_dir=sim_dir,
                             fold_mfe=True, fold_max=1,
                             num_cpus=1, tmp_pfx=tmp_pfx)

                # 3. Sim params: substitutions only, no deletions, no low-quality
                # bases. Setting each base's substitution fractions to sum to
                # 1.0 (one-third each) ensures the deletion probability is 0.
                # "loq" set to 0.0 ensures no low-quality (N) base calls.
                # fold writes to params/{ref}/full/simulated.ct ("full" is the region).
                param_dir = sim_dir / "params" / self._REF / "full"
                pmut = [
                    ("am", 0.05), ("cm", 0.05), ("gm", 0.05), ("tm", 0.05),
                    ("loq", 0.0),
                    # Substitution fractions sum to 1.0 → deletion rate = 0.
                    ("ac", 1/3), ("ag", 1/3), ("at", 1/3),
                    ("ca", 1/3), ("cg", 1/3), ("ct", 1/3),
                    ("ga", 1/3), ("gc", 1/3), ("gt", 1/3),
                    ("ta", 1/3), ("tc", 1/3), ("tg", 1/3),
                ]
                run_sim_params(ct_file=[param_dir / "simulated.ct"],
                               pmut_paired=pmut, pmut_unpaired=pmut,
                               probe=PROBE_DMS,
                               seed=42)

                # 4. Simulate 500 single-end reads with read_length=50.
                # read_length < reflen (200) so every read fits entirely within
                # the reference; no adapter padding is needed, so all reads
                # are clean reference-derived sequence.  Single-end avoids the
                # factor-of-2 counting difference between rf-count (counts each
                # BAM record) and relate (counts each paired read as one unit).
                fastqs = run_sim_fastq(input_path=(), param_dir=(param_dir,),
                                       sample=self._SAMPLE, num_reads=500,
                                       paired_end=False, read_length=50,
                                       seed=42)
                samples_dir = fastqs[0].parent.parent.parent

                # 5. Align → sorted, indexed BAM.
                # run_align returns output directories, not BAM file paths.
                # Use dmfastqz for single-end demultiplexed reads.
                align_dirs = run_align(fasta, dmfastqz=[samples_dir],
                                       out_dir=out_dir,
                                       min_mapq=0, min_reads=1,
                                       num_cpus=1, force=True,
                                       tmp_pfx=tmp_pfx)
                bam = next(p for d in align_dirs for p in d.rglob("*.bam"))

                # 6. rf-count → MM binary file.
                # -q 0 / -mq 0: accept all base and mapping qualities so that
                # neither threshold diverges from relate's min_phred=0, min_mapq=0.
                subprocess.run(
                    ["rf-count", "-m", "-mm",
                     "-f", str(fasta),
                     "-o", str(rf_out),
                     "-q", "0", "-mq", "0",
                     "-ow", str(bam)],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                mm_file = next(rf_out.rglob("*.mm"))

                # 7. Import MM → relate-format output with position and read tables.
                tmp_mm = root / "tmp_mm"
                tmp_mm.mkdir()
                import_mm(
                    mm_file,
                    sample=self._SAMPLE,
                    out_dir=out_importmm,
                    tmp_dir=tmp_mm,
                    branch="",
                    min_reads=1,
                    batch_size=1000,
                    insert3=True,
                    write_read_names=False,
                    relate_pos_table=True,
                    relate_read_table=True,
                    brotli_level=5,
                    force=False,
                )

                # 8. Relate on the same BAM.
                # clip_end5/3=0: rf-count does not clip read ends.
                # min_phred=0:   sim reads have perfect quality, matching -q 0 above.
                # ambindel=False: consistent with importmm sentinel value.
                run_relate(fasta, align_dirs,
                           out_dir=out_dir,
                           clip_end5=0, clip_end3=0,
                           insert3=True, ambindel=False, overhangs=True,
                           min_mapq=0, min_phred=0, min_reads=1,
                           relate_pos_table=True, relate_read_table=True,
                           relate_cx=False,
                           num_cpus=1, brotli_level=5,
                           force=False, tmp_pfx=tmp_pfx)

                # 9. Compare Mutated counts element-by-element in both tables.
                #
                # rf-count's MM format stores only reads that have ≥1 mutation;
                # perfect-match reads are absent.  As a result, importmm's
                # Covered and Matched totals will be lower than relate's (which
                # processes all aligned reads), and read-table row counts will
                # differ.  Only the Mutated column is directly comparable: both
                # tools count the exact same mutations from the same BAM file.
                #
                # Position table: both tables share the same (Position, Base)
                # index, so the Mutated value at each index must be identical.
                mm_pos = next(out_importmm.rglob("*-position-table.csv"))
                rel_pos = next(out_dir.rglob("*-position-table.csv"))
                mm_pos_tbl = pd.read_csv(mm_pos, index_col=[0, 1], header=0)
                rel_pos_tbl = pd.read_csv(rel_pos, index_col=[0, 1], header=0)
                self.assertEqual(mm_pos_tbl.index.tolist(),
                                 rel_pos_tbl.index.tolist(),
                                 msg="Position tables have different indexes")
                mm_pos_muts = mm_pos_tbl["Mutated"].astype(int).tolist()
                rel_pos_muts = rel_pos_tbl["Mutated"].astype(int).tolist()
                self.assertEqual(len(mm_pos_muts), reflen)
                self.assertListEqual(mm_pos_muts, rel_pos_muts,
                                 msg="Position table Mutated counts differ "
                                     "element-by-element")

                # Read table: importmm contains only the mutated reads (those
                # rf-count put in the MM file); relate contains all reads.
                # Read names and row order differ, so compare sorted per-read
                # Mutated counts: importmm's full list vs relate's mutated-only
                # subset (Mutated > 0), both sorted ascending.
                mm_read = next(out_importmm.rglob("*-read-table.csv.gz"))
                rel_read = next(out_dir.rglob("*-read-table.csv.gz"))
                mm_read_tbl = pd.read_csv(mm_read, index_col=0, header=0)
                rel_read_tbl = pd.read_csv(rel_read, index_col=0, header=0)
                mm_read_muts = sorted(mm_read_tbl["Mutated"].astype(int))
                rel_read_muts = sorted(
                    rel_read_tbl.loc[rel_read_tbl["Mutated"] > 0,
                                     "Mutated"].astype(int)
                )
                self.assertListEqual(mm_read_muts, rel_read_muts,
                                 msg="Read table per-read Mutated counts "
                                     "differ (sorted)")
        finally:
            set_config()


if __name__ == "__main__":
    ut.main(verbosity=2)
