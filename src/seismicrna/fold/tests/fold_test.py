import tempfile
import unittest as ut
from datetime import datetime
from pathlib import Path

import numpy as np

from seismicrna.core.arg import (
    FOLD_BACKEND_AUTO,
    FOLD_BACKEND_RNASTRUCTURE,
    FOLD_BACKEND_VIENNARNA,
    FOLD_ENERGY_METHOD_AUTO,
    FOLD_ENERGY_METHOD_DEIGAN,
    FOLD_ENERGY_METHOD_CORDERO,
    FOLD_ENERGY_METHOD_EDDY,
    PROBE_DMS,
    PROBE_SHAPE,
    PROBE_ETC,
    PROBE_NONE,
)
from seismicrna.core.extern import (
    RNASTRUCTURE_FOLD_CMD,
    RNASTRUCTURE_SHAPEKNOTS_CMD,
    VIENNA_RNAFOLD_CMD,
    VIENNA_RNASUBOPT_CMD,
    dependency_exists,
)
from seismicrna.core.logs import Level, set_config
from seismicrna.core.report import (
    FoldBackendF,
    FoldDryRunF,
    FoldEnergyMethodF,
    FoldIsolatedF,
    FoldMinFreeEnergyF,
    FoldTempF,
)
from seismicrna.core.seq.region import Region
from seismicrna.core.seq.xna import DNA
from seismicrna.fold.main import run as run_fold
from seismicrna.fold.report import FoldReport
from seismicrna.fold.rnastructure import guess_data_path
from seismicrna.filter.main import run as run_filter
from seismicrna.idmut.io import IDmutBatchIO, ReadNamesBatchIO, RefseqIO
from seismicrna.idmut.report import IDmutReport


REF = "ref"
REF_SEQ = DNA("GGGCGCAAAGCGCCCAAAAGGGCGCAAAGCGCCC")
SAMPLE = "sample"
N_READS = 6
READ_NAMES = [f"Read{i}" for i in range(N_READS)]
END5S = [[1] for _ in range(N_READS)]
END3S = [[len(REF_SEQ)] for _ in range(N_READS)]
# {position: {relation_byte: [read_indices]}}. Byte 32 is a substitution
# to C in the relationship encoding; the actual code does not matter for
# fold tests as long as the position table is non-trivial.
MUTS = {
    3: {32: [0, 3]},
    7: {32: [1, 4]},
    12: {32: [2, 5]},
    18: {32: [0, 1]},
    22: {32: [2, 3]},
    28: {32: [4, 5]},
}


def write_idmut(out_dir: Path) -> Path:
    """Build a minimal idmut report covering one batch of reads."""
    branches = dict()
    began = datetime.now()
    refseq = RefseqIO(sample=SAMPLE, branches=branches, ref=REF, refseq=REF_SEQ)
    _, refseq_checksum = refseq.save(out_dir)
    muts = {
        pos: {
            rel: np.array(reads, dtype=int) for rel, reads in MUTS.get(pos, {}).items()
        }
        for pos in range(1, len(REF_SEQ) + 1)
    }
    idmut_batch = IDmutBatchIO(
        sample=SAMPLE,
        branches=branches,
        region=Region(REF, REF_SEQ),
        batch=0,
        seg_end5s=np.array(END5S),
        seg_end3s=np.array(END3S),
        muts=muts,
    )
    _, idmut_checksum = idmut_batch.save(out_dir)
    name_batch = ReadNamesBatchIO(
        sample=SAMPLE, branches=branches, ref=REF, batch=0, names=np.array(READ_NAMES)
    )
    _, names_checksum = name_batch.save(out_dir)
    report = IDmutReport(
        sample=SAMPLE,
        branches=branches,
        ref=REF,
        min_mapq=0,
        min_phred=0,
        phred_enc=33,
        overhangs=True,
        insert3=True,
        ambindel=False,
        clip_end5=0,
        clip_end3=0,
        min_reads=0,
        n_reads_xam=N_READS,
        n_reads_rel=N_READS,
        n_batches=1,
        checksums={
            ReadNamesBatchIO.btype(): [names_checksum],
            IDmutBatchIO.btype(): [idmut_checksum],
        },
        refseq_checksum=refseq_checksum,
        began=began,
        ended=datetime.now(),
    )
    return report.save(out_dir)


def _datapath_ok() -> bool:
    try:
        guess_data_path()
    except Exception:
        return False
    return True


def _backend_available(backend: str, pseudoknots: bool) -> bool:
    if backend == FOLD_BACKEND_VIENNARNA:
        if pseudoknots:
            return False
        return dependency_exists(VIENNA_RNAFOLD_CMD) and dependency_exists(
            VIENNA_RNASUBOPT_CMD
        )
    if backend == FOLD_BACKEND_RNASTRUCTURE:
        if pseudoknots:
            return dependency_exists(RNASTRUCTURE_SHAPEKNOTS_CMD) and _datapath_ok()
        return dependency_exists(RNASTRUCTURE_FOLD_CMD) and _datapath_ok()
    return False


def _resolve(probe: str, backend: str, method: str) -> tuple[str, str]:
    if backend == FOLD_BACKEND_AUTO:
        backend = (
            FOLD_BACKEND_RNASTRUCTURE if probe == PROBE_DMS else FOLD_BACKEND_VIENNARNA
        )
    if method == FOLD_ENERGY_METHOD_AUTO:
        method = (
            FOLD_ENERGY_METHOD_CORDERO
            if probe == PROBE_DMS
            else FOLD_ENERGY_METHOD_EDDY
        )
    return backend, method


def _valid_combo(backend: str, method: str, pseudoknots: bool) -> bool:
    if pseudoknots and backend != FOLD_BACKEND_RNASTRUCTURE:
        return False
    if method == FOLD_ENERGY_METHOD_EDDY:
        return backend == FOLD_BACKEND_VIENNARNA
    if method == FOLD_ENERGY_METHOD_CORDERO:
        return backend == FOLD_BACKEND_RNASTRUCTURE
    return True


BACKEND_METHOD_COMBOS = [
    (FOLD_BACKEND_AUTO, FOLD_ENERGY_METHOD_AUTO, False),
    (FOLD_BACKEND_RNASTRUCTURE, FOLD_ENERGY_METHOD_DEIGAN, False),
    (FOLD_BACKEND_RNASTRUCTURE, FOLD_ENERGY_METHOD_CORDERO, False),
    (FOLD_BACKEND_RNASTRUCTURE, FOLD_ENERGY_METHOD_DEIGAN, True),
    (FOLD_BACKEND_RNASTRUCTURE, FOLD_ENERGY_METHOD_CORDERO, True),
    (FOLD_BACKEND_VIENNARNA, FOLD_ENERGY_METHOD_DEIGAN, False),
    (FOLD_BACKEND_VIENNARNA, FOLD_ENERGY_METHOD_EDDY, False),
]

# Defaults + each boolean flipped individually.
BOOLEAN_VARIANTS = [
    {"fold_table_region": False, "fold_mfe": False, "fold_isolated": False},
    {"fold_table_region": True, "fold_mfe": False, "fold_isolated": False},
    {"fold_table_region": False, "fold_mfe": True, "fold_isolated": False},
    {"fold_table_region": False, "fold_mfe": False, "fold_isolated": True},
]


class FoldCombinationsBase(ut.TestCase):
    """Run ``seismic fold`` over the parameter combination matrix."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._out_dir = Path(self._tmp.name)
        set_config(verbosity=Level.FATAL, exit_on_error=True)
        idmut_report = write_idmut(self._out_dir)
        self._filter_dirs = {
            probe: run_filter(
                [idmut_report], probe=probe, branch=probe, filter_pos_table=True
            )
            for probe in [PROBE_DMS, PROBE_SHAPE, PROBE_ETC, PROBE_NONE]
        }

    def tearDown(self):
        self._tmp.cleanup()
        self._tmp = None
        self._out_dir = None
        self._filter_dirs = {}
        set_config()

    def test_fold_combinations(self):
        for probe, filter_dirs in self._filter_dirs.items():
            for backend, method, pseudoknots in BACKEND_METHOD_COMBOS:
                resolved_backend, resolved_method = _resolve(probe, backend, method)
                if not _valid_combo(resolved_backend, resolved_method, pseudoknots):
                    continue
                if not _backend_available(resolved_backend, pseudoknots):
                    continue
                for booleans in BOOLEAN_VARIANTS:
                    for dry_run in (True, False):
                        with self.subTest(
                            probe=probe,
                            fold_backend=backend,
                            pseudoknots=pseudoknots,
                            fold_energy_method=method,
                            fold_dry_run=dry_run,
                            **booleans,
                        ):
                            report_files = run_fold(
                                filter_dirs,
                                fold_backend=backend,
                                pseudoknots=pseudoknots,
                                fold_energy_method=method,
                                fold_dry_run=dry_run,
                                tmp_pfx=self._tmp.name,
                                force=True,
                                **booleans,
                            )
                            self.assertGreaterEqual(len(report_files), 1)
                            report = FoldReport.load(report_files[0])
                            self.assertEqual(
                                report.get_field(FoldBackendF), resolved_backend
                            )
                            self.assertEqual(
                                report.get_field(FoldEnergyMethodF), resolved_method
                            )
                            self.assertEqual(report.get_field(FoldDryRunF), dry_run)
                            self.assertEqual(
                                report.get_field(FoldMinFreeEnergyF),
                                booleans["fold_mfe"],
                            )
                            self.assertEqual(
                                report.get_field(FoldIsolatedF),
                                booleans["fold_isolated"],
                            )
                            self.assertEqual(report.get_field(FoldTempF), 37.0)


if __name__ == "__main__":
    ut.main(verbosity=2)
