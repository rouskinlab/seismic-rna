import tempfile
import unittest as ut
from datetime import datetime
from pathlib import Path

import numpy as np

from seismicrna.core.arg import (FOLD_BACKEND_AUTO,
                                 FOLD_BACKEND_FOLD,
                                 FOLD_BACKEND_SHAPEKNOTS,
                                 FOLD_BACKEND_RNAFOLD,
                                 FOLD_ENERGY_METHOD_AUTO,
                                 FOLD_ENERGY_METHOD_DEIGAN,
                                 FOLD_ENERGY_METHOD_CORDERO,
                                 FOLD_ENERGY_METHOD_EDDY,
                                 PROBE_DMS,
                                 PROBE_SHAPE)
from seismicrna.core.extern import (RNASTRUCTURE_FOLD_CMD,
                                    RNASTRUCTURE_SHAPEKNOTS_CMD,
                                    VIENNA_RNAFOLD_CMD,
                                    VIENNA_RNASUBOPT_CMD,
                                    dependency_exists)
from seismicrna.core.logs import Level, set_config
from seismicrna.core.report import (FoldBackendF,
                                    FoldDryRunF,
                                    FoldEnergyMethodF,
                                    FoldIsolatedF,
                                    FoldMinFreeEnergyF,
                                    FoldTempF)
from seismicrna.core.seq.region import Region
from seismicrna.core.seq.xna import DNA
from seismicrna.fold.main import run as run_fold
from seismicrna.fold.report import FoldReport
from seismicrna.fold.rnastructure import guess_data_path
from seismicrna.mask.main import run as run_mask
from seismicrna.relate.io import RelateBatchIO, ReadNamesBatchIO, RefseqIO
from seismicrna.relate.report import RelateReport


REF = "ref"
REF_SEQ = DNA("GGGCGCAAAGCGCCCAAAAGGGCGCAAAGCGCCC")
SAMPLE = "sample"
N_READS = 6
READ_NAMES = [f"Read{i}" for i in range(N_READS)]
END5S = [[1] for _ in range(N_READS)]
END3S = [[len(REF_SEQ)] for _ in range(N_READS)]
# {position: {relation_byte: [read_indices]}}.  Byte 32 is a substitution
# to C in the relate-flag encoding; the actual code does not matter for
# fold tests as long as the position table is non-trivial.
MUTS = {3: {32: [0, 3]},
        7: {32: [1, 4]},
        12: {32: [2, 5]},
        18: {32: [0, 1]},
        22: {32: [2, 3]},
        28: {32: [4, 5]}}


def write_relate(out_dir: Path) -> Path:
    """ Build a minimal relate report covering one batch of reads. """
    branches = dict()
    began = datetime.now()
    refseq = RefseqIO(sample=SAMPLE,
                      branches=branches,
                      ref=REF,
                      refseq=REF_SEQ)
    _, refseq_checksum = refseq.save(out_dir)
    muts = {pos: {rel: np.array(reads, dtype=int)
                  for rel, reads in MUTS.get(pos, {}).items()}
            for pos in range(1, len(REF_SEQ) + 1)}
    relate_batch = RelateBatchIO(sample=SAMPLE,
                                 branches=branches,
                                 region=Region(REF, REF_SEQ),
                                 batch=0,
                                 seg_end5s=np.array(END5S),
                                 seg_end3s=np.array(END3S),
                                 muts=muts)
    _, relate_checksum = relate_batch.save(out_dir)
    name_batch = ReadNamesBatchIO(sample=SAMPLE,
                                  branches=branches,
                                  ref=REF,
                                  batch=0,
                                  names=np.array(READ_NAMES))
    _, names_checksum = name_batch.save(out_dir)
    report = RelateReport(sample=SAMPLE,
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
                          checksums={ReadNamesBatchIO.btype(): [names_checksum],
                                     RelateBatchIO.btype(): [relate_checksum]},
                          refseq_checksum=refseq_checksum,
                          began=began,
                          ended=datetime.now())
    return report.save(out_dir)


def _datapath_ok() -> bool:
    try:
        guess_data_path()
    except Exception:
        return False
    return True


def _backend_available(backend: str) -> bool:
    if backend == FOLD_BACKEND_RNAFOLD:
        return (dependency_exists(VIENNA_RNAFOLD_CMD)
                and dependency_exists(VIENNA_RNASUBOPT_CMD))
    if backend == FOLD_BACKEND_FOLD:
        return dependency_exists(RNASTRUCTURE_FOLD_CMD) and _datapath_ok()
    if backend == FOLD_BACKEND_SHAPEKNOTS:
        return (dependency_exists(RNASTRUCTURE_SHAPEKNOTS_CMD)
                and _datapath_ok())
    return False


def _resolve(probe: str, backend: str, method: str) -> tuple[str, str]:
    if backend == FOLD_BACKEND_AUTO:
        backend = (FOLD_BACKEND_FOLD if probe == PROBE_DMS
                   else FOLD_BACKEND_RNAFOLD)
    if method == FOLD_ENERGY_METHOD_AUTO:
        method = (FOLD_ENERGY_METHOD_CORDERO if probe == PROBE_DMS
                  else FOLD_ENERGY_METHOD_EDDY)
    return backend, method


def _valid_combo(backend: str, method: str) -> bool:
    if method == FOLD_ENERGY_METHOD_EDDY:
        return backend == FOLD_BACKEND_RNAFOLD
    if method == FOLD_ENERGY_METHOD_CORDERO:
        return backend in (FOLD_BACKEND_FOLD, FOLD_BACKEND_SHAPEKNOTS)
    return True


BACKEND_METHOD_COMBOS = [
    (FOLD_BACKEND_AUTO, FOLD_ENERGY_METHOD_AUTO),
    (FOLD_BACKEND_FOLD, FOLD_ENERGY_METHOD_DEIGAN),
    (FOLD_BACKEND_FOLD, FOLD_ENERGY_METHOD_CORDERO),
    (FOLD_BACKEND_SHAPEKNOTS, FOLD_ENERGY_METHOD_DEIGAN),
    (FOLD_BACKEND_SHAPEKNOTS, FOLD_ENERGY_METHOD_CORDERO),
    (FOLD_BACKEND_RNAFOLD, FOLD_ENERGY_METHOD_DEIGAN),
    (FOLD_BACKEND_RNAFOLD, FOLD_ENERGY_METHOD_EDDY),
]

# Defaults + each boolean flipped individually.
BOOLEAN_VARIANTS = [
    {"fold_full": True, "fold_mfe": False, "fold_isolated": False},
    {"fold_full": False, "fold_mfe": False, "fold_isolated": False},
    {"fold_full": True, "fold_mfe": True, "fold_isolated": False},
    {"fold_full": True, "fold_mfe": False, "fold_isolated": True},
]


class FoldCombinationsBase(ut.TestCase):
    """ Run ``seismic fold`` over the parameter combination matrix. """

    # Override in concrete subclasses.
    PROBE: str | None = None

    def setUp(self):
        if self.PROBE is None:
            self.skipTest("abstract base class")
        self._tmp = tempfile.TemporaryDirectory()
        self._out_dir = Path(self._tmp.name)
        set_config(verbosity=Level.FATAL, exit_on_error=True)
        self._tmp_pfx = self._out_dir / "tmp"
        relate_report = write_relate(self._out_dir)
        self._mask_dirs = run_mask([relate_report],
                                   probe=self.PROBE,
                                   mask_pos_table=True,
                                   tmp_pfx=self._tmp_pfx)
        if not self._mask_dirs:
            self.skipTest("mask produced no output")

    def tearDown(self):
        if self._tmp is not None:
            self._tmp.cleanup()
        self._tmp = None
        self._out_dir = None
        self._mask_dirs = None
        set_config()

    def test_fold_combinations(self):
        for backend, method in BACKEND_METHOD_COMBOS:
            resolved_backend, resolved_method = _resolve(self.PROBE,
                                                         backend,
                                                         method)
            if not _valid_combo(resolved_backend, resolved_method):
                continue
            if not _backend_available(resolved_backend):
                continue
            for booleans in BOOLEAN_VARIANTS:
                for dry_run in (True, False):
                    with self.subTest(probe=self.PROBE,
                                      fold_backend=backend,
                                      fold_energy_method=method,
                                      fold_dry_run=dry_run,
                                      **booleans):
                        report_files = run_fold(
                            self._mask_dirs,
                            fold_backend=backend,
                            fold_energy_method=method,
                            fold_dry_run=dry_run,
                            tmp_pfx=self._tmp_pfx,
                            force=True,
                            **booleans,
                        )
                        self.assertGreaterEqual(len(report_files), 1)
                        report = FoldReport.load(report_files[0])
                        self.assertEqual(report.get_field(FoldBackendF),
                                         resolved_backend)
                        self.assertEqual(report.get_field(FoldEnergyMethodF),
                                         resolved_method)
                        self.assertEqual(report.get_field(FoldDryRunF),
                                         dry_run)
                        self.assertEqual(report.get_field(FoldMinFreeEnergyF),
                                         booleans["fold_mfe"])
                        self.assertEqual(report.get_field(FoldIsolatedF),
                                         booleans["fold_isolated"])
                        self.assertEqual(report.get_field(FoldTempF), 37.0)


class TestFoldDMS(FoldCombinationsBase):
    PROBE = PROBE_DMS


class TestFoldSHAPE(FoldCombinationsBase):
    PROBE = PROBE_SHAPE


if __name__ == "__main__":
    ut.main(verbosity=2)
