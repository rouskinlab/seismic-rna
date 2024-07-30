import os
import shutil
import unittest as ut
from pathlib import Path

from seismicrna.core.arg.cli import opt_out_dir, opt_sim_dir
from seismicrna.core.logs import get_config, set_config
from seismicrna.sim.total import run as sim_total_run
from seismicrna.wf import run as wf_run

STEPS = ["align", "relate", "mask", "cluster", "table", "fold", "graph"]


class TestWorkflow(ut.TestCase):
    OUT_DIR = Path(opt_out_dir.default).absolute()
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    LOG_DIR = Path("log").absolute()
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"

    def setUp(self):
        self.SIM_DIR.mkdir()
        self.OUT_DIR.mkdir()
        self.LOG_DIR.mkdir(exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        shutil.rmtree(self.OUT_DIR)
        shutil.rmtree(self.LOG_DIR)

    def test_wf_sim_paired_20000reads_2clusts(self):
        # Suppress warnings.
        config = get_config()
        set_config(verbose=0, quiet=1)
        # Simulate the data to be processed with wf.
        fastqs = sim_total_run(sim_dir=str(self.SIM_DIR),
                               sample=self.SAMPLE,
                               ref=self.REF,
                               refs=self.REFS,
                               fold_max=2,
                               num_reads=20000,
                               ends_var=0.001)
        sample_dir = self.SIM_DIR.joinpath("samples", self.SAMPLE)
        for fastq, mate in zip(fastqs, [1, 2], strict=True):
            self.assertEqual(fastq,
                             sample_dir.joinpath(f"{self.REF}_R{mate}.fq.gz"))
            self.assertTrue(os.path.isfile(fastq))
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        self.assertTrue(fasta.is_file())
        # Process the data with wf.
        wf_run(fasta=fasta,
               input_path=[],
               dmfastqx=fastqs,
               cluster=True,
               fold=True,
               quantile=0.95,
               export=True,
               out_dir=self.OUT_DIR)
        for step in STEPS:
            step_dir = self.OUT_DIR.joinpath(self.SAMPLE, step)
            self.assertTrue(step_dir.is_dir())
        # Restore the original logging configuration.
        set_config(*config)


if __name__ == "__main__":
    ut.main()
