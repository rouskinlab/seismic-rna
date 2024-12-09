import os
import shutil
import unittest as ut
from pathlib import Path

from seismicrna.core.arg.cli import opt_out_dir, opt_sim_dir
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.sim.fastq import run as run_sim_fastq
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.sim.total import run as run_sim_total
from seismicrna.wf import run as wf_run

STEPS = ["align", "relate", "mask", "cluster", "fold", "graph"]


class TestWorkflow(ut.TestCase):
    OUT_DIR = Path(opt_out_dir.default).absolute()
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None

    def setUp(self):
        self.SIM_DIR.mkdir()
        self.OUT_DIR.mkdir()
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   raise_on_error=True)

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        shutil.rmtree(self.OUT_DIR)
        set_config(**self._config._asdict())

    def test_wf_sim_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        fastqs = run_sim_total(sim_dir=str(self.SIM_DIR),
                               sample=self.SAMPLE,
                               ref=self.REF,
                               refs=self.REFS,
                               fold_max=2,
                               num_reads=20000,
                               center_fvar=0.001,
                               length_fvar=0.001)
        sample_dir = self.SIM_DIR.joinpath("samples", self.SAMPLE)
        for fastq, mate in zip(fastqs, [1, 2], strict=True):
            self.assertEqual(fastq,
                             sample_dir.joinpath(f"{self.REF}_R{mate}.fq.gz"))
            self.assertTrue(os.path.isfile(fastq))
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        self.assertTrue(fasta.is_file())
        # Process the data with wf.
        wf_run(fasta=fasta,
               input_path=(),
               dmfastqx=fastqs,
               cluster=True,
               fold=True,
               quantile=0.95,
               export=True,
               out_dir=self.OUT_DIR)
        for step in STEPS:
            step_dir = self.OUT_DIR.joinpath(self.SAMPLE, step)
            self.assertTrue(step_dir.is_dir())

    def test_wf_sim_2samples_2refs_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        samples = ["sample1", "sample2"]
        refs = ["refA", "refB"]
        samples_dir = self.SIM_DIR.joinpath("samples")
        all_fastas = list()
        for ref in refs:
            run_sim_ref(refs=ref, ref=ref)
            fasta = self.SIM_DIR.joinpath("refs", f"{ref}.fa")
            all_fastas.append(fasta)
            run_sim_fold(fasta, fold_max=2)
            param_dir = self.SIM_DIR.joinpath("params", ref, "full")
            ct_file = param_dir.joinpath("simulated.ct")
            run_sim_params(ct_file=(ct_file,))
            for sample in samples:
                fastqs = run_sim_fastq(input_path=(),
                                       param_dir=(param_dir,),
                                       sample=sample,
                                       num_reads=10000)
                sample_dir = samples_dir.joinpath(sample)
                for fastq, mate in zip(fastqs, [1, 2], strict=True):
                    self.assertEqual(
                        fastq,
                        sample_dir.joinpath(f"{ref}_R{mate}.fq.gz")
                    )
                    self.assertTrue(os.path.isfile(fastq))
        # Merge the FASTA files for all references.
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        with open(fasta, "x") as f:
            for ref in all_fastas:
                with open(ref) as r:
                    f.write(r.read())
        # Process the data with wf.
        wf_run(fasta=fasta,
               input_path=(),
               dmfastqx=(samples_dir,),
               cluster=True,
               fold=True,
               quantile=0.95,
               export=True,
               out_dir=self.OUT_DIR)
        for sample in samples:
            for step in STEPS:
                step_dir = self.OUT_DIR.joinpath(sample, step)
                self.assertTrue(step_dir.is_dir())


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
