import tempfile
import unittest as ut
from pathlib import Path
from unittest.mock import patch

from seismicrna.align.write import (calc_flags_sep_strands,
                                    fq_pipeline,
                                    FastqUnit)
from seismicrna.align.xamops import run_bowtie2_build
from seismicrna.core.logs import Level, get_config, set_config, restore_config
from seismicrna.core.ngs import (FLAG_PAIRED,
                                 FLAG_PROPER,
                                 FLAG_FIRST,
                                 FLAG_SECOND,
                                 FLAG_REVERSE)

from seismicrna.sim.total import run as run_sim_total


class TestCalcFlags(ut.TestCase):

    @restore_config
    def test_f1r2_paired(self):
        set_config(verbosity=Level.ERROR)
        for bt2_mixed in [False, True]:
            expect = (([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND | FLAG_REVERSE,
                        FLAG_FIRST]),
                      ([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND,
                        FLAG_FIRST | FLAG_REVERSE]))
            result = calc_flags_sep_strands(True, True, bt2_mixed)
            self.assertEqual(result, expect)

    def test_f1r2_single(self):
        expect = (([0],
                   [FLAG_REVERSE | FLAG_PAIRED]),
                  ([FLAG_REVERSE],
                   [FLAG_PAIRED]))
        for mixed in [False, True]:
            result = calc_flags_sep_strands(True, False, mixed)
            self.assertEqual(result, expect)

    @restore_config
    def test_f2r1_paired(self):
        set_config(verbosity=Level.ERROR)
        for bt2_mixed in [False, True]:
            expect = (([FLAG_FIRST | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND,
                        FLAG_FIRST | FLAG_REVERSE]),
                      ([FLAG_FIRST | FLAG_PAIRED | FLAG_PROPER,
                        FLAG_SECOND | FLAG_REVERSE | FLAG_PAIRED | FLAG_PROPER],
                       [FLAG_SECOND | FLAG_REVERSE,
                        FLAG_FIRST]))
            result = calc_flags_sep_strands(False, True, bt2_mixed)
            self.assertEqual(result, expect)

    def test_f2r1_single(self):
        expect = (([FLAG_REVERSE],
                   [FLAG_PAIRED]),
                  ([0],
                   [FLAG_REVERSE | FLAG_PAIRED]))
        for mixed in [False, True]:
            result = calc_flags_sep_strands(False, False, mixed)
            self.assertEqual(result, expect)


class TestFqPipeline(ut.TestCase):

    def setUp(self):
        self.maxDiff = 10000
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   exit_on_error=True)

    def tearDown(self):
        set_config(*self._config)

    def test_fq_pipeline(self):
        for sep_strands in [False, True]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                sim_dir = tmpdir / 'sim'
                sim_dir.mkdir()
                ref_name = 'testref'
                sample_name = 'sample1'
                profile_name = 'simulated'
                # Run the total simulation pipeline
                fastqs = run_sim_total(
                    sim_dir=sim_dir,
                    tmp_pfx=tmpdir,
                    sample=sample_name,
                    refs=ref_name,
                    ref=ref_name,
                    reflen=50,
                    profile_name=profile_name,
                    fold_coords=[],
                    fold_primers=[],
                    fold_regions_file=None,
                    fold_constraint=None,
                    fold_temp=310.15,
                    fold_md=0,
                    fold_mfe=True,
                    fold_max=1,
                    fold_percent=10.0,
                    pmut_paired=[],
                    pmut_unpaired=[],
                    vmut_paired=0.01,
                    vmut_unpaired=0.01,
                    center_fmean=0.5,
                    center_fvar=0.0,
                    length_fmean=1.0,
                    length_fvar=0.0,
                    clust_conc=1.0,
                    paired_end=False,
                    read_length=50,
                    reverse_fraction=0.0,
                    probe="SHAPE",
                    min_mut_gap=None,
                    mut_collisions='auto',
                    fq_gzip=False,
                    num_reads=1024,
                    keep_tmp=False,
                    force=True,
                    num_cpus=1,
                    seed=1
                )
                self.assertEqual(len(fastqs), 1)
                fastq_path = Path(fastqs[0])
                # Prepare other required arguments
                out_dir = tmpdir / "out"
                tmp_dir = tmpdir / "tmp"
                out_dir.mkdir()
                tmp_dir.mkdir()
                # Create FastqUnit
                fq_unit = FastqUnit(fastqz=fastq_path,
                                    phred_enc=33,
                                    one_ref=True)
                # Build the Bowtie2 index.
                fasta = sim_dir / "refs" / f"{ref_name}.fa"
                bowtie2_index = fasta.with_suffix("")
                run_bowtie2_build(fasta, bowtie2_index)
                # Run fq_pipeline
                result = fq_pipeline(
                    fq_inp=fq_unit,
                    fasta=fasta,
                    bowtie2_index=bowtie2_index,
                    out_dir=out_dir,
                    tmp_dir=tmp_dir,
                    keep_tmp=False,
                    branches={},
                    fastp=False,
                    fastp_5=False,
                    fastp_3=False,
                    fastp_w=4,
                    fastp_m=20,
                    fastp_poly_g='',
                    fastp_poly_g_min_len=10,
                    fastp_poly_x=False,
                    fastp_poly_x_min_len=10,
                    fastp_adapter_trimming=False,
                    fastp_adapter_1='',
                    fastp_adapter_2='',
                    fastp_adapter_fasta=None,
                    fastp_detect_adapter_for_pe=False,
                    fastp_min_length=15,
                    bt2_local=False,
                    bt2_discordant=False,
                    bt2_mixed=False,
                    bt2_dovetail=False,
                    bt2_contain=False,
                    bt2_score_min_e2e="L,-1,-0.8",
                    bt2_score_min_loc="L,1,0.8",
                    bt2_i=1,
                    bt2_x=600,
                    bt2_gbar=4,
                    bt2_l=22,
                    bt2_s="L,1,0.1",
                    bt2_d=15,
                    bt2_r=2,
                    bt2_dpad=15,
                    bt2_orient="fr",
                    bt2_un=False,
                    seed=1,
                    min_mapq=10,
                    min_reads=1,
                    sep_strands=sep_strands,
                    f1r2_fwd=False,
                    rev_label="-minus",
                    num_cpus=1
                )
                self.assertIsNotNone(result)


if __name__ == "__main__":
    ut.main(verbosity=2)
