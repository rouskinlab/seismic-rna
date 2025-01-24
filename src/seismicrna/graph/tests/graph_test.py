import shutil
import unittest as ut
from pathlib import Path

from seismicrna.core.arg.cli import opt_sim_dir, GROUP_BY_K, KEY_PEARSON
from seismicrna.core.logs import Level, set_config, get_config
from seismicrna.graph.corroll import RollingCorrelationRunner
from seismicrna.graph.delprof import DeltaProfileRunner
from seismicrna.graph.histpos import PositionHistogramRunner
from seismicrna.graph.scatter import ScatterRunner
from seismicrna.graph.snrroll import RollingSNRRunner
from seismicrna.sim import (ref as ref_mod,
                            fold as fold_mod,
                            params as params_mod,
                            relate as relate_mod)
from seismicrna.table import run as run_table
from seismicrna.wf import run as run_wf


class GraphTests(ut.TestCase):
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    REFS = "test_refs"
    REF = "test_ref"
    REF_LEN = 90
    NUM_STRUCTS = 2
    SAMPLE1 = "test_sample1"
    SAMPLE2 = "test_sample2"
    SAMPLES = [SAMPLE1, SAMPLE2]

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   raise_on_error=True)
        self.SIM_DIR.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        set_config(*self._config)

    def test_all_graphs(self):
        fasta = ref_mod.run(sim_dir=self.SIM_DIR,
                            refs=self.REFS,
                            ref=self.REF,
                            reflen=self.REF_LEN)
        ct_file, = fold_mod.run(fasta,
                                sim_dir=self.SIM_DIR,
                                fold_max=self.NUM_STRUCTS)
        params_mod.run(ct_file=[ct_file],
                       length_fmean=1.)
        relate_dirs = [
            relate_dir
            for sample in self.SAMPLES
            for relate_dir in relate_mod.run(param_dir=[ct_file.parent],
                                             sample=sample,
                                             paired_end=False,
                                             read_length=self.REF_LEN,
                                             num_reads=(2 ** 15))
        ]
        run_table(relate_dirs)
        graph_kwargs = dict(cgroup=GROUP_BY_K,
                            csv=True,
                            html=True,
                            svg=True,
                            pdf=True,
                            png=True,
                            max_procs=1,
                            force=False)
        rel_graph_kwargs = graph_kwargs | dict(rels=("m",),
                                               use_ratio=True,
                                               quantile=0.0)
        pair_graph_kwargs = rel_graph_kwargs | dict(out_dir=self.SIM_DIR,
                                                    comppair=True,
                                                    compself=False)
        run_wf(fasta,
               tuple(relate_dirs),
               out_dir=str(self.SIM_DIR),
               cluster=True,
               max_clusters=self.NUM_STRUCTS,
               max_em_iter=30,
               jackpot=False,
               fold=True,
               quantile=0.95,
               graph_mprof=True,
               graph_tmprof=True,
               graph_ncov=True,
               graph_mhist=True,
               graph_abundance=True,
               graph_giniroll=True,
               graph_roc=True,
               graph_aucroll=True,
               graph_poscorr=True,
               graph_mutdist=True,
               mutdist_null=True,
               **graph_kwargs)
        RollingSNRRunner.run([self.SIM_DIR],
                             window=21,
                             winmin=7,
                             **rel_graph_kwargs)
        PositionHistogramRunner.run([self.SIM_DIR],
                                    hist_bins=37,
                                    hist_margin=0.01,
                                    **rel_graph_kwargs)
        RollingCorrelationRunner.run([self.SIM_DIR],
                                     metric=KEY_PEARSON,
                                     window=21,
                                     winmin=7,
                                     **pair_graph_kwargs)
        DeltaProfileRunner.run([self.SIM_DIR],
                               **pair_graph_kwargs)
        ScatterRunner.run([self.SIM_DIR],
                          metric=KEY_PEARSON,
                          **pair_graph_kwargs)


if __name__ == "__main__":
    ut.main()
