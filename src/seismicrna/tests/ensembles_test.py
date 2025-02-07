import shutil
import unittest as ut
from itertools import product
from pathlib import Path

import numpy as np

from seismicrna.cluster.data import JoinClusterMutsDataset
from seismicrna.core.arg.cli import opt_sim_dir
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.rna.convert import db_to_ct
from seismicrna.core.seq.fasta import write_fasta
from seismicrna.core.seq.region import FULL_NAME
from seismicrna.core.seq.xna import DNA
from seismicrna.ensembles import (calc_regions,
                                  run as run_ensembles)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.relate import run as sim_relate

rng = np.random.default_rng()


class TestCalcRegions(ut.TestCase):

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR)

    def tearDown(self):
        set_config(self._config)

    def test_region_min_overlap_25(self):
        result = calc_regions(41, 145, 60, 0.25)
        expect = [(41, 100), (86, 145)]
        self.assertEqual(result, expect)

    def test_region_min_overlap_75(self):
        result = calc_regions(41, 145, 60, 0.75)
        expect = [(41, 100), (56, 115), (71, 130), (86, 145)]
        self.assertEqual(result, expect)

    def test_region_not_divisible(self):
        result = calc_regions(41, 170, 60, 0.25)
        expect = [(41, 100), (76, 135), (111, 170)]
        self.assertEqual(result, expect)

    def test_region_length_larger(self):
        result = calc_regions(41, 49, 52, 0.9)
        expect = [(41, 49)]
        self.assertEqual(result, expect)

    def test_total_region_length_1(self):
        result = calc_regions(41, 41, 52, 0.9)
        expect = [(41, 41)]
        self.assertEqual(result, expect)


@ut.skip("Takes a very long time")
class TestEnsembles(ut.TestCase):
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    SAMPLE = "test_sample"
    REFS = "test_refs"
    REF = "test_ref"
    PROFILE = "ensembles"

    # Folding modules of the reference sequence (each 60 nt).
    MODULES = [
        ("TGACGAACAACGTGTTTGTGAACCATATAGGTAAACGCTGAATGCGTTCGCGCGGAGGGT",
         ["..(((((((...)))))))..(((......((.(((((.....))))).))......)))",
          "...(((...(((((((((((.(((.....)))...))).))))))))))).((.....))"]),
        ("TTTGCAGGAAGATGGTCAACTCTACACCTAGTTTTTACCAGTCCACAAGAGTTTGAACTG",
         [".(((..(((...((((.((..(((....)))..)).)))).))).))).(((....)))."]),
        ("GTGCCTTAACCTGAGTACGCCCATATCATGGGAGACATTACAACTCAAATTCTAGGTGTG",
         ["..((((.....(((((...(((((...)))))..........)))))......))))...",
          "((((.(((...)))))))...((((((.((((..................))))))))))"]),
    ]

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.DETAIL,
                   log_file_path=None,
                   exit_on_error=True)
        self.SIM_DIR.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        set_config(*self._config)

    @classmethod
    def sim_data(cls, module_nums: list[int], read_length: int):
        # Assemble and write the reference sequence.
        modules = dict(cls.MODULES[m] for m in module_nums)
        refseq = DNA("".join(modules.keys()))
        refs_dir = cls.SIM_DIR.joinpath("refs")
        refs_dir.mkdir()
        fasta = refs_dir.joinpath(f"{cls.REFS}.fa")
        write_fasta(fasta, [(cls.REF, refseq)])
        # Assemble and write the secondary structures.
        structures = list(map("".join, product(*modules.values())))
        param_dir = cls.SIM_DIR.joinpath("params", cls.REF, FULL_NAME)
        param_dir.mkdir(parents=True)
        db_file = param_dir.joinpath(f"{cls.PROFILE}.db")
        with open(db_file, "x") as f:
            for i, struct in enumerate(structures):
                if i == 0:
                    f.write(f">structure0\n{refseq.tr()}\n{structures[0]}\n")
                else:
                    f.write(f">structure{i}\n{structures[i]}\n")
        ct_file = db_to_ct(db_file)
        # Simulate data.
        sim_params(ct_file=[ct_file],
                   # Make pmut_unpaired for A and C large (17%) so that
                   # most reads get at least two mutations despite being
                   # short and are thus useful for clustering.
                   pmut_unpaired=[("am", 1 / 6), ("cm", 1 / 6)],
                   # Make all reads the same length.
                   length_fmean=(read_length / len(refseq)),
                   length_fvar=0.,
                   # Make clust_conc very large so that the proportion
                   # of each cluster is approximately equal, which makes
                   # clustering easier.
                   clust_conc=1000.)
        relate_dirs = sim_relate(param_dir=[param_dir],
                                 sample=cls.SAMPLE,
                                 profile_name=cls.PROFILE,
                                 num_reads=200000,
                                 paired_end=False,
                                 brotli_level=0)
        return relate_dirs

    def run_ensembles(self,
                      relate_dirs: list[Path],
                      region_length: int,
                      region_min_overlap: float,
                      expect_ks: list[int],
                      expect_regions: list[list[str]],
                      **kwargs):
        joined = f"ensembles-{region_length}_"
        join_dirs = run_ensembles(relate_dirs,
                                  joined=joined,
                                  region_length=region_length,
                                  region_min_overlap=region_min_overlap,
                                  mask_pos_table=False,
                                  mask_read_table=False,
                                  # Calculating the cluster tables here
                                  # makes group_clusters() run faster.
                                  cluster_pos_table=True,
                                  cluster_abundance_table=True,
                                  jackpot=False,
                                  brotli_level=0,
                                  **kwargs)
        cluster_dirs = sorted(d for d in join_dirs
                              if d.parent.parent.name == "cluster")
        expect_dirs = [self.SIM_DIR.joinpath("samples",
                                             self.SAMPLE,
                                             "cluster",
                                             self.REF,
                                             f"{joined}{i}")
                       for i in range(1, len(expect_ks) + 1)]
        self.assertListEqual(cluster_dirs,
                             list(expect_dirs))
        for expect_dir, expect_k, expect_regs in zip(expect_dirs,
                                                     expect_ks,
                                                     expect_regions,
                                                     strict=True):
            report_file = expect_dir.joinpath("cluster-report.json")
            dataset = JoinClusterMutsDataset(report_file)
            self.assertListEqual(dataset.ks, [expect_k])
            self.assertListEqual(sorted(dataset.region_names),
                                 sorted(expect_regs))

    def test_modules012_read180(self):
        relate_dirs = self.sim_data([0, 1, 2], 180)
        # The regions are 1-60, 21-80, 41-100, 61-120, 81-140, 101-160,
        # and 121-180.
        # Regions 1-60, 21-80, and 41-100 overlap module 0 (1-60) and
        # should form 2 structures.
        # Region 61-120 coincides exactly with module 1 and should form
        # 1 structure.
        # Regions 81-140, 101-160, and 121-180 overlap module 2 (121-180)
        # and should form 2 structures.
        self.run_ensembles(relate_dirs, 60, (2 / 3),
                           [2, 1, 2],
                           [["1-60", "21-80", "41-100"],
                            ["61-120"],
                            ["81-140", "101-160", "121-180"]])
        # The regions are 1-90, 31-120, 61-150, and 91-180.
        # Every region includes a part of module 0 (1-60) and module 1
        # (121-180), but not both, so every region forms 2 structures.
        self.run_ensembles(relate_dirs, 90, (2 / 3),
                           [2],
                           [["1-90", "31-120", "61-150", "91-180"]])
        # The regions are 1-120, 31-150, and 61-180.
        # Regions 1-120 and 61-180 include module 0 (1-60) and module 2
        # (121-180), respectively, which each form 2 structures.
        # Region 31-150 overlaps both modules 0 and 2, so it can form
        # 2 x 2 = 4 structures.
        self.run_ensembles(relate_dirs, 120, 0.75,
                           [2, 4, 2],
                           [["1-120"],
                            ["31-150"],
                            ["61-180"]])
        # The only region is 1-180.
        # This region includes module 0 (1-60) and module 2 (121-180),
        # which each form 2 structures, so the entire RNA can form 4.
        self.run_ensembles(relate_dirs, 180, 0.5,
                           [4],
                           [["1-180"]])

    def test_modules012_read120(self):
        # Similar to test_modules012_read180 but with 120 nt reads.
        relate_dirs = self.sim_data([0, 1, 2], 120)
        self.run_ensembles(relate_dirs, 60, (2 / 3),
                           [2, 1, 2],
                           [["1-60", "21-80", "41-100"],
                            ["61-120"],
                            ["81-140", "101-160", "121-180"]])
        self.run_ensembles(relate_dirs, 90, (2 / 3),
                           [2],
                           [["1-90", "31-120", "61-150", "91-180"]])
        self.run_ensembles(relate_dirs, 120, 0.75,
                           [2, 4, 2],
                           [["1-120"],
                            ["31-150"],
                            ["61-180"]])
        self.run_ensembles(relate_dirs, 180, 0.5,
                           [4],
                           [["1-180"]])

    def test_modules012_read60(self):
        # Similar to test_modules012_read180 but with 60 nt reads.
        relate_dirs = self.sim_data([0, 1, 2], 60)
        self.run_ensembles(relate_dirs, 60, (2 / 3),
                           [2, 1, 2],
                           [["1-60", "21-80", "41-100"],
                            ["61-120"],
                            ["81-140", "101-160", "121-180"]])
        self.run_ensembles(relate_dirs, 90, (2 / 3),
                           [2],
                           [["1-90", "31-120", "61-150", "91-180"]])
        # With 60 nt reads, the reads are not long enough to provide
        # information on both modules 0 and 2 simultaneously, so the
        # algorithm cannot tell that together they form 4 clusters.
        self.run_ensembles(relate_dirs, 120, 0.75,
                           [2],
                           [["1-120", "31-150", "61-180"]])
        self.run_ensembles(relate_dirs, 180, 0.5,
                           [2],
                           [["1-180"]])

    def test_modules02_read60(self):
        relate_dirs = self.sim_data([0, 2], 60)
        # The regions are 1-60, 31-90, and 61-120.
        # Region 1-60 coincides exactly with module 0 and should form
        # 2 structures.
        # Region 31-90 overlaps module 0 (1-60) and module 2 (61-120),
        # which each form 2 structures, so this region can form 2 x 2
        # = 4 structures.
        # Region 61-120 coincides exactly with module 2 and should form
        # 2 structures.
        self.run_ensembles(relate_dirs, 60, 0.5,
                           [2, 4, 2],
                           [["1-60"],
                            ["31-90"],
                            ["61-120"]],
                           max_marcd_join=0.0)
        # Now rerun while limiting every region to at most 2 clusters.
        # The regions should not be joined together because region 31-90
        # with 2 clusters will not be sufficiently similar to 1-60 or to
        # 61-120 to be able to join with them.
        self.run_ensembles(relate_dirs, 60, 0.5,
                           [2, 2, 2],
                           [["1-60"],
                            ["31-90"],
                            ["61-120"]],
                           max_marcd_join=0.0,
                           max_clusters=2,
                           force=True)
        # Now rerun while tolerating larger differences between clusters
        # and confirm that all three regions are now joined.
        self.run_ensembles(relate_dirs, 60, 0.5,
                           [2],
                           [["1-60", "31-90", "61-120"]],
                           max_marcd_join=1.0,
                           max_clusters=2,
                           force=True)


if __name__ == "__main__":
    ut.main(verbosity=2)
