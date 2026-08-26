import tempfile
import unittest as ut
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from seismicrna.core.depend import dependency_exists
from seismicrna.core.header import make_header
from seismicrna.core.logs import Level, set_config
from seismicrna.core.rna.io import from_ct
from seismicrna.core.rna.state import RNAState
from seismicrna.core.rna.struct import RNAStructure
from seismicrna.core.seq.region import Region, SEQ_INDEX_NAMES
from seismicrna.core.seq.xna import DNA
from seismicrna.core.shell import VIENNA_RNACOFOLD_CMD
from seismicrna.core.table.base import COVER_REL, INFOR_REL, MUTAT_REL, TABLE_RELS
from seismicrna.filter.main import run as run_filter
from seismicrna.fold.main import run as run_fold
from seismicrna.idmut.io import IDmutBatchIO, ReadNamesBatchIO, RefseqIO
from seismicrna.idmut.report import IDmutReport
from seismicrna.duplex.dataset import LINKER, LINKER_LENGTH
from seismicrna.core.error import IncompatibleOptionsError, IncompatibleValuesError
from seismicrna.duplex.main import (
    Strand,
    iter_duplex_pairs,
    make_duplex,
    run as run_duplex,
    strand_uses_table_region,
)
from seismicrna.duplex.table import DuplexPositionTable, DuplexPositionTableLoader
from seismicrna.duplex.profile import RNACofoldProfile
from seismicrna.fold.main import resolve_duplex_energy_method


class DuplexEnergyMethodTest(ut.TestCase):
    """Cofolding uses SHAPE flags directly (Deigan), DMS pseudo-energies
    (Cordero), and rejects any method RNAcofold cannot express."""

    def _profile(self, method: str) -> RNACofoldProfile:
        seq_a, seq_b = DNA("GGGAAACCCAAAGGGAAACCC"), DNA("GGGTTTCCCAAAGGGTTTCCC")
        fused = seq_a + DNA(LINKER) + seq_b
        region = Region("A__B", fused)
        mus = pd.Series(np.linspace(0.01, 0.2, len(fused)), index=region.range)
        return RNACofoldProfile(
            region=region,
            sample="s",
            branches=dict(),
            mus_reg="full",
            mus_name="x",
            mus=mus,
            cut=len(seq_a),
            fold_temp=37.0,
            fold_energy_method=method,
            fold_quantile=0.95,
            deigan_slope=1.8,
            deigan_intercept=-0.6,
        )

    def test_resolve_auto(self):
        self.assertEqual(resolve_duplex_energy_method("DMS", "auto"), "Cordero")
        self.assertEqual(resolve_duplex_energy_method("SHAPE", "auto"), "Deigan")

    def test_resolve_passthrough(self):
        self.assertEqual(resolve_duplex_energy_method("DMS", "Deigan"), "Deigan")
        self.assertEqual(resolve_duplex_energy_method("DMS", "Cordero"), "Cordero")

    def test_resolve_rejects_unsupported(self):
        # Eddy and any other unimplemented strategy are rejected.
        for method in ("Eddy", "madeup"):
            self.assertRaises(
                IncompatibleOptionsError, resolve_duplex_energy_method, "DMS", method
            )

    def test_deigan_uses_shape_flags(self):
        """Deigan feeds reactivities with the user's slope/intercept."""
        p = self._profile("Deigan")
        self.assertIs(p._fold_data, p.mus_normalized)
        self.assertEqual(p.cofold_shape_method, "Dm1.8b-0.6")

    def test_cordero_uses_pseudomus(self):
        """Cordero (DMS) feeds pseudo-mutation rates, not reactivities."""
        p = self._profile("Cordero")
        self.assertIs(p._fold_data, p.pseudomus)
        self.assertTrue(p.cofold_shape_method.startswith("Dm"))

    def test_profile_rejects_unsupported(self):
        p = self._profile("Eddy")
        self.assertRaises(IncompatibleOptionsError, lambda: p.cofold_shape_method)
        self.assertRaises(IncompatibleOptionsError, lambda: p._fold_data)

    def test_cordero_fit_is_per_strand(self):
        """Each strand's Cordero pseudoenergies are fit independently, so
        one strand's reactivities do not affect the other's (a pooled fit
        would couple them through a shared scale factor)."""
        seq = DNA("ACACACACACACACAC")

        def five_prime_pe(hi3: float):
            fused = seq + DNA(LINKER) + seq
            region = Region("A__B", fused)
            lo = np.linspace(0.01, 0.05, len(seq))  # 5' strand (fixed)
            hi = np.full(len(seq), hi3)  # 3' strand (varied)
            linker = np.full(len(LINKER), np.nan)
            mus = pd.Series(np.concatenate([lo, linker, hi]), index=region.range)
            p = RNACofoldProfile(
                region=region,
                sample="s",
                branches=dict(),
                mus_reg="full",
                mus_name="x",
                mus=mus,
                cut=len(seq),
                fold_temp=37.0,
                fold_energy_method="Cordero",
                fold_quantile=0.95,
                deigan_slope=1.8,
                deigan_intercept=-0.6,
            )
            return p.pseudoenergies.iloc[: len(seq)].values

        np.testing.assert_allclose(five_prime_pe(0.10), five_prime_pe(0.40))


SAMPLE = "s"
N_READS = 6
NAMES = np.array([f"r{i}" for i in range(N_READS)])
SEQ_A = DNA("GGGAAACCCAAAAAGGGAAACCC")
SEQ_B = DNA("GGGTTTCCCAAAAAGGGTTTCCC")
MUTS_A = {5: {32: [0, 1, 2]}, 15: {32: [0, 3]}}
MUTS_B = {5: {32: [1, 4]}, 16: {32: [2, 3]}}


def write_idmut(out_dir: Path, ref: str, seq: DNA, muts: dict) -> Path:
    branches = dict()
    began = datetime.now()
    _, rc = RefseqIO(sample=SAMPLE, branches=branches, ref=ref, refseq=seq).save(
        out_dir
    )
    m = {
        pos: {
            rel: np.array(reads, dtype=int) for rel, reads in muts.get(pos, {}).items()
        }
        for pos in range(1, len(seq) + 1)
    }
    _, bc = IDmutBatchIO(
        sample=SAMPLE,
        branches=branches,
        region=Region(ref, seq),
        batch=0,
        seg_end5s=np.array([[1]] * N_READS),
        seg_end3s=np.array([[len(seq)]] * N_READS),
        muts=m,
    ).save(out_dir)
    _, nc = ReadNamesBatchIO(
        sample=SAMPLE, branches=branches, ref=ref, batch=0, names=NAMES
    ).save(out_dir)
    return IDmutReport(
        sample=SAMPLE,
        branches=branches,
        ref=ref,
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
        checksums={ReadNamesBatchIO.btype(): [nc], IDmutBatchIO.btype(): [bc]},
        refseq_checksum=rc,
        began=began,
        ended=datetime.now(),
    ).save(out_dir)


def make_cluster_strand(ref: str, seq: DNA, best_k: int) -> Strand:
    """Build a Strand backed by a synthetic clustered position table
    with clusters 1..best_k, each carrying distinct mutation rates."""
    header = make_header(rels=TABLE_RELS, ks=[best_k])
    index = pd.MultiIndex.from_arrays(
        [list(range(1, len(seq) + 1)), list(str(seq))], names=SEQ_INDEX_NAMES
    )
    data = pd.DataFrame(0.0, index=index, columns=header.index)
    rng = np.random.default_rng(len(ref))
    for rel, k, clust in header.index:
        if rel in (COVER_REL, INFOR_REL):
            data[(rel, k, clust)] = 100.0
        elif rel == MUTAT_REL:
            # Distinct rates per cluster so the profiles differ.
            data[(rel, k, clust)] = rng.uniform(0.0, 20.0, len(seq)) * clust
        else:  # MATCH_REL
            data[(rel, k, clust)] = 90.0
    table = SimpleNamespace(
        region=Region(ref, seq),
        refseq=seq,
        ref=ref,
        reg="full",
        sample=SAMPLE,
        branches=dict(),
        data=data,
        header=header,
        _dataset=SimpleNamespace(probe="DMS", best_k=best_k),
    )
    return Strand(
        seq=seq,
        ref=ref,
        reg="full",
        sample=SAMPLE,
        branches=dict(),
        table=table,
        source="",
    )


class DuplexClusterProductTest(ut.TestCase):
    """A duplex of clustered sources is the cross-product of their clusters."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)
        set_config(verbosity=Level.ERROR, exit_on_error=True)

    def tearDown(self):
        self._tmp.cleanup()
        set_config()

    def _duplex(self, k1: int, k2: int) -> DuplexPositionTableLoader:
        s1 = make_cluster_strand("molA", SEQ_A, k1)
        s2 = make_cluster_strand("molB", SEQ_B, k2)
        make_duplex(s1, s2, self.out, branch="", force=True)
        csv = next(self.out.rglob("*duplex-position-table.csv"))
        return DuplexPositionTableLoader(csv)

    def test_propagates_parent_branch(self):
        """A parent (source) branch is carried into the duplex path."""
        s1 = make_cluster_strand("molA", SEQ_A, 1)
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        s1.branches = {"filter": "myfork"}
        s2.branches = {"filter": "myfork"}
        make_duplex(s1, s2, self.out, branch="", force=True)
        csv = next(self.out.rglob("*duplex-position-table.csv"))
        self.assertIn("duplex_myfork", str(csv))

    def test_rejects_mismatched_parent_branches(self):
        """Two data sources on different branches cannot be duplexed."""
        s1 = make_cluster_strand("molA", SEQ_A, 1)
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        s1.branches = {"filter": "a"}
        s2.branches = {"filter": "b"}
        self.assertRaises(
            IncompatibleOptionsError,
            make_duplex,
            s1,
            s2,
            self.out,
            branch="",
            force=True,
        )

    def test_rejects_mismatched_probes(self):
        """Two data sources probed with different chemicals cannot be
        duplexed (the duplex folds with one energy method)."""
        s1 = make_cluster_strand("molA", SEQ_A, 1)
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        s2.table._dataset.probe = "SHAPE"
        self.assertRaises(
            IncompatibleOptionsError,
            make_duplex,
            s1,
            s2,
            self.out,
            branch="",
            force=True,
        )

    def test_cross_product_count(self):
        """K1 x K2 clusters give K1*K2 fused profiles over the duplex."""
        t = self._duplex(2, 3)
        self.assertEqual(t._dataset.best_k, 6)
        profs = list(t.iter_profiles(fold_table_region=True))
        self.assertEqual(len(profs), 6)
        length = len(SEQ_A) + LINKER_LENGTH + len(SEQ_B)
        for _, prof in profs:
            self.assertEqual(prof.region.length, length)

    def test_cross_product_maps_source_clusters(self):
        """Combination n=(i,j) carries strand-1 cluster i on the 5' side
        and strand-2 cluster j on the 3' side (product order)."""
        s1 = make_cluster_strand("molA", SEQ_A, 2)
        s2 = make_cluster_strand("molB", SEQ_B, 2)
        make_duplex(s1, s2, self.out, branch="", force=True)
        t = DuplexPositionTableLoader(
            next(self.out.rglob("*duplex-position-table.csv"))
        )
        cut = len(SEQ_A)
        # Expected 5' mutation rates for each source-1 cluster.
        exp1 = {
            i: s1.table.data[(MUTAT_REL, 2, i)].values
            / s1.table.data[(INFOR_REL, 2, i)].values
            for i in (1, 2)
        }
        # product(clusters1, clusters2) -> combos 1..4 map to
        # (i,j) = (1,1),(1,2),(2,1),(2,2).
        expected_i = {1: 1, 2: 1, 3: 2, 4: 2}
        for n, i in expected_i.items():
            mus = t.fetch_ratio(k=4, clust=n, rel=MUTAT_REL, squeeze=True)
            # The stored counts are rounded to PRECISION=1 decimal, so
            # the recovered ratio matches only to ~5e-4.
            np.testing.assert_allclose(mus.values[:cut], exp1[i], rtol=0, atol=2e-3)

    def test_average_source_stays_unclustered(self):
        """Combining two unclustered strands yields one average profile."""
        s1 = make_cluster_strand("molA", SEQ_A, 1)
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        make_duplex(s1, s2, self.out, branch="", force=True)
        t = DuplexPositionTableLoader(
            next(self.out.rglob("*duplex-position-table.csv"))
        )
        self.assertEqual(t._dataset.duplex_ks, [])
        self.assertFalse(t.header.get_is_clustered())
        self.assertEqual(len(list(t.iter_profiles(fold_table_region=True))), 1)

    @ut.skipUnless(dependency_exists(VIENNA_RNACOFOLD_CMD), "RNAcofold not installed")
    def test_fallback_table_without_a_best_k(self):
        """A clustered table whose cluster step found no best number of
        clusters (best_k = 0, its fallback) is duplexed using the one
        number of clusters it contains, rather than producing none."""
        s1 = make_cluster_strand("molA", SEQ_A, 1)
        # The cluster step records best_k = 0 when no K passed its filters.
        s1.table._dataset.best_k = 0
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        s2.table._dataset.best_k = 0
        make_duplex(s1, s2, self.out, branch="", force=True)
        table = DuplexPositionTableLoader(
            next(self.out.rglob("*duplex-position-table.csv"))
        )
        # One profile, as for any 1 x 1 duplex; not zero.
        self.assertEqual(len(list(table.iter_profiles(fold_table_region=True))), 1)

    def test_no_best_k_and_several_ks_is_an_error(self):
        """If several numbers of clusters are available and none is best,
        which to duplex is ambiguous, so it raises rather than guessing."""
        s1 = make_cluster_strand("molA", SEQ_A, 2)
        s1.table._dataset.best_k = 0
        s1.table.header = make_header(rels=TABLE_RELS, ks=[1, 2])
        s2 = make_cluster_strand("molB", SEQ_B, 1)
        self.assertRaisesRegex(
            IncompatibleValuesError,
            "no best number of clusters",
            make_duplex,
            s1,
            s2,
            self.out,
            branch="",
            force=True,
        )

    def test_cross_product_folds(self):
        """Each fused cluster profile folds into its own duplex."""
        s1 = make_cluster_strand("molA", SEQ_A, 2)
        s2 = make_cluster_strand("molB", SEQ_B, 2)
        report = make_duplex(s1, s2, self.out, branch="", force=True)
        reports = run_fold(
            [str(report.parent)],
            branch="",
            fold_coords=(),
            fold_primers=(),
            fold_regions_file=None,
            fold_table_region=True,
            fold_dry_run=False,
            fold_backend="ViennaRNA",
            pseudoknots=False,
            fold_energy_method="Deigan",
            deigan_slope=1.8,
            deigan_intercept=-0.6,
            fold_temp=37.0,
            fold_quantile=0.95,
            fold_fpaired=0.5,
            fold_constraint=None,
            fold_commands=None,
            eddy_prior_paired_file=None,
            eddy_prior_unpaired_file=None,
            fold_md=0,
            fold_mfe=True,
            fold_max=1,
            fold_percent=0.0,
            fold_edelta=1.0,
            fold_isolated=False,
            tmp_pfx=str(self.out / "tmp"),
            keep_tmp=False,
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        # One CT per fused cluster (2 x 2 = 4).
        cts = list(Path(reports[0]).parent.glob("*.ct"))
        self.assertEqual(len(cts), 4)

    def test_cluster_product_with_no_data_partner(self):
        """A clustered strand duplexed with a data-less partner gives one
        fused profile per cluster of the clustered strand, each of which
        cofolds on its own reactivities."""
        partner_seq = DNA("GGGTTTCCC")
        s1 = make_cluster_strand("molA", SEQ_A, 3)
        s2 = Strand.from_sequence(partner_seq, "ASO")
        report = make_duplex(s1, s2, self.out, branch="", force=True)
        table = DuplexPositionTableLoader(
            next(self.out.rglob("*duplex-position-table.csv"))
        )
        self.assertEqual(table.ref, "molA__ASO")
        self.assertEqual(table.cut, len(SEQ_A))
        self.assertEqual(
            table.region.length, len(SEQ_A) + LINKER_LENGTH + len(partner_seq)
        )
        # Crossing 3 clusters with a single data-less strand gives 3.
        self.assertEqual(table._dataset.duplex_ks, [3])
        self.assertEqual(len(list(table.iter_profiles(fold_table_region=True))), 3)
        # The partner carries no data; the clustered strand does.
        self.assertTrue(
            table.data.iloc[len(SEQ_A) + LINKER_LENGTH :].isna().all().all()
        )
        self.assertTrue(table.data.iloc[: len(SEQ_A)].notna().any().any())
        if not dependency_exists(VIENNA_RNACOFOLD_CMD):
            self.skipTest("RNAcofold not installed")
        reports = run_fold(
            [str(report.parent)],
            branch="",
            fold_coords=(),
            fold_primers=(),
            fold_regions_file=None,
            fold_table_region=True,
            fold_dry_run=False,
            fold_backend="ViennaRNA",
            pseudoknots=False,
            fold_energy_method="Deigan",
            deigan_slope=1.8,
            deigan_intercept=-0.6,
            fold_temp=37.0,
            fold_quantile=0.95,
            fold_fpaired=0.5,
            fold_constraint=None,
            fold_commands=None,
            eddy_prior_paired_file=None,
            eddy_prior_unpaired_file=None,
            fold_md=0,
            fold_mfe=True,
            fold_max=1,
            fold_percent=0.0,
            fold_edelta=1.0,
            fold_isolated=False,
            tmp_pfx=str(self.out / "tmp"),
            keep_tmp=False,
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        # One CT per cluster, each named for its cluster of the duplex.
        cts = sorted(p.name for p in Path(reports[0]).parent.glob("*.ct"))
        self.assertEqual(cts, [f"full__cluster-3-{clust}.ct" for clust in range(1, 4)])


class DuplexTest(ut.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)
        set_config(verbosity=Level.ERROR, exit_on_error=True)
        ra = write_idmut(self.out, "molA", SEQ_A, MUTS_A)
        rb = write_idmut(self.out, "molB", SEQ_B, MUTS_B)
        self.dms_a = run_filter([ra], probe="DMS", branch="", filter_pos_table=True)
        self.dms_b = run_filter([rb], probe="DMS", branch="", filter_pos_table=True)
        run_duplex(
            list(self.dms_a) + list(self.dms_b),
            branch="",
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        self.duplex_csv = next(self.out.rglob("*duplex-position-table.csv"))

    def tearDown(self):
        self._tmp.cleanup()
        set_config()

    def test_loads_via_api(self):
        """A duplex table loads like any position table and reconstructs the
        fused identity, the strand break, and its two sources."""
        t = DuplexPositionTableLoader(self.duplex_csv)
        self.assertIsInstance(t, DuplexPositionTable)
        self.assertEqual(t.ref, "molA__molB")
        self.assertEqual(t.cut, len(SEQ_A))
        self.assertEqual(t.region.length, len(SEQ_A) + LINKER_LENGTH + len(SEQ_B))
        # refseq comes from the brickle and matches the CSV sequence.
        self.assertEqual(str(t.refseq), "".join(t.data.index.get_level_values("Base")))
        # The two sources are reachable, with read-level data.
        self.assertEqual(t.table1.ref, "molA")
        self.assertEqual(t.table2.ref, "molB")
        self.assertEqual(t.data1.shape, t.table1.data.shape)
        self.assertGreaterEqual(t.table1._dataset.num_batches, 1)
        # One fused profile spans the whole duplex.
        profs = list(t.iter_profiles(fold_table_region=True))
        self.assertEqual(len(profs), 1)
        self.assertEqual(profs[0][1].region.length, t.region.length)

    def test_duplex_sequence_partner(self):
        """A data-less partner sequence becomes the 3' strand: the duplex
        loads, its cut is the 5' strand's length, and the partner's
        positions carry no data."""
        partner = DNA("GGGTTTCCC")
        run_duplex(
            list(self.dms_a),
            duplex_sequence=(("partner", partner),),
            branch="",
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        csv = next(
            p
            for p in self.out.rglob("*duplex-position-table.csv")
            if "molA__partner" in str(p)
        )
        t = DuplexPositionTableLoader(csv)
        self.assertEqual(t.ref, "molA__partner")
        self.assertEqual(t.cut, len(SEQ_A))
        self.assertEqual(t.region.length, len(SEQ_A) + LINKER_LENGTH + len(partner))
        # The 5' strand has a source; the 3' strand has none.
        self.assertEqual(t.table1.ref, "molA")
        self.assertIsNone(t.table2)
        self.assertIsNone(t.data2)
        # The partner's rows (after the linker) are empty.
        self.assertTrue(t.data.iloc[len(SEQ_A) + LINKER_LENGTH :].isna().all().all())

    def test_duplex_named_sequence_partner(self):
        """A (name, sequence) partner names the 3' strand."""
        run_duplex(
            list(self.dms_a),
            duplex_sequence=(("ASO", DNA("GGGTTTCCC")),),
            branch="",
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        csv = next(
            p
            for p in self.out.rglob("*duplex-position-table.csv")
            if "molA__ASO" in str(p)
        )
        t = DuplexPositionTableLoader(csv)
        self.assertEqual(t.ref, "molA__ASO")
        self.assertEqual(t.table1.ref, "molA")
        self.assertIsNone(t.table2)

    def test_per_reference_region_override(self):
        """Per-reference overrides beat the global default in both
        directions, and a conflict is rejected."""
        r = strand_uses_table_region
        # Global default applies when a reference is not overridden.
        self.assertFalse(r("A", default=False, region_refs=set(), full_refs=set()))
        self.assertTrue(r("A", default=True, region_refs=set(), full_refs=set()))
        # An override flips that reference against the default.
        self.assertTrue(r("A", default=False, region_refs={"A"}, full_refs=set()))
        self.assertFalse(r("A", default=True, region_refs=set(), full_refs={"A"}))
        # One reference full, another region, under the same default.
        self.assertFalse(r("A", default=False, region_refs={"B"}, full_refs=set()))
        self.assertTrue(r("B", default=False, region_refs={"B"}, full_refs=set()))
        # Naming a reference in both lists is an error.
        self.assertRaises(
            IncompatibleOptionsError,
            r,
            "A",
            default=False,
            region_refs={"A"},
            full_refs={"A"},
        )

    def test_duplex_dimer(self):
        """--dimer combines a table with itself into a homodimer."""
        run_duplex(
            list(self.dms_a),
            dimer=True,
            branch="",
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        csv = next(
            p
            for p in self.out.rglob("*duplex-position-table.csv")
            if "molA__molA" in str(p)
        )
        t = DuplexPositionTableLoader(csv)
        # A homodimer's reference stays joined (molA__molA) so its fold is
        # never confused with a monomer fold of molA.
        self.assertEqual(t.ref, "molA__molA")
        self.assertEqual(t.cut, len(SEQ_A))
        self.assertEqual(t.region.length, 2 * len(SEQ_A) + LINKER_LENGTH)
        self.assertEqual(t.table1.ref, "molA")
        self.assertEqual(t.table2.ref, "molA")
        self.assertTrue(np.allclose(t.data1.values, t.data2.values, equal_nan=True))

    @ut.skipUnless(dependency_exists(VIENNA_RNACOFOLD_CMD), "RNAcofold not installed")
    def test_fold_duplex_is_graphable(self):
        """Folding a duplex produces a duplex structure (via RNAcofold) that
        pairs with the duplex profile as one entity."""
        reports = run_fold(
            [str(self.duplex_csv.parent)],
            branch="",
            fold_coords=(),
            fold_primers=(),
            fold_regions_file=None,
            fold_table_region=True,
            fold_dry_run=False,
            fold_backend="ViennaRNA",
            pseudoknots=False,
            fold_energy_method="Deigan",
            deigan_slope=1.8,
            deigan_intercept=-0.6,
            fold_temp=37.0,
            fold_quantile=0.95,
            fold_fpaired=0.5,
            fold_constraint=None,
            fold_commands=None,
            eddy_prior_paired_file=None,
            eddy_prior_unpaired_file=None,
            fold_md=0,
            fold_mfe=True,
            fold_max=1,
            fold_percent=0.0,
            fold_edelta=1.0,
            fold_isolated=False,
            tmp_pfx=str(self.out / "tmp"),
            keep_tmp=False,
            verify_times=False,
            num_cpus=1,
            force=True,
        )
        self.assertGreaterEqual(len(reports), 1)
        ct = next(Path(reports[0]).parent.glob("*.ct"))
        struct = next(from_ct(ct))
        self.assertEqual(struct.region.length, len(SEQ_A) + LINKER_LENGTH + len(SEQ_B))
        # The duplex table's profile pairs with the duplex structure directly
        # (region matches -> one combined entity, no per-strand alignment).
        t = DuplexPositionTableLoader(self.duplex_csv)
        _, prof = next(t.iter_profiles(fold_table_region=True))
        state = RNAState.from_struct_profile(
            RNAStructure.from_db_string(
                struct.db_string,
                prof.region.seq,
                seq5=prof.region.end5,
                ref=prof.region.ref,
                reg=prof.region.name,
                title="t",
            ),
            prof,
        )
        self.assertEqual(state.region.length, t.region.length)


class IterDuplexPairsTest(ut.TestCase):
    """Duplex candidates are grouped by branches (what make_duplex
    requires of two strands) and paired in a deterministic order."""

    @staticmethod
    def _table(ref: str, branches: dict | None = None):
        # iter_duplex_pairs reads only these four attributes.
        return SimpleNamespace(
            ref=ref, reg="full", sample=SAMPLE, branches=branches or dict()
        )

    def test_every_pair(self):
        a, b, c = map(self._table, ["molA", "molB", "molC"])
        pairs = list(iter_duplex_pairs([a, b, c]))
        self.assertEqual(
            [(t1.ref, t2.ref) for t1, t2 in pairs],
            [("molA", "molB"), ("molA", "molC"), ("molB", "molC")],
        )

    def test_order_independent(self):
        a, b = self._table("molA"), self._table("molB")
        forward = [(t1.ref, t2.ref) for t1, t2 in iter_duplex_pairs([a, b])]
        reverse = [(t1.ref, t2.ref) for t1, t2 in iter_duplex_pairs([b, a])]
        self.assertEqual(forward, [("molA", "molB")])
        self.assertEqual(forward, reverse)

    def test_groups_by_branches(self):
        # Two tables on different branch chains never pair, since
        # make_duplex would reject them; two on the same chain do.
        a = self._table("molA")
        b = self._table("molB", dict(filter="alt"))
        self.assertEqual(list(iter_duplex_pairs([a, b])), [])
        c = self._table("molC", dict(filter="alt"))
        self.assertEqual(
            [(t1.ref, t2.ref) for t1, t2 in iter_duplex_pairs([a, b, c])],
            [("molB", "molC")],
        )

    def test_groups_by_step_not_flattened_branches(self):
        """Tables from different steps do not pair, even on default
        branches: flattening the branches would drop every step and every
        empty branch name, making a filter table and a duplex table look
        alike, but make_duplex compares the branches themselves."""
        filtered = self._table("molA", {"filter": ""})
        duplexed = self._table("molB__molC", {"filter": "", "duplex": ""})
        self.assertEqual(list(iter_duplex_pairs([filtered, duplexed])), [])

    def test_one_table_makes_no_pair(self):
        self.assertEqual(list(iter_duplex_pairs([self._table("molA")])), [])


class DuplexPairingTest(ut.TestCase):
    """Every way of choosing the 3' strand contributes its own duplexes,
    so they compose rather than override one another."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)
        set_config(verbosity=Level.ERROR, exit_on_error=True)
        ra = write_idmut(self.out, "molA", SEQ_A, MUTS_A)
        rb = write_idmut(self.out, "molB", SEQ_B, MUTS_B)
        self.tables = list(
            run_filter([ra], probe="DMS", branch="", filter_pos_table=True)
        ) + list(run_filter([rb], probe="DMS", branch="", filter_pos_table=True))

    def tearDown(self):
        self._tmp.cleanup()
        set_config()

    def _duplex_refs(self, **kwargs):
        """Run duplex and return the reference of every duplex made."""
        run_duplex(
            self.tables, branch="", verify_times=False, num_cpus=1, force=True, **kwargs
        )
        return {
            DuplexPositionTableLoader(csv).ref
            for csv in self.out.rglob("*duplex-position-table.csv")
        }

    def test_pairwise_is_the_default(self):
        self.assertEqual(self._duplex_refs(), {"molA__molB"})

    def test_all_sources_compose(self):
        """Pairwise, --dimer, and a partner sequence all contribute at
        once; none of them suppresses the others."""
        self.assertEqual(
            self._duplex_refs(dimer=True, duplex_sequence=(("ASO", DNA("GGGTTTCCC")),)),
            {"molA__molB", "molA__molA", "molB__molB", "molA__ASO", "molB__ASO"},
        )

    def test_no_duplex_pair_leaves_partners(self):
        """--no-duplex-pair drops only the pairwise combinations."""
        self.assertEqual(
            self._duplex_refs(
                duplex_pair=False, duplex_sequence=(("ASO", DNA("GGGTTTCCC")),)
            ),
            {"molA__ASO", "molB__ASO"},
        )

    def test_nothing_to_duplex(self):
        """With every source of pairs turned off, nothing is made."""
        self.assertEqual(self._duplex_refs(duplex_pair=False), set())

    def test_rerunning_does_not_duplex_duplexes(self):
        """Running duplex again over a directory that already holds its
        own output makes the same duplexes, not chimeras of them: a
        duplex table among the inputs is not a strand of a new duplex."""
        first = self._duplex_refs()
        # The output directory now holds a duplex table as well as the
        # two filter tables, and is given back as the input.
        self.tables = [str(self.out)]
        self.assertEqual(self._duplex_refs(), first)
        self.assertEqual(first, {"molA__molB"})


if __name__ == "__main__":
    ut.main()
