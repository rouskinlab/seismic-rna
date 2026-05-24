import json
import tempfile
import unittest as ut
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd

from seismicrna.core import path
from seismicrna.core.error import NoDataError
from seismicrna.core.header import format_clust_name
from seismicrna.core.logs import Level, set_config
from seismicrna.core.rna import RNAStructure, to_ct
from seismicrna.core.seq.region import BASE_NAME, POS_NAME, Region
from seismicrna.core.seq.xna import DNA
from seismicrna.core.table import (
    COVER_REL,
    INFOR_REL,
    MUTAT_REL,
    SUBST_REL,
    SUB_A_REL,
    SUB_C_REL,
    SUB_G_REL,
    SUB_T_REL,
    DELET_REL,
    INSRT_REL,
)
from seismicrna.export.web import (
    META_SYMBOL,
    REF_SEQ,
    REG_END5,
    REG_END3,
    REG_POS,
    STRUCTURE,
    COVER_COUNT,
    INFOR_COUNT,
    SUBST_COUNT,
    SUB_A_COUNT,
    SUB_C_COUNT,
    SUB_G_COUNT,
    SUB_T_COUNT,
    DELET_COUNT,
    INSRT_COUNT,
)
from seismicrna.sim.abstract import (
    _calc_ratios,
    _calc_ratio_stats,
    _accumulate_ratios,
    _format_param_tokens,
    abstract_seismicgraph_file,
    abstract_table,
    new_parameter_dict,
)
from seismicrna.core.arg import (
    opt_pmut_paired,
    opt_pmut_unpaired,
    opt_vmut_paired,
    opt_vmut_unpaired,
)
from seismicrna.sim.muts import (
    _make_pmut_means_kwargs,
    _make_vmut_kwargs,
    make_pmut_means_paired,
    make_pmut_means_unpaired,
    make_vmut_paired,
    make_vmut_unpaired,
)


REL_COLUMNS = [
    COVER_REL,
    INFOR_REL,
    MUTAT_REL,
    SUBST_REL,
    SUB_A_REL,
    SUB_C_REL,
    SUB_G_REL,
    SUB_T_REL,
    DELET_REL,
    INSRT_REL,
]


def _build_counts_df(
    refseq: DNA, per_position: dict[int, dict[str, int]], cov: int = 10, info: int = 8
) -> pd.DataFrame:
    """Build a counts DataFrame indexed by (Position, Base).

    Every position starts with cov=cov, info=info and zero in every
    other column.  `per_position` overlays specific values per 1-indexed
    position (e.g. ``{2: {SUB_C_REL: 3, MUTAT_REL: 3, SUBST_REL: 3}}``).
    """
    length = len(refseq)
    positions = list(range(1, length + 1))
    index = pd.MultiIndex.from_tuples(
        [(pos, refseq[pos - 1]) for pos in positions], names=[POS_NAME, BASE_NAME]
    )
    data = {col: np.zeros(length, dtype=int) for col in REL_COLUMNS}
    data[COVER_REL][:] = cov
    data[INFOR_REL][:] = info
    df = pd.DataFrame(data, index=index)
    for pos, overrides in per_position.items():
        for col, value in overrides.items():
            df.at[(pos, refseq[pos - 1]), col] = value
    return df


def _write_ct(
    ct_path: Path,
    ref: str,
    refseq: DNA,
    reg_name: str,
    db_string: str,
    profile_name: str,
) -> Path:
    """Write a single-structure CT file to ``ct_path``."""
    struct = RNAStructure.from_db_string(
        db_string=db_string, seq=refseq, ref=ref, reg=reg_name, title=profile_name
    )
    ct_path.parent.mkdir(parents=True, exist_ok=True)
    to_ct([struct], ct_path, force=True)
    return ct_path


def _fold_ct_path(
    top: Path, sample: str, ref: str, reg: str, profile_name: str
) -> Path:
    """Build the canonical fold output CT path for a profile."""
    return top / sample / "fold" / ref / reg / f"{profile_name}.ct"


def _profile_name(reg: str, k: int, clust: int) -> str:
    """Construct the profile name fold/filter would use for a cluster."""
    mus_name = path.fill_whitespace(format_clust_name(k, clust), fill="-")
    return f"{reg}__{mus_name}"


class FakeHeader:
    def __init__(self, clusts):
        self.clusts = list(clusts)


class FakeFilterPositionTable:
    """Duck-typed FilterPositionTable for abstract_table tests.

    Only the attributes that ``abstract_table`` reads are implemented.
    """

    def __init__(
        self,
        *,
        top: Path,
        sample: str,
        ref: str,
        refseq: DNA,
        reg: str,
        clusts: list[tuple[int, int]],
        counts_per_clust: dict[tuple[int, int], pd.DataFrame],
    ):
        self.top = top
        self.sample = sample
        self.ref = ref
        self.region = Region(ref, refseq, name=reg)
        self.header = FakeHeader(clusts)
        self._counts = counts_per_clust

    @property
    def reg(self) -> str:
        return self.region.name

    def fetch_count(self, *, k: int, clust: int) -> pd.DataFrame:
        return self._counts[(k, clust)].copy()


def _write_seismicgraph_json(json_path: Path, sample: str, entries: list[dict]) -> Path:
    """Serialise a SEISMICgraph JSON file.

    Each ``entries`` item is a dict shaped like::

        {
            "ref": str,
            "refseq": DNA,
            "regions": [
                {
                    "name": str,
                    "end5": int,
                    "end3": int,
                    "positions": list[int],  # unmasked
                    "profiles": {
                        profile_name: {
                            COVER_COUNT: list[int],
                            INFOR_COUNT: list[int],
                            SUBST_COUNT: list[int],
                            SUB_A_COUNT/SUB_C/SUB_G/SUB_T: list[int],
                            DELET_COUNT: list[int],
                            INSRT_COUNT: list[int],
                            STRUCTURE: dot_bracket,
                        },
                        ...
                    }
                }, ...
            ]
        }
    """
    data = {f"{META_SYMBOL}sample": sample}
    for entry in entries:
        ref_block = {f"{META_SYMBOL}{REF_SEQ}": str(entry["refseq"])}
        for region in entry["regions"]:
            reg_block = {
                f"{META_SYMBOL}{REG_END5}": region["end5"],
                f"{META_SYMBOL}{REG_END3}": region["end3"],
                f"{META_SYMBOL}{REG_POS}": list(region["positions"]),
            }
            for profile_name, profile_data in region["profiles"].items():
                reg_block[profile_name] = profile_data
            ref_block[region["name"]] = reg_block
        data[entry["ref"]] = ref_block
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(data, f)
    return json_path


def _profile_data(
    length: int,
    *,
    cov: int = 10,
    info: int = 8,
    per_position: dict[int, dict[str, int]] | None = None,
    structure: str = "",
) -> dict:
    """Build a per-profile count dict for SEISMICgraph JSON."""
    profile: dict = {
        COVER_COUNT: [cov] * length,
        INFOR_COUNT: [info] * length,
        SUBST_COUNT: [0] * length,
        SUB_A_COUNT: [0] * length,
        SUB_C_COUNT: [0] * length,
        SUB_G_COUNT: [0] * length,
        SUB_T_COUNT: [0] * length,
        DELET_COUNT: [0] * length,
        INSRT_COUNT: [0] * length,
        STRUCTURE: structure,
    }
    if per_position:
        for pos_1, overrides in per_position.items():
            idx = pos_1 - 1
            for key, value in overrides.items():
                profile[key][idx] = value
    return profile


class AbstractTestBase(ut.TestCase):
    """Silence logging during tests; subclasses get tmpdirs as needed."""

    @classmethod
    def setUpClass(cls):
        # exit_on_error=False keeps logger.error from re-raising; the code
        # under test relies on find_files_chain swallowing FS errors via
        # logger.error so it can return an empty result.
        set_config(verbosity=Level.FATAL, exit_on_error=False)

    @classmethod
    def tearDownClass(cls):
        set_config()


# ---------------------------------------------------------------------------
# _calc_ratios
# ---------------------------------------------------------------------------


class TestCalcRatios(AbstractTestBase):
    """Verify _calc_ratios computes per-position ratios analytically."""

    REFSEQ = DNA("ACGTACGT")

    def _is_paired_mask(self, paired_positions: list[int]) -> np.ndarray:
        positions = range(1, len(self.REFSEQ) + 1)
        return np.array([pos in set(paired_positions) for pos in positions])

    def test_loq_paired_vs_unpaired(self):
        # cov=10 for every position; info=5 at positions 1,3,5,7 and 8 at 2,4,6,8.
        counts = _build_counts_df(self.REFSEQ, per_position={})
        for pos in [1, 3, 5, 7]:
            counts.at[(pos, self.REFSEQ[pos - 1]), INFOR_REL] = 5
        # Pair every other position (1,3,5,7 are paired).
        is_paired = self._is_paired_mask([1, 3, 5, 7])
        paired, unpaired = _calc_ratios(counts, is_paired)
        # Every paired position has info=5, cov=10 → loq = 1 - 0.5 = 0.5
        np.testing.assert_array_almost_equal(
            paired["loq"], np.array([0.5, 0.5, 0.5, 0.5])
        )
        # Every unpaired position has info=8, cov=10 → loq = 1 - 0.8 = 0.2
        np.testing.assert_array_almost_equal(
            unpaired["loq"], np.array([0.2, 0.2, 0.2, 0.2])
        )

    def test_per_base_substitution_ratios(self):
        # A positions are at 1 and 5 in "ACGTACGT".
        # Position 1 (A, paired): mutat=4, sub_C=2, sub_G=1, sub_T=1
        # Position 5 (A, paired): mutat=4, sub_C=0, sub_G=2, sub_T=2
        per_position = {
            1: {MUTAT_REL: 4, SUB_C_REL: 2, SUB_G_REL: 1, SUB_T_REL: 1},
            5: {MUTAT_REL: 4, SUB_C_REL: 0, SUB_G_REL: 2, SUB_T_REL: 2},
        }
        counts = _build_counts_df(self.REFSEQ, per_position=per_position)
        is_paired = self._is_paired_mask([1, 2, 3, 4, 5, 6, 7, 8])
        paired, _ = _calc_ratios(counts, is_paired)
        # Expected per-A ratios: paired['ac'], paired['ag'], paired['at']
        # For positions 1 and 5 (both paired):
        np.testing.assert_array_almost_equal(paired["ac"], np.array([2 / 4, 0 / 4]))
        np.testing.assert_array_almost_equal(paired["ag"], np.array([1 / 4, 2 / 4]))
        np.testing.assert_array_almost_equal(paired["at"], np.array([1 / 4, 2 / 4]))

    def test_per_base_mutation_rate_xm(self):
        # C positions are at 2 and 6 in "ACGTACGT".
        per_position = {
            2: {MUTAT_REL: 2},  # info default = 8
            6: {MUTAT_REL: 4},
        }
        counts = _build_counts_df(self.REFSEQ, per_position=per_position)
        is_paired = self._is_paired_mask([2, 6])
        paired, unpaired = _calc_ratios(counts, is_paired)
        # cm = MUTAT / INFOR for C positions in the paired mask.
        np.testing.assert_array_almost_equal(paired["cm"], np.array([2 / 8, 4 / 8]))
        # Unpaired had no C positions → empty array.
        self.assertEqual(unpaired["cm"].size, 0)

    def test_unpaired_collects_only_unpaired_positions(self):
        # G positions are at 3 and 7.  Make 3 paired and 7 unpaired.
        per_position = {3: {MUTAT_REL: 1}, 7: {MUTAT_REL: 5}}
        counts = _build_counts_df(self.REFSEQ, per_position=per_position)
        is_paired = self._is_paired_mask([3])
        paired, unpaired = _calc_ratios(counts, is_paired)
        np.testing.assert_array_almost_equal(paired["gm"], np.array([1 / 8]))
        np.testing.assert_array_almost_equal(unpaired["gm"], np.array([5 / 8]))


# ---------------------------------------------------------------------------
# _calc_ratio_stats
# ---------------------------------------------------------------------------


class TestCalcRatioStats(AbstractTestBase):
    MARGIN = 1.0e-6

    def _full_ratios(self, **overrides) -> dict:
        ratios = new_parameter_dict()
        ratios.update(overrides)
        return ratios

    def test_mean_clipped_to_margin_when_zero(self):
        ratios = self._full_ratios(am=np.array([0.0, 0.0]))
        means, _ = _calc_ratio_stats(ratios, margin=self.MARGIN)
        self.assertEqual(means["am"], self.MARGIN)

    def test_mean_clipped_to_upper_bound(self):
        # Values very near 1.0 should clip to 1 - margin.
        ratios = self._full_ratios(cm=np.array([1.0, 1.0, 1.0]))
        means, _ = _calc_ratio_stats(ratios, margin=self.MARGIN)
        self.assertAlmostEqual(means["cm"], 1.0 - self.MARGIN)

    def test_fvar_known_values(self):
        values = np.array([0.01, 0.02, 0.03])
        ratios = self._full_ratios(am=values)
        _, fvar = _calc_ratio_stats(ratios, margin=self.MARGIN)
        expected_mean = float(np.nanmean(values))
        expected_var = float(np.nanvar(values))
        expected_fvar = expected_var / (expected_mean * (1.0 - expected_mean))
        self.assertAlmostEqual(fvar["a"], expected_fvar, places=10)

    def test_fvar_falls_back_to_margin_when_too_few_values(self):
        ratios = self._full_ratios(am=np.array([0.01]))
        _, fvar = _calc_ratio_stats(ratios, margin=self.MARGIN)
        self.assertEqual(fvar["a"], self.MARGIN)

    def test_fvar_omits_n_base(self):
        ratios = self._full_ratios(nm=np.array([0.1, 0.2, 0.3]))
        _, fvar = _calc_ratio_stats(ratios, margin=self.MARGIN)
        self.assertNotIn("n", fvar)

    def test_fvar_keys_are_acgt(self):
        ratios = self._full_ratios(
            am=np.array([0.01, 0.02]),
            cm=np.array([0.01, 0.02]),
            gm=np.array([0.01, 0.02]),
            tm=np.array([0.01, 0.02]),
        )
        _, fvar = _calc_ratio_stats(ratios, margin=self.MARGIN)
        self.assertEqual(set(fvar), {"a", "c", "g", "t"})


# ---------------------------------------------------------------------------
# _accumulate_ratios
# ---------------------------------------------------------------------------


class TestAccumulateRatios(AbstractTestBase):
    def _pair_with_am(self, values: np.ndarray) -> tuple[dict, dict]:
        paired = new_parameter_dict()
        unpaired = new_parameter_dict()
        paired["am"] = values
        unpaired["am"] = values * 2.0
        return paired, unpaired

    def test_concatenates_arrays_per_key(self):
        p1 = self._pair_with_am(np.array([0.1, 0.2]))
        p2 = self._pair_with_am(np.array([0.3]))
        paired, unpaired = _accumulate_ratios(iter([p1, p2]))
        np.testing.assert_array_almost_equal(paired["am"], np.array([0.1, 0.2, 0.3]))
        np.testing.assert_array_almost_equal(unpaired["am"], np.array([0.2, 0.4, 0.6]))

    def test_accumulates_from_empty(self):
        paired, unpaired = _accumulate_ratios(iter([]))
        baseline = new_parameter_dict()
        self.assertEqual(set(paired), set(baseline))
        for key in baseline:
            self.assertEqual(paired[key].size, 0)
            self.assertEqual(unpaired[key].size, 0)

    def test_three_inputs_compound(self):
        a = self._pair_with_am(np.array([0.1]))
        b = self._pair_with_am(np.array([0.2, 0.3]))
        c = self._pair_with_am(np.array([0.4, 0.5, 0.6]))
        paired, _ = _accumulate_ratios(iter([a, b, c]))
        self.assertEqual(paired["am"].size, 6)
        np.testing.assert_array_almost_equal(
            paired["am"], np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        )


# ---------------------------------------------------------------------------
# abstract_seismicgraph_file
# ---------------------------------------------------------------------------


class TestAbstractSeismicgraphFile(AbstractTestBase):
    REF = "ref"
    REFSEQ = DNA("ACGTACGT")
    REG = "ACGTACGT"  # default hyphenless name when end5=1, end3=8 → 'full'
    # Use an explicit name to avoid Region name defaulting.

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _profile_full_mutated_at(
        self,
        positions: list[int],
        sub_col: str = SUB_C_COUNT,
        *,
        structure: str | None = None,
    ) -> dict:
        """Profile where each ``positions`` entry contributes one C-sub."""
        length = len(self.REFSEQ)
        if structure is None:
            structure = "." * length
        per_position = {p: {sub_col: 1, SUBST_COUNT: 1} for p in positions}
        return _profile_data(length, per_position=per_position, structure=structure)

    def test_single_profile_single_region(self):
        # Put a substitution at A positions 1 (paired) and 5 (unpaired).
        # SUBST → MUTAT (via the function's MUTAT = SUBST+DELET+INSRT line).
        structure = "((....))"  # positions 1,2,7,8 paired; 3-6 unpaired
        profile = self._profile_full_mutated_at([1, 5], structure=structure)
        json_path = _write_seismicgraph_json(
            self.tmpdir / "data.json",
            sample="s",
            entries=[
                {
                    "ref": self.REF,
                    "refseq": self.REFSEQ,
                    "regions": [
                        {
                            "name": "reg1",
                            "end5": 1,
                            "end3": len(self.REFSEQ),
                            "positions": list(range(1, len(self.REFSEQ) + 1)),
                            "profiles": {"reg1__average": profile},
                        }
                    ],
                }
            ],
        )
        paired, unpaired = abstract_seismicgraph_file(json_path)
        # Position 1 is A and paired; MUTAT=1, INFOR=8 → am paired = 1/8.
        np.testing.assert_array_almost_equal(paired["am"], np.array([1 / 8]))
        # Position 5 is A and unpaired; MUTAT=1, INFOR=8 → am unpaired = 1/8.
        np.testing.assert_array_almost_equal(unpaired["am"], np.array([1 / 8]))
        # ac for position 1: SUB_C=1, MUTAT=1 → 1.0
        np.testing.assert_array_almost_equal(paired["ac"], np.array([1.0]))

    def test_multiple_profiles_within_region(self):
        structure = "((....))"
        prof_a = self._profile_full_mutated_at([1], structure=structure)
        prof_b = self._profile_full_mutated_at([5], structure=structure)
        json_path = _write_seismicgraph_json(
            self.tmpdir / "data.json",
            sample="s",
            entries=[
                {
                    "ref": self.REF,
                    "refseq": self.REFSEQ,
                    "regions": [
                        {
                            "name": "reg1",
                            "end5": 1,
                            "end3": len(self.REFSEQ),
                            "positions": list(range(1, len(self.REFSEQ) + 1)),
                            "profiles": {
                                "reg1__cluster-2-1": prof_a,
                                "reg1__cluster-2-2": prof_b,
                            },
                        }
                    ],
                }
            ],
        )
        paired, unpaired = abstract_seismicgraph_file(json_path)
        # Profile A: position 1 (A, paired) substituted → paired am has 1/8.
        # Profile B: position 5 (A, unpaired) substituted → unpaired am has 1/8.
        # Profile A's position 5 has no substitution → unpaired am 0/8 = 0.
        # Profile B's position 1 has no substitution → paired am 0/8 = 0.
        # Each profile contributes both A positions to the paired/unpaired
        # buckets that match its structure.
        self.assertEqual(paired["am"].size, 2)
        self.assertEqual(unpaired["am"].size, 2)
        # Sums should match: paired am total = 1/8 + 0 from both profiles.
        self.assertAlmostEqual(float(paired["am"].sum()), 1 / 8)
        self.assertAlmostEqual(float(unpaired["am"].sum()), 1 / 8)

    def test_multiple_regions(self):
        structure_full = "." * len(self.REFSEQ)
        # Region 1 (full length): substitution at position 5 (A, unpaired).
        prof1 = self._profile_full_mutated_at([5], structure=structure_full)
        # Region 2 spans positions 1-4 only: substitution at position 2 (C).
        refseq_short_len = 4
        per_position = {2: {SUB_A_COUNT: 1, SUBST_COUNT: 1}}
        prof2 = _profile_data(
            refseq_short_len,
            per_position=per_position,
            structure="." * refseq_short_len,
        )
        json_path = _write_seismicgraph_json(
            self.tmpdir / "data.json",
            sample="s",
            entries=[
                {
                    "ref": self.REF,
                    "refseq": self.REFSEQ,
                    "regions": [
                        {
                            "name": "reg_full",
                            "end5": 1,
                            "end3": len(self.REFSEQ),
                            "positions": list(range(1, len(self.REFSEQ) + 1)),
                            "profiles": {"reg_full__average": prof1},
                        },
                        {
                            "name": "reg_short",
                            "end5": 1,
                            "end3": refseq_short_len,
                            "positions": list(range(1, refseq_short_len + 1)),
                            "profiles": {"reg_short__average": prof2},
                        },
                    ],
                }
            ],
        )
        paired, unpaired = abstract_seismicgraph_file(json_path)
        # reg_full unpaired A positions: 1 and 5 (2 entries); reg_short
        # unpaired A position: 1 (1 entry).  Total = 3.
        self.assertEqual(unpaired["am"].size, 3)
        # reg_short contributes a single C subbed to A: unpaired ca = 1.0
        np.testing.assert_array_almost_equal(unpaired["ca"], np.array([1.0]))

    def test_min_aucroc_skips_profile(self):
        # A profile whose mutations sit on a paired position has AUC < 0.5
        # for paired-detection (mutation rate higher when paired).
        structure = "((....))"  # 1,2,7,8 paired
        # Substitute every paired position only.
        prof = self._profile_full_mutated_at([1, 2, 7, 8], structure=structure)
        json_path = _write_seismicgraph_json(
            self.tmpdir / "data.json",
            sample="s",
            entries=[
                {
                    "ref": self.REF,
                    "refseq": self.REFSEQ,
                    "regions": [
                        {
                            "name": "reg1",
                            "end5": 1,
                            "end3": len(self.REFSEQ),
                            "positions": list(range(1, len(self.REFSEQ) + 1)),
                            "profiles": {"reg1__average": prof},
                        }
                    ],
                }
            ],
        )
        paired_skip, unpaired_skip = abstract_seismicgraph_file(
            json_path, min_aucroc=0.99
        )
        # All accumulators should be empty since the only profile is skipped.
        for key, values in paired_skip.items():
            self.assertEqual(values.size, 0, msg=f"paired[{key!r}] not empty: {values}")
        for key, values in unpaired_skip.items():
            self.assertEqual(
                values.size, 0, msg=f"unpaired[{key!r}] not empty: {values}"
            )

    def test_metadata_keys_ignored(self):
        # Extra '#'-prefixed keys at every level must not be treated as data.
        structure = "((....))"
        prof = self._profile_full_mutated_at([1], structure=structure)
        # Manually craft JSON because _write_seismicgraph_json doesn't emit
        # extra metadata keys.
        data = {
            f"{META_SYMBOL}sample": "s",
            f"{META_SYMBOL}custom": "ignore-me",
            self.REF: {
                f"{META_SYMBOL}{REF_SEQ}": str(self.REFSEQ),
                f"{META_SYMBOL}extra": "ignore-me",
                "reg1": {
                    f"{META_SYMBOL}{REG_END5}": 1,
                    f"{META_SYMBOL}{REG_END3}": len(self.REFSEQ),
                    f"{META_SYMBOL}{REG_POS}": list(range(1, len(self.REFSEQ) + 1)),
                    f"{META_SYMBOL}note": "ignore-me",
                    "reg1__average": prof,
                },
            },
        }
        json_path = self.tmpdir / "data.json"
        with open(json_path, "w") as f:
            json.dump(data, f)
        paired, _ = abstract_seismicgraph_file(json_path)
        # Should still process the single A position; not raise on '#' keys.
        np.testing.assert_array_almost_equal(paired["am"], np.array([1 / 8]))


# ---------------------------------------------------------------------------
# abstract_table
# ---------------------------------------------------------------------------


class TestAbstractTable(AbstractTestBase):
    REF = "ref"
    REFSEQ = DNA("ACGTACGT")
    REG = "myreg"
    SAMPLE = "sample"
    # Structure: position 1 paired, 5 unpaired (both A's).
    DB_STRING = "((....))"

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _make_table(
        self,
        clusts: list[tuple[int, int]],
        counts_per_clust: dict[tuple[int, int], pd.DataFrame],
    ) -> FakeFilterPositionTable:
        return FakeFilterPositionTable(
            top=self.tmpdir,
            sample=self.SAMPLE,
            ref=self.REF,
            refseq=self.REFSEQ,
            reg=self.REG,
            clusts=clusts,
            counts_per_clust=counts_per_clust,
        )

    def _seed_counts(
        self, mutated_paired: int = 1, mutated_unpaired: int = 1
    ) -> pd.DataFrame:
        # Position 1 is A & paired; position 5 is A & unpaired.
        per_position = {
            1: {
                MUTAT_REL: mutated_paired,
                SUB_C_REL: mutated_paired,
                SUBST_REL: mutated_paired,
            },
            5: {
                MUTAT_REL: mutated_unpaired,
                SUB_G_REL: mutated_unpaired,
                SUBST_REL: mutated_unpaired,
            },
        }
        return _build_counts_df(self.REFSEQ, per_position=per_position)

    def _write_ct(self, profile_name: str):
        ct_path = _fold_ct_path(
            self.tmpdir, self.SAMPLE, self.REF, self.REG, profile_name
        )
        return _write_ct(
            ct_path, self.REF, self.REFSEQ, self.REG, self.DB_STRING, profile_name
        )

    def test_non_clustered_auto_detect(self):
        counts = self._seed_counts()
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        self._write_ct(_profile_name(self.REG, 0, 0))  # "myreg__average"
        result = abstract_table(table, struct_file=None)
        self.assertIsNotNone(result)
        paired, unpaired = result
        # Paired am: position 1 has MUTAT=1, INFOR=8 → 1/8.
        np.testing.assert_array_almost_equal(paired["am"], np.array([1 / 8]))
        np.testing.assert_array_almost_equal(unpaired["am"], np.array([1 / 8]))
        # Paired ac at position 1: SUB_C=1 / MUTAT=1 = 1.0
        np.testing.assert_array_almost_equal(paired["ac"], np.array([1.0]))
        # Unpaired ag at position 5: SUB_G=1 / MUTAT=1 = 1.0
        np.testing.assert_array_almost_equal(unpaired["ag"], np.array([1.0]))

    def test_non_clustered_with_explicit_struct_file(self):
        counts = self._seed_counts(mutated_paired=2, mutated_unpaired=3)
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        ct_path = self._write_ct(_profile_name(self.REG, 0, 0))
        result = abstract_table(table, struct_file=[ct_path])
        self.assertIsNotNone(result)
        paired, unpaired = result
        np.testing.assert_array_almost_equal(paired["am"], np.array([2 / 8]))
        np.testing.assert_array_almost_equal(unpaired["am"], np.array([3 / 8]))

    def test_clustered_multiple_ct_files(self):
        counts_c1 = self._seed_counts(mutated_paired=1, mutated_unpaired=0)
        counts_c2 = self._seed_counts(mutated_paired=3, mutated_unpaired=2)
        table = self._make_table(
            clusts=[(2, 1), (2, 2)],
            counts_per_clust={(2, 1): counts_c1, (2, 2): counts_c2},
        )
        # Wrap fetch_count to mimic clustered multi-level columns.
        original_fetch = table.fetch_count

        def fetch_clustered(*, k, clust):
            df = original_fetch(k=k, clust=clust)
            df.columns = pd.MultiIndex.from_tuples(
                [(k, clust, col) for col in df.columns], names=["k", "clust", "rel"]
            )
            return df

        table.fetch_count = fetch_clustered  # type: ignore[assignment]
        self._write_ct(_profile_name(self.REG, 2, 1))
        self._write_ct(_profile_name(self.REG, 2, 2))
        result = abstract_table(table, struct_file=None)
        self.assertIsNotNone(result)
        paired, unpaired = result
        # Paired am should have one entry per cluster: 1/8 and 3/8.
        self.assertEqual(paired["am"].size, 2)
        self.assertEqual(set(paired["am"].tolist()), {1 / 8, 3 / 8})
        # Unpaired am: 0/8 and 2/8.
        self.assertEqual(unpaired["am"].size, 2)
        self.assertEqual(set(unpaired["am"].tolist()), {0.0, 2 / 8})

    def test_struct_file_ref_mismatch_raises(self):
        counts = self._seed_counts()
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        # Write a CT file under a different ref/reg directory tree.
        other_ct = _fold_ct_path(
            self.tmpdir,
            self.SAMPLE,
            "other_ref",
            self.REG,
            _profile_name(self.REG, 0, 0),
        )
        _write_ct(
            other_ct,
            "other_ref",
            self.REFSEQ,
            self.REG,
            self.DB_STRING,
            _profile_name(self.REG, 0, 0),
        )
        with self.assertRaises(ValueError) as cm:
            abstract_table(table, struct_file=[other_ct])
        self.assertIn("reference", str(cm.exception))

    def test_struct_file_reg_mismatch_raises(self):
        counts = self._seed_counts()
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        other_ct = _fold_ct_path(
            self.tmpdir, self.SAMPLE, self.REF, "other_reg", "other_reg__average"
        )
        _write_ct(
            other_ct,
            self.REF,
            self.REFSEQ,
            "other_reg",
            self.DB_STRING,
            "other_reg__average",
        )
        with self.assertRaises(ValueError) as cm:
            abstract_table(table, struct_file=[other_ct])
        self.assertIn("region", str(cm.exception))

    def test_no_ct_files_raises_no_data_error(self):
        counts = self._seed_counts()
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        # No CT file written.
        with self.assertRaises(NoDataError):
            abstract_table(table, struct_file=None)

    def test_missing_cluster_ct_file_uses_only_available(self):
        counts_c1 = self._seed_counts(mutated_paired=1, mutated_unpaired=0)
        counts_c2 = self._seed_counts(mutated_paired=3, mutated_unpaired=2)
        table = self._make_table(
            clusts=[(2, 1), (2, 2)],
            counts_per_clust={(2, 1): counts_c1, (2, 2): counts_c2},
        )
        original_fetch = table.fetch_count

        def fetch_clustered(*, k, clust):
            df = original_fetch(k=k, clust=clust)
            df.columns = pd.MultiIndex.from_tuples(
                [(k, clust, col) for col in df.columns], names=["k", "clust", "rel"]
            )
            return df

        table.fetch_count = fetch_clustered  # type: ignore[assignment]
        # Only write CT for cluster-2-1.
        self._write_ct(_profile_name(self.REG, 2, 1))
        result = abstract_table(table, struct_file=None)
        self.assertIsNotNone(result)
        paired, _ = result
        # Only cluster (2,1) contributes; paired am should be 1/8 only.
        np.testing.assert_array_almost_equal(paired["am"], np.array([1 / 8]))

    def test_min_aucroc_skips_cluster(self):
        # Build a single cluster where mutations are placed only on paired
        # positions (anti-correlates with unpaired → AUC < 0.5 for paired).
        per_position = {
            1: {MUTAT_REL: 1, SUB_C_REL: 1, SUBST_REL: 1},
            2: {MUTAT_REL: 1, SUB_A_REL: 1, SUBST_REL: 1},
            7: {MUTAT_REL: 1, SUB_A_REL: 1, SUBST_REL: 1},
            8: {MUTAT_REL: 1, SUB_C_REL: 1, SUBST_REL: 1},
        }
        counts = _build_counts_df(self.REFSEQ, per_position=per_position)
        table = self._make_table(clusts=[(0, 0)], counts_per_clust={(0, 0): counts})
        self._write_ct(_profile_name(self.REG, 0, 0))
        # AUC for paired ~ 1.0 here (mutations on paired), so set threshold
        # very high to force skipping.
        result = abstract_table(table, struct_file=None, min_aucroc=1.5)
        # All clusters skipped → returns None.
        self.assertIsNone(result)


# ---------------------------------------------------------------------------
# End-to-end composition: multiple tables and multiple SEISMICgraph files
# ---------------------------------------------------------------------------


class TestEndToEndComposition(AbstractTestBase):
    REF = "ref"
    REFSEQ = DNA("ACGTACGT")
    DB_STRING = "((....))"  # paired: 1,2,7,8; unpaired: 3,4,5,6

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _make_table(
        self, sample: str, reg: str, mutated_paired: int, mutated_unpaired: int
    ):
        per_position = {
            1: {
                MUTAT_REL: mutated_paired,
                SUB_C_REL: mutated_paired,
                SUBST_REL: mutated_paired,
            },
            5: {
                MUTAT_REL: mutated_unpaired,
                SUB_G_REL: mutated_unpaired,
                SUBST_REL: mutated_unpaired,
            },
        }
        counts = _build_counts_df(self.REFSEQ, per_position=per_position)
        table = FakeFilterPositionTable(
            top=self.tmpdir,
            sample=sample,
            ref=self.REF,
            refseq=self.REFSEQ,
            reg=reg,
            clusts=[(0, 0)],
            counts_per_clust={(0, 0): counts},
        )
        profile = _profile_name(reg, 0, 0)
        _write_ct(
            _fold_ct_path(self.tmpdir, sample, self.REF, reg, profile),
            self.REF,
            self.REFSEQ,
            reg,
            self.DB_STRING,
            profile,
        )
        return table

    def _write_json(
        self, name: str, reg: str, mutated_paired: int, mutated_unpaired: int
    ):
        # Place mutations at positions 1 (paired A) and 5 (unpaired A).
        per_position = {
            1: {SUB_C_COUNT: mutated_paired, SUBST_COUNT: mutated_paired},
            5: {SUB_G_COUNT: mutated_unpaired, SUBST_COUNT: mutated_unpaired},
        }
        profile = _profile_data(
            len(self.REFSEQ), per_position=per_position, structure=self.DB_STRING
        )
        return _write_seismicgraph_json(
            self.tmpdir / f"{name}.json",
            sample="s",
            entries=[
                {
                    "ref": self.REF,
                    "refseq": self.REFSEQ,
                    "regions": [
                        {
                            "name": reg,
                            "end5": 1,
                            "end3": len(self.REFSEQ),
                            "positions": list(range(1, len(self.REFSEQ) + 1)),
                            "profiles": {f"{reg}__average": profile},
                        }
                    ],
                }
            ],
        )

    def _final_stats(self, table_results, json_results):
        paired, unpaired = _accumulate_ratios(
            chain(filter(None, table_results), json_results)
        )
        paired_means, paired_fvar = _calc_ratio_stats(paired)
        unpaired_means, unpaired_fvar = _calc_ratio_stats(unpaired)
        return (
            paired,
            unpaired,
            paired_means,
            paired_fvar,
            unpaired_means,
            unpaired_fvar,
        )

    def test_one_table_one_json(self):
        table = self._make_table("s1", "r1", 1, 2)
        json_path = self._write_json("export1", "r1", 3, 4)
        results = self._final_stats(
            [abstract_table(table, struct_file=None)],
            [abstract_seismicgraph_file(json_path)],
        )
        paired, unpaired, p_means, p_fvar, u_means, u_fvar = results
        # paired am array: table contributes 1/8; JSON contributes 3/8.
        np.testing.assert_array_almost_equal(
            np.sort(paired["am"]), np.sort(np.array([1 / 8, 3 / 8]))
        )
        # unpaired am array: table 2/8; JSON 4/8.
        np.testing.assert_array_almost_equal(
            np.sort(unpaired["am"]), np.sort(np.array([2 / 8, 4 / 8]))
        )
        # mean and fvar match analytical.
        expected_mean = float(np.nanmean(paired["am"]))
        expected_var = float(np.nanvar(paired["am"]))
        self.assertAlmostEqual(p_means["am"], expected_mean, places=10)
        self.assertAlmostEqual(
            p_fvar["a"],
            expected_var / (expected_mean * (1.0 - expected_mean)),
            places=10,
        )
        self.assertEqual(set(p_fvar), {"a", "c", "g", "t"})

    def test_multiple_tables_and_jsons(self):
        tables = [
            self._make_table("s1", "r1", 1, 2),
            self._make_table("s2", "r2", 5, 6),
        ]
        jsons = [self._write_json("e1", "r1", 3, 4), self._write_json("e2", "r2", 7, 8)]
        table_results = [abstract_table(t, struct_file=None) for t in tables]
        json_results = [abstract_seismicgraph_file(j) for j in jsons]
        paired, unpaired, p_means, p_fvar, u_means, u_fvar = self._final_stats(
            table_results, json_results
        )
        # paired am: 1/8, 3/8, 5/8, 7/8
        np.testing.assert_array_almost_equal(
            np.sort(paired["am"]), np.sort(np.array([1 / 8, 3 / 8, 5 / 8, 7 / 8]))
        )
        # unpaired am: 2/8, 4/8, 6/8, 8/8
        np.testing.assert_array_almost_equal(
            np.sort(unpaired["am"]), np.sort(np.array([2 / 8, 4 / 8, 6 / 8, 8 / 8]))
        )
        # Verify the per-base variance derived from the combined sources.
        combined = paired["am"]
        expected_mean = float(np.nanmean(combined))
        expected_var = float(np.nanvar(combined))
        self.assertAlmostEqual(p_means["am"], expected_mean, places=10)
        self.assertAlmostEqual(
            p_fvar["a"],
            expected_var / (expected_mean * (1.0 - expected_mean)),
            places=10,
        )

    def test_only_one_source_type(self):
        # Single table, no JSON.
        table = self._make_table("s1", "r1", 1, 2)
        paired, unpaired, p_means, p_fvar, _, _ = self._final_stats(
            [abstract_table(table, struct_file=None)], []
        )
        np.testing.assert_array_almost_equal(paired["am"], np.array([1 / 8]))
        # Single value → fvar falls back to margin.
        self.assertEqual(p_fvar["a"], 1.0e-6)
        # Single JSON, no table.
        json_path = self._write_json("e1", "r1", 3, 4)
        paired, _, p_means, p_fvar, _, _ = self._final_stats(
            [], [abstract_seismicgraph_file(json_path)]
        )
        np.testing.assert_array_almost_equal(paired["am"], np.array([3 / 8]))
        self.assertEqual(p_fvar["a"], 1.0e-6)


# ---------------------------------------------------------------------------
# _format_param_tokens — ensures sim abstract output round-trips to sim total
# ---------------------------------------------------------------------------


class TestFormatParamTokens(AbstractTestBase):
    """Guard against regressions where sim abstract emits parameters
    that ``make_pmut_means`` does not accept (e.g. ``pnd``)."""

    def _full_means(self) -> dict[str, float]:
        # Same key set _calc_ratio_stats returns, with arbitrary values.
        means = {"loq": 0.01, "nm": 0.001, "nd": 0.5}
        for base in "acgt":
            means[f"{base}m"] = 0.002
            for other in "acgt":
                if other != base:
                    means[f"{base}{other}"] = 0.30
        return means

    def _full_fvar(self) -> dict[str, float]:
        return {"a": 0.001, "c": 0.001, "g": 0.001, "t": 0.001}

    def test_nd_token_not_emitted(self):
        tokens = _format_param_tokens(
            self._full_means(), self._full_fvar(), self._full_means(), self._full_fvar()
        )
        joined = " ".join(tokens)
        self.assertNotIn(
            " nd ", joined, msg=f"'nd' leaked into params output: {joined}"
        )

    def test_emitted_tokens_round_trip_to_make_pmut_means(self):
        # Build pmut-paired/unpaired tuple lists from the (filtered) means
        # exactly as click would after parsing the printed tokens, then
        # confirm make_pmut_means_paired / make_pmut_means_unpaired accept
        # them without raising.
        paired_means = self._full_means()
        unpaired_means = self._full_means()
        tokens = _format_param_tokens(
            paired_means, self._full_fvar(), unpaired_means, self._full_fvar()
        )
        flag_to_bucket = {
            opt_pmut_paired.opts[-1]: [],
            opt_pmut_unpaired.opts[-1]: [],
            opt_vmut_paired.opts[-1]: [],
            opt_vmut_unpaired.opts[-1]: [],
        }
        for token in tokens:
            flag, key, value = token.split(" ")
            if flag not in flag_to_bucket:
                self.fail(f"Unexpected flag: {flag!r}")
            flag_to_bucket[flag].append((key, float(value)))
        # These should not raise.
        make_pmut_means_paired(
            "SHAPE", **_make_pmut_means_kwargs(flag_to_bucket[opt_pmut_paired.opts[-1]])
        )
        make_pmut_means_unpaired(
            "SHAPE",
            **_make_pmut_means_kwargs(flag_to_bucket[opt_pmut_unpaired.opts[-1]]),
        )
        make_vmut_paired(
            "SHAPE", **_make_vmut_kwargs(flag_to_bucket[opt_vmut_paired.opts[-1]])
        )
        make_vmut_unpaired(
            "SHAPE", **_make_vmut_kwargs(flag_to_bucket[opt_vmut_unpaired.opts[-1]])
        )

    def test_all_other_keys_still_emitted(self):
        # Every means key except "nd" must appear, on both paired and
        # unpaired sides.
        means = self._full_means()
        tokens = _format_param_tokens(
            means, self._full_fvar(), means, self._full_fvar()
        )
        joined = " ".join(tokens)
        pmut_p = opt_pmut_paired.opts[-1]
        pmut_u = opt_pmut_unpaired.opts[-1]
        vmut_p = opt_vmut_paired.opts[-1]
        vmut_u = opt_vmut_unpaired.opts[-1]
        for key in means:
            if key == "nd":
                continue
            self.assertIn(f"{pmut_p} {key} ", joined)
            self.assertIn(f"{pmut_u} {key} ", joined)
        for base in self._full_fvar():
            self.assertIn(f"{vmut_p} {base} ", joined)
            self.assertIn(f"{vmut_u} {base} ", joined)


if __name__ == "__main__":
    ut.main(verbosity=2)
