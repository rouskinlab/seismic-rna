import brotli
import json
import os
import pickle
import shutil
from pathlib import Path
from typing import Iterable

from click import command

from .core import path
from .core.arg import (
    CMD_MIGRATE,
    DEFAULT_MUT_COLLISIONS,
    MUT_COLLISIONS_DROP,
    arg_input_path,
    opt_out_dir,
    opt_num_cpus,
)
from .core.logs import logger
from .core.run import run_func
from .core.task import as_list_of_tuples, dispatch

# Field title renames for idmut (was relate) reports.
_IDMUT_FIELD_RENAMES = {
    "Number of reads processed by relate": "Number of reads processed successfully"
}

# Field title renames for filter (was mask) reports.
_FILTER_FIELD_RENAMES = {
    # Option-field title renames.
    "Mask paired-end reads with discontiguous mates": "Drop paired-end reads with discontiguous mates",  # noqa: E501
    "Mask reads with fewer than this many bases covering the region": "Drop reads with fewer than this many bases covering the region",  # noqa: E501
    "Mask reads with less than this fraction of informative base calls": "Drop reads with less than this fraction of informative base calls",  # noqa: E501
    "Mask reads with more than this fraction of mutated base calls": "Drop reads with more than this fraction of mutated base calls",  # noqa: E501
    "Mask reads with two mutations separated by fewer than this many bases": "Filter out mutations separated by fewer than this many bases",  # noqa: E501
    "Stop masking after this many iterations (0 for no limit)": "Stop the filter step after this many iterations (0 for no limit)",  # noqa: E501
    # Result-field title renames.
    "Positions kept after masking": "Positions kept after filtering",
    "Number of positions kept after masking": "Number of positions kept after filtering",  # noqa: E501
    "Total number of reads before masking": "Total number of reads before filtering",
    "Number of reads masked from a list": "Number of reads dropped from a list",
    "Number of reads with too few bases covering the region": "Number of reads dropped due to too few bases covering the region",  # noqa: E501
    "Number of reads with discontiguous mates": "Number of reads dropped due to discontiguous mates",  # noqa: E501
    "Number of reads with too few informative base calls": "Number of reads dropped due to too few informative base calls",  # noqa: E501
    "Number of reads with too many mutations": "Number of reads dropped due to too many mutations",  # noqa: E501
    "Number of reads with two mutations too close": "Number of reads dropped due to two mutations too close",  # noqa: E501
    "Number of reads kept after masking": "Number of reads kept after filtering",
}


def _load_refseq_str(refseq_file: Path) -> str | None:
    """Load a refseq.brickle file and return the DNA sequence as a plain
    uppercase string, or None if the file is missing or unreadable."""
    if not refseq_file.is_file():
        return None
    try:
        with open(refseq_file, "rb") as f:
            data = f.read()
        state = pickle.loads(brotli.decompress(data))
        if isinstance(state, dict):
            compressed = state.get("_s")
        else:
            compressed = getattr(state, "_s", None)
        if compressed is None:
            return None
        from .core.seq.xna import decompress

        return str(decompress(compressed))
    except Exception as error:
        logger.warning(
            "Could not read reference sequence from {}: {}", refseq_file, error
        )
        return None


def _find_refseq_file(ref: str, sample_dir: Path) -> Path | None:
    """Locate the refseq.brickle for ref, first inside sample_dir/idmut/<ref>/
    and then inside sibling sample directories (needed for pool samples that
    do not keep their own refseq copy)."""
    local = sample_dir / "idmut" / ref / "refseq.brickle"
    if local.is_file():
        return local
    for sibling in sample_dir.parent.iterdir():
        if sibling == sample_dir:
            continue
        candidate = sibling / "idmut" / ref / "refseq.brickle"
        if candidate.is_file():
            return candidate
    return None


def _rename_step_dir(sample_dir: Path, old_step: str, new_step: str):
    """Rename sample_dir/<old_step>/ → sample_dir/<new_step>/,
    recursively renaming files whose names begin with <old_step>- to
    begin with <new_step>-.  Idempotent: skips if new_step/ exists."""
    old_dir = sample_dir / old_step
    new_dir = sample_dir / new_step
    if new_dir.exists():
        return
    if not old_dir.is_dir():
        return
    for entry in list(old_dir.rglob("*")):
        if entry.is_file() and entry.name.startswith(f"{old_step}-"):
            entry.rename(
                entry.with_name(entry.name.replace(f"{old_step}-", f"{new_step}-", 1))
            )
    old_dir.rename(new_dir)


def _update_filter_report(data: dict, report_file: Path, sample_dir: Path):
    """Apply filter-specific JSON transformations to a filter-report dict
    in place.  Handles both the straightforward field-title renames from the
    mask→filter rename and the structural differences present in v0.24.4
    reports (combined G/U masking, missing per-base fields).
    JoinFilterReport files are identified by the presence of "Joined regions"
    and get only the option-field updates (not the per-base result fields)."""
    is_join = "Joined regions" in data
    # Rename field titles.
    for old_title, new_title in _FILTER_FIELD_RENAMES.items():
        if old_title in data:
            data[new_title] = data.pop(old_title)
    # Result fields below are absent from JoinFilterReport.
    if is_join:
        return
    # Option fields added after v0.24.4 (not present in JoinFilterReport).
    data.setdefault("Use the default options for this chemical probe", "DMS")
    data.setdefault("Drop reads covering less than this fraction of the region", 0.0)
    # mut_collisions must be stored as a resolved value ("drop" or "merge"),
    # never "auto" — the filter step resolves "auto" before writing the report.
    probe = data.get("Use the default options for this chemical probe", "DMS")
    resolved_mut_collisions = DEFAULT_MUT_COLLISIONS.get(probe, MUT_COLLISIONS_DROP)
    data.setdefault(
        "If two mutations are closer than --min-mut-gap positions, MERGE "
        "the mutations, DROP the read, or AUTO-select based on the probe.",
        resolved_mut_collisions,
    )
    # v0.24.4: convert combined "Mask G and U bases" flag to separate per-base
    # flags.  setdefault preserves individual flags if already present.
    if "Mask G and U bases" in data:
        gu_val = data.pop("Mask G and U bases")
        data.setdefault("Mask positions with base G", gu_val)
        data.setdefault("Mask positions with base U", gu_val)
        data.setdefault("Mask positions with base A", False)
        data.setdefault("Mask positions with base C", False)
    # v0.24.4: split combined "Positions with G or U bases" into separate
    # per-base position lists using the reference sequence.
    if "Positions with G or U bases" in data:
        combined_pos: list[int] = data.pop("Positions with G or U bases")
        data.pop("Number of positions with G or U bases", None)
        ref = report_file.parent.parent.name
        refseq_file = _find_refseq_file(ref, sample_dir)
        refseq = _load_refseq_str(refseq_file) if refseq_file else None
        if refseq is not None:
            g_pos = [
                p for p in combined_pos if 0 < p <= len(refseq) and refseq[p - 1] == "G"
            ]
            u_pos = [
                p for p in combined_pos if 0 < p <= len(refseq) and refseq[p - 1] == "T"
            ]
        else:
            # Fallback when no refseq is available.
            g_pos, u_pos = combined_pos, []
        data.setdefault("Number of positions masked for having base G", len(g_pos))
        data.setdefault("Number of positions masked for having base U", len(u_pos))
        data.setdefault("Positions masked for having base G", g_pos)
        data.setdefault("Positions masked for having base U", u_pos)
    # Add required fields that were absent in v0.24.4 (introduced later or
    # produced by per-base masking options that did not yet exist).
    data.setdefault("Number of positions masked for having base A", 0)
    data.setdefault("Number of positions masked for having base C", 0)
    data.setdefault("Number of positions masked for having base N", 0)
    data.setdefault("Positions masked for having base A", [])
    data.setdefault("Positions masked for having base C", [])
    data.setdefault("Positions masked for having base N", [])
    # opt_min_fcov_read was added after v0.24.4 (commit f6acbc02).
    data.setdefault(
        "Number of reads dropped due to covering too small a fraction of the region", 0
    )


def _update_fold_report(data: dict):
    """Apply fold-specific JSON transformations to a fold-report dict
    in place.  Handles title renames and temperature unit conversion
    (Kelvin → Celsius) present in v0.24.4 reports."""
    # v0.24.4 used "ratios" and "(0.0 disables)" in the quantile title.
    old_quantile = "Normalize and winsorize ratios to this quantile (0.0 disables)"
    new_quantile = "Normalize and winsorize reactivities to this quantile for folding"
    if old_quantile in data:
        data[new_quantile] = data.pop(old_quantile)
    # v0.24.4 stored temperature in Kelvin; v0.25 uses Celsius.
    old_temp = "Predict structures at this temperature (Kelvin)"
    new_temp = "Predict structures at this temperature (Celsius)"
    if old_temp in data:
        data[new_temp] = round(data.pop(old_temp) - 273.15, 10)
    # v0.24.4 only supported RNAstructure; record that explicitly.
    data.setdefault(
        "Model RNA structures using RNAstructure (Fold/ShapeKnots) or "
        "ViennaRNA (RNAfold/RNAsubopt); auto selects RNAstructure for DMS "
        "and ViennaRNA for other probes",
        "RNAstructure",
    )
    # v0.24.4 did not support pseudoknots.
    data.setdefault(
        "Predict pseudoknotted structures (requires "
        "--fold-backend=RNAstructure; uses ShapeKnots when set, "
        "Fold otherwise)",
        False,
    )
    # Fields added after v0.24.4.
    data.setdefault(
        "Only generate the fold command and input files; do not run folding", False
    )
    data.setdefault(
        "Use this method to incorporate reactivities into folding energies. "
        "auto selects Cordero for DMS and Eddy for other probes; Eddy requires "
        "--fold-backend=ViennaRNA; Cordero requires "
        "--fold-backend=RNAstructure",
        "auto",
    )
    data.setdefault(
        "Slope (kcal/mol) for SHAPE reactivities using Deigan method; "
        "used only with --fold-energy-method=Deigan",
        1.8,
    )
    data.setdefault(
        "Intercept (kcal/mol) for SHAPE reactivities using Deigan method; "
        "used only with --fold-energy-method=Deigan",
        -0.6,
    )
    data.setdefault(
        "Maximum absolute energy difference (kcal/mol) from the MFE for "
        "suboptimal structures output by RNAsubopt (overriden by --fold-mfe)",
        1.0,
    )
    data.setdefault("Allow isolated (non-stacked) base pairs when folding", False)
    data.setdefault("Checksum of the ViennaRNA command file (SHA-512)", "")
    data.setdefault("Checksum of the constraints file (SHA-512)", "")
    data.setdefault("Checksum of the Eddy paired-prior file (SHA-512)", "")
    data.setdefault("Checksum of the Eddy unpaired-prior file (SHA-512)", "")


def _update_report_json(report_file: Path, sample_dir: Path):
    """Load a *-report.json file, apply all v0.24→v0.25 field migrations,
    and write it back in place."""
    with open(report_file) as f:
        data = json.load(f)
    # Update Branches dict: relate→idmut, mask→filter.
    branches = data.get("Branches", {})
    if "relate" in branches:
        branches["idmut"] = branches.pop("relate")
    if "mask" in branches:
        branches["filter"] = branches.pop("mask")
    # Update batch-checksum dict (relevant only for idmut/filter reports).
    checksums = data.get("Checksums of batches (SHA-512)", {})
    if "relate" in checksums:
        checksums["idmut"] = checksums.pop("relate")
    if "mask" in checksums:
        checksums["filter"] = checksums.pop("mask")
    # Step-specific field transformations.
    stem = report_file.stem
    if stem == "fold-report" or stem.endswith("__fold-report"):
        _update_fold_report(data)
    elif stem == "idmut-report":
        for old_title, new_title in _IDMUT_FIELD_RENAMES.items():
            if old_title in data:
                data[new_title] = data.pop(old_title)
        # PoolReport fields added after v0.24.4.
        if "Pooled samples" in data:
            data.setdefault(
                "Pool samples only if their Pearson correlation is at least this large",
                0.0,
            )
            data.setdefault(
                "Pool samples only if their mean arsince distance is at most this", 1.0
            )
    elif stem == "filter-report":
        _update_filter_report(data, report_file, sample_dir)
    elif stem == "align-report":
        # SeedF is required (default=None) and was absent in v0.24.4.
        data.setdefault("Seed for the random number generator", 0)
    elif stem == "cluster-report":
        is_join = "Joined regions" in data
        if not is_join:
            # SeedF is required (default=None) and was absent in v0.24.4.
            data.setdefault("Seed for the random number generator", 0)
            # Fields added after v0.24.4.
            data.setdefault(
                "Maximum number of simulations to compute the jackpotting quotient", 12
            )
            data.setdefault(
                "Skip calculating the jackpotting quotient if reads × "
                "positions exceeds this limit",
                268435456,
            )
    with open(report_file, "w") as f:
        json.dump(data, f, indent=4)


def _update_all_reports(sample_dir: Path):
    """Update every *-report.json under sample_dir with v0.24→v0.25 field
    migrations."""
    for report_file in sample_dir.rglob("*-report.json"):
        _update_report_json(report_file, sample_dir)


def migrate_sample_dir(sample_dir: Path, num_cpus: int):
    _rename_step_dir(sample_dir, "relate", "idmut")
    _rename_step_dir(sample_dir, "mask", "filter")
    _update_all_reports(sample_dir)


_STEP_NAMES = ("relate", "mask", "idmut", "filter", "cluster", "fold", "table")


def _is_sample_dir(d: Path) -> bool:
    return d.is_dir() and any((d / step).is_dir() for step in _STEP_NAMES)


def migrate_out_dir(out_dir: Path, num_cpus: int):
    sample_dirs = [d for d in out_dir.iterdir() if _is_sample_dir(d)]
    if not sample_dirs:
        raise ValueError(
            f"{out_dir} contains no sample directories. "
            "Pass the seismic output directory (whose immediate children are "
            "sample directories), not a parent or project root directory."
        )
    dispatch(
        migrate_sample_dir,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        ordered=False,
        raise_on_error=False,
        as_list=True,
        args=as_list_of_tuples(sample_dirs),
    )


@run_func(CMD_MIGRATE)
def run(input_path: Iterable[str | Path], *, out_dir: str | Path, num_cpus: int):
    """Migrate output directories from v0.24 to v0.25"""
    input_path = list(input_path)
    if len(input_path) != 1:
        raise ValueError(
            "seismic migrate can process 1 directory at a time, "
            f"but got {len(input_path)}"
        )
    input_dir = path.sanitize(input_path[0], strict=True)
    out_dir = path.sanitize(out_dir, strict=False)
    if out_dir.exists():
        if os.path.samefile(input_dir, out_dir):
            message = (
                "For safety, seismic migrate refuses to overwrite "
                "existing directories, so the output directory must "
                "be different from the input directory, but got "
                f"input_dir={input_dir}, out_dir={out_dir}. Please "
                "specify a different output directory that does not "
                "yet exist with -o"
            )
        else:
            message = (
                "For safety, seismic migrate refuses to overwrite "
                "existing directories, so the output directory must "
                f"not exist, but got out_dir={out_dir}, which exists. "
                "Please specify a different output directory that "
                "does not yet exist with -o"
            )
        raise FileExistsError(message)
    try:
        shutil.copytree(input_dir, out_dir)
        migrate_out_dir(out_dir, num_cpus=num_cpus)
    except Exception as error:
        if out_dir.exists():
            shutil.rmtree(out_dir, ignore_errors=True)
        raise error


params = [arg_input_path, opt_out_dir, opt_num_cpus]


@command(CMD_MIGRATE, params=params)
def cli(*args, **kwargs):
    """Migrate output directories from v0.24 to v0.25"""
    run(*args, **kwargs)
