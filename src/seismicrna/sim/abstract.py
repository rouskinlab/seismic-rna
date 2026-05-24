import json
import os
import sys
from functools import cache
from itertools import chain
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from click import command

from ..core import path
from ..core.arg import (
    arg_input_path,
    opt_branch,
    opt_min_aucroc,
    opt_struct_file,
    opt_verify_times,
    opt_num_cpus,
    opt_pmut_paired,
    opt_pmut_unpaired,
    opt_vmut_paired,
    opt_vmut_unpaired,
)
from ..core.header import format_clust_name
from ..core.logs import logger
from ..core.rna import UNPAIRED_MARK, from_ct
from ..core.rna.roc import compute_auc_roc
from ..core.run import run_func
from ..core.seq import DNA, Region, BASE_NAME, BASEN
from ..core.table import (
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
from ..core.task import as_list_of_tuples, dispatch
from ..core.validate import require_equal, require_between, require_isinstance
from ..export.web import (
    META_SYMBOL,
    REF_SEQ,
    REG_END5,
    REG_END3,
    REG_POS,
    STRUCTURE,
    POS_DATA,
)
from ..filter.table import FilterPositionTableLoader

COMMAND = __name__.split(os.path.extsep)[-1]


@cache
def get_acgt_parameters():
    return {
        "m": (MUTAT_REL, INFOR_REL),
        "a": (SUB_A_REL, MUTAT_REL),
        "c": (SUB_C_REL, MUTAT_REL),
        "g": (SUB_G_REL, MUTAT_REL),
        "t": (SUB_T_REL, MUTAT_REL),
    }


@cache
def get_other_parameters():
    return {
        "loq": (INFOR_REL, COVER_REL),
        "nm": (MUTAT_REL, INFOR_REL),
        "nd": (DELET_REL, MUTAT_REL),
    }


def new_parameter_dict():
    parameters = {code: np.array([], dtype=float) for code in get_other_parameters()}
    for base in DNA.four():
        lower_base = base.lower()
        for code in get_acgt_parameters():
            if lower_base != code:
                parameters[f"{lower_base}{code}"] = np.array([], dtype=float)
    return parameters


def _calc_ratios(counts: pd.DataFrame, is_paired: np.ndarray):
    message = "\n".join(
        map(
            str, ["calculating ratios from", "counts:", counts, "is_paired:", is_paired]
        )
    )
    logger.routine(f"Began {message}")
    is_unpaired = ~is_paired

    def calc_ratio(num: str, den: str, rows):
        return (counts.loc[rows, num] / counts.loc[rows, den]).dropna().values

    # Get fraction of low-quality bases.
    key = "loq"
    numer, denom = get_other_parameters()[key]
    paired_ratios = {key: 1.0 - calc_ratio(numer, denom, is_paired)}
    logger.detail(f"paired[{key}] = {paired_ratios[key]}")
    unpaired_ratios = {key: 1.0 - calc_ratio(numer, denom, is_unpaired)}
    logger.detail(f"unpaired[{key}] = {unpaired_ratios[key]}")
    # Get fraction of mutated Ns.
    is_n = counts.index.get_level_values(BASE_NAME) == BASEN
    is_n_paired = is_n & is_paired
    is_n_unpaired = is_n & is_unpaired
    for key, (numer, denom) in get_other_parameters().items():
        if key.startswith(BASEN.lower()):
            paired_ratios[key] = calc_ratio(numer, denom, is_n_paired)
            logger.detail(f"paired[{key}] = {paired_ratios[key]}")
            unpaired_ratios[key] = calc_ratio(numer, denom, is_n_unpaired)
            logger.detail(f"unpaired[{key}] = {unpaired_ratios[key]}")
    # Get parameters of the four DNA bases.
    for base in DNA.four():
        lower_base = base.lower()
        is_base = counts.index.get_level_values(BASE_NAME) == base
        is_base_paired = is_base & is_paired
        is_base_unpaired = is_base & is_unpaired
        for code, (numer, denom) in get_acgt_parameters().items():
            if lower_base != code:
                key = f"{lower_base}{code}"
                paired_ratios[key] = calc_ratio(numer, denom, is_base_paired)
                logger.detail(f"paired[{key}] = {paired_ratios[key]}")
                unpaired_ratios[key] = calc_ratio(numer, denom, is_base_unpaired)
                logger.detail(f"unpaired[{key}] = {unpaired_ratios[key]}")
    logger.routine(f"Ended {message}")
    return paired_ratios, unpaired_ratios


def _accumulate_ratios(paired_unpaired_ratios: Iterable[tuple[dict, dict]]):
    logger.routine("Began accumulating groups of ratios")
    all_paired_ratios = new_parameter_dict()
    all_unpaired_ratios = new_parameter_dict()
    count = 0
    for paired_ratios, unpaired_ratios in paired_unpaired_ratios:
        assert paired_ratios.keys() == all_paired_ratios.keys()
        for key, ratios in paired_ratios.items():
            all_paired_ratios[key] = np.concatenate([all_paired_ratios[key], ratios])
            logger.detail(f"accum_paired[{key}] = {all_paired_ratios[key]}")
        assert unpaired_ratios.keys() == all_unpaired_ratios.keys()
        for key, ratios in unpaired_ratios.items():
            all_unpaired_ratios[key] = np.concatenate(
                [all_unpaired_ratios[key], ratios]
            )
            logger.detail(f"accum_unpaired[{key}] = {all_unpaired_ratios[key]}")
        count += 1
    logger.routine(f"Ended accumulating {count} group(s) of ratios")
    return all_paired_ratios, all_unpaired_ratios


def abstract_table(
    table: FilterPositionTableLoader,
    struct_file: list[str | Path] | None = None,
    branch: str = "",
    min_aucroc: float = 0.0,
):
    fold_branches = path.add_branch(path.FOLD_STEP, branch, table.branches)
    # Build a mapping from profile name to CT file when explicit files are given.
    profile_to_ct: dict[str, Path] = {}
    if struct_file is not None:
        for sf in path.find_files_chain(struct_file, path.CT_FILE_LAST_SEGS):
            fields = path.parse(sf, path.CT_FILE_LAST_SEGS)
            if fields[path.REF] == table.ref and fields[path.REG] == table.reg:
                profile_to_ct[fields[path.PROFILE]] = sf
        if not profile_to_ct:
            logger.warning(
                f"No CT files in struct_file matched reference {table.ref!r}, "
                f"region {table.reg!r}"
            )
    # Compute ratios for each cluster, pairing it with its CT file.
    all_ratios = []
    for hk, hc in table.header.clusts:
        mus_name = path.fill_whitespace(format_clust_name(hk, hc), fill="-")
        profile_name = f"{table.reg}__{mus_name}"
        # First, get the CT file from the struct_file list, if given.
        ct_file = profile_to_ct.get(profile_name)
        if ct_file is None:
            # If no CT file was given for this table via struct_file, then
            # determine the path of the CT file from the fold step.
            ct_file = path.build(
                path.CT_FILE_ALL_SEGS,
                {
                    path.TOP: table.top,
                    path.SAMPLE: table.sample,
                    path.STEP: path.FOLD_STEP,
                    path.BRANCHES: fold_branches,
                    path.REF: table.ref,
                    path.REG: table.reg,
                    path.PROFILE: profile_name,
                    path.EXT: path.CT_EXT,
                },
            )
            if not ct_file.is_file():
                # If no CT file exists from the fold step either, then stop.
                logger.warning(
                    f"No CT file for {table.ref!r}/{table.reg!r} profile "
                    f"{profile_name!r}: {ct_file}"
                )
                continue
        structures = iter(from_ct(ct_file))
        try:
            struct = next(structures)
        except StopIteration:
            logger.warning(f"{ct_file} contains 0 structures: skipping")
            continue
        try:
            next(structures)
        except StopIteration:
            pass
        else:
            logger.warning(f"{ct_file} contains > 1 structure; using the first")
        require_equal("struct.ref", struct.ref, table.ref, "table.ref", classes=str)
        require_equal(
            "struct.region.name",
            struct.region.name,
            table.region.name,
            "table.region.name",
            classes=str,
        )
        require_equal(
            "struct.region.seq",
            struct.region.seq,
            table.region.seq,
            "table.region.seq",
            classes=DNA,
        )
        counts = table.fetch_count(k=hk, clust=hc)
        counts = counts.set_axis(counts.columns.get_level_values(-1), axis=1)
        if min_aucroc > 0.0:
            auc = compute_auc_roc(
                counts[MUTAT_REL] / counts[INFOR_REL], struct.is_paired
            )
            if not (auc >= min_aucroc):
                logger.detail(
                    f"Skipping {table} cluster {hk}-{hc}: AUC-ROC {auc} < {min_aucroc}"
                )
                continue
        all_ratios.append(_calc_ratios(counts, struct.is_paired.values))
    if not all_ratios:
        return None
    return _accumulate_ratios(iter(all_ratios))


def _abstract_seismicgraph_file(seismicgraph_file: Path, min_aucroc: float = 0.0):
    with open(seismicgraph_file) as f:
        sample_data = json.load(f)
    require_isinstance("sample_data", sample_data, dict)
    for ref in sample_data:
        if ref.startswith(META_SYMBOL):
            # Skip metadata.
            continue
        ref_data = sample_data[ref]
        require_isinstance("ref_data", ref_data, dict)
        refseq = DNA(ref_data[f"{META_SYMBOL}{REF_SEQ}"])
        for reg in ref_data:
            if reg.startswith(META_SYMBOL):
                # Skip metadata.
                continue
            reg_data = ref_data[reg]
            require_isinstance("reg_data", reg_data, dict)
            end5 = reg_data[f"{META_SYMBOL}{REG_END5}"]
            end3 = reg_data[f"{META_SYMBOL}{REG_END3}"]
            positions = reg_data[f"{META_SYMBOL}{REG_POS}"]
            region = Region(ref, refseq, end5=end5, end3=end3)
            region.add_mask("mask", positions, complement=True)
            for profile in reg_data:
                if profile.startswith(META_SYMBOL):
                    # Skip metadata.
                    continue
                profile_data = reg_data[profile]
                require_isinstance("profile_data", profile_data, dict)
                counts = pd.DataFrame.from_dict(
                    {
                        rel: dict(zip(region.unmasked, profile_data[key], strict=True))
                        for key, rel in POS_DATA.items()
                    }
                ).reindex(region.range)
                counts[MUTAT_REL] = (
                    counts[SUBST_REL] + counts[DELET_REL] + counts[INSRT_REL]
                )
                dot_bracket = profile_data[STRUCTURE]
                require_equal(
                    "len(dot_bracket)",
                    len(dot_bracket),
                    region.length,
                    "region.length",
                    classes=int,
                )
                is_paired = np.array([x != UNPAIRED_MARK for x in dot_bracket])
                if min_aucroc > 0.0:
                    auc = compute_auc_roc(
                        counts[MUTAT_REL] / counts[INFOR_REL],
                        pd.Series(is_paired, index=counts.index),
                    )
                    if not (auc >= min_aucroc):
                        logger.detail(
                            f"Skipping profile {profile!r}: "
                            f"AUC-ROC {auc} < {min_aucroc}"
                        )
                        continue
                yield _calc_ratios(counts, is_paired)


def abstract_seismicgraph_file(seismicgraph_file: Path, min_aucroc: float = 0.0):
    return _accumulate_ratios(
        _abstract_seismicgraph_file(seismicgraph_file, min_aucroc)
    )


# Keys present in the means dicts that are NOT accepted by
# `make_pmut_means`.  "nd" (N-base deletion-given-mutation) is computed
# for sanity-checking but the model assumes pnd == 1.0, so emitting it
# would crash `sim total` / `sim muts` with an unexpected-kwarg error.
_NON_MODEL_PMUT_KEYS = frozenset({"nd"})


def _format_param_tokens(
    paired_means: dict, paired_fvar: dict, unpaired_means: dict, unpaired_fvar: dict
) -> list[str]:
    """Format the parameter tokens emitted by ``sim abstract``.

    Skips keys that are not parameters of ``make_pmut_means`` so the
    output can be fed straight back into ``sim total`` / ``sim muts``.
    """
    tokens: list[str] = []
    for key, mean in paired_means.items():
        if key not in _NON_MODEL_PMUT_KEYS:
            tokens.append(f"{opt_pmut_paired.opts[-1]} {key} {mean}")
    for base, var in paired_fvar.items():
        tokens.append(f"{opt_vmut_paired.opts[-1]} {base} {var}")
    for key, mean in unpaired_means.items():
        if key not in _NON_MODEL_PMUT_KEYS:
            tokens.append(f"{opt_pmut_unpaired.opts[-1]} {key} {mean}")
    for base, var in unpaired_fvar.items():
        tokens.append(f"{opt_vmut_unpaired.opts[-1]} {base} {var}")
    return tokens


def _calc_ratio_stats(ratios: dict[str, np.ndarray], margin: float = 1.0e-6):
    require_between("margin", margin, 0.0, 1.0, classes=float, inclusive=False)
    max_value = 1.0 - margin
    means = {
        key: np.clip(np.nanmean(values) if values.size > 0 else 0.0, margin, max_value)
        for key, values in ratios.items()
    }
    fvar = {}
    for key, values in ratios.items():
        if key.endswith("m") and not key.startswith(BASEN.lower()):
            base = key[0]
            if values.size >= 2:
                base_mean = np.clip(np.nanmean(values), margin, max_value)
                fvar[base] = float(np.nanvar(values) / (base_mean * (1.0 - base_mean)))
            else:
                fvar[base] = margin
    return means, fvar


@run_func(COMMAND, default=None)
def run(
    input_path: Iterable[str | Path],
    *,
    struct_file: Iterable[str | Path],
    branch: str,
    min_aucroc: float,
    print_params: bool = True,
    verify_times: bool,
    num_cpus: int,
):
    """Abstract simulation parameters from existing datasets."""
    input_path = list(input_path)
    struct_file_arg = list(struct_file) or None
    # Accumulate ratios from table files.
    args = as_list_of_tuples(
        FilterPositionTableLoader.load_tables(input_path, verify_times=verify_times)
    )
    table_ratios = dispatch(
        abstract_table,
        num_cpus=num_cpus,
        pass_num_cpus=False,
        as_list=False,
        ordered=False,
        raise_on_error=False,
        args=args,
        kwargs=dict(struct_file=struct_file_arg, branch=branch, min_aucroc=min_aucroc),
    )
    # Accumulate ratios from SEISMICgraph files.
    seismicgraph_files = path.find_files_chain(input_path, [path.WebAppFileSeg])
    seismicgraph_ratios = dispatch(
        abstract_seismicgraph_file,
        num_cpus=num_cpus,
        pass_num_cpus=False,
        as_list=False,
        ordered=False,
        raise_on_error=True,
        args=as_list_of_tuples(seismicgraph_files),
        kwargs=dict(min_aucroc=min_aucroc),
    )
    # Accumulate all ratios and calculate statistics.
    paired_ratios, unpaired_ratios = _accumulate_ratios(
        chain(filter(None, table_ratios), seismicgraph_ratios)
    )
    paired_means, paired_fvar = _calc_ratio_stats(paired_ratios)
    unpaired_means, unpaired_fvar = _calc_ratio_stats(unpaired_ratios)
    if print_params:
        params_text = " ".join(
            _format_param_tokens(
                paired_means, paired_fvar, unpaired_means, unpaired_fvar
            )
        )
        sys.stdout.write(f"{params_text}\n")
    return (paired_means, paired_fvar), (unpaired_means, unpaired_fvar)


params = [
    arg_input_path,
    opt_struct_file,
    opt_branch,
    opt_min_aucroc,
    opt_verify_times,
    opt_num_cpus,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Abstract simulation parameters from existing datasets."""
    run(*args, **kwargs)
