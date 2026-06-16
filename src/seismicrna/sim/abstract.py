import os
import sys
from collections import defaultdict
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
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_regions_file,
    opt_fold_table_region,
    opt_min_aucroc,
    opt_struct_file,
    opt_verify_times,
    opt_num_cpus,
    opt_pmut_paired,
    opt_pmut_unpaired,
    opt_vmut_paired,
    opt_vmut_unpaired,
    optional_path,
)
from ..core.logs import logger
from ..core.rna import RNAState, from_ct
from ..core.run import run_func
from ..core.seq import DNA, RefRegions, Region, BASE_NAME, BASEN
from ..core.table import (
    COVER_REL,
    INFOR_REL,
    MUTAT_REL,
    SUB_A_REL,
    SUB_C_REL,
    SUB_G_REL,
    SUB_T_REL,
)
from ..core.task import dispatch
from ..core.validate import require_between
from ..filter.table import FilterPositionTableLoader, PartialPositionTable
from ..cluster.data import ClusterPositionTableLoader

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
    return {"loq": (INFOR_REL, COVER_REL), "nm": (MUTAT_REL, INFOR_REL)}


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
    with logger.debug.single_context("{}", message):
        is_unpaired = ~is_paired

        def calc_ratio(num: str, den: str, rows):
            return (counts.loc[rows, num] / counts.loc[rows, den]).dropna().values

        # Get fraction of low-quality bases.
        key = "loq"
        numer, denom = get_other_parameters()[key]
        paired_ratios = {key: 1.0 - calc_ratio(numer, denom, is_paired)}
        logger.trace("paired[{}] = {}", key, paired_ratios[key])
        unpaired_ratios = {key: 1.0 - calc_ratio(numer, denom, is_unpaired)}
        logger.trace("unpaired[{}] = {}", key, unpaired_ratios[key])
        # Get fraction of mutated Ns.
        is_n = counts.index.get_level_values(BASE_NAME) == BASEN
        is_n_paired = is_n & is_paired
        is_n_unpaired = is_n & is_unpaired
        for key, (numer, denom) in get_other_parameters().items():
            if key.startswith(BASEN.lower()):
                paired_ratios[key] = calc_ratio(numer, denom, is_n_paired)
                logger.trace("paired[{}] = {}", key, paired_ratios[key])
                unpaired_ratios[key] = calc_ratio(numer, denom, is_n_unpaired)
                logger.trace("unpaired[{}] = {}", key, unpaired_ratios[key])
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
                    logger.trace("paired[{}] = {}", key, paired_ratios[key])
                    unpaired_ratios[key] = calc_ratio(numer, denom, is_base_unpaired)
                    logger.trace("unpaired[{}] = {}", key, unpaired_ratios[key])
    return paired_ratios, unpaired_ratios


def _accumulate_ratios(paired_unpaired_ratios: Iterable[tuple[dict, dict]]):
    logger.debug("Began accumulating groups of ratios")
    all_paired_ratios = new_parameter_dict()
    all_unpaired_ratios = new_parameter_dict()
    count = 0
    for paired_ratios, unpaired_ratios in paired_unpaired_ratios:
        assert paired_ratios.keys() == all_paired_ratios.keys()
        for key, ratios in paired_ratios.items():
            all_paired_ratios[key] = np.concatenate([all_paired_ratios[key], ratios])
            logger.trace("accum_paired[{}] = {}", key, all_paired_ratios[key])
        assert unpaired_ratios.keys() == all_unpaired_ratios.keys()
        for key, ratios in unpaired_ratios.items():
            all_unpaired_ratios[key] = np.concatenate(
                [all_unpaired_ratios[key], ratios]
            )
            logger.trace("accum_unpaired[{}] = {}", key, all_unpaired_ratios[key])
        count += 1
    logger.debug("Ended accumulating {} group(s) of ratios", count)
    return all_paired_ratios, all_unpaired_ratios


def abstract_table(
    table: PartialPositionTable,
    regions: Iterable[Region] | None,
    ct_files: dict[str, Path],
    *,
    fold_table_region: bool,
    branch: str,
    min_aucroc: float,
):
    # Compute ratios for each cluster, pairing it with its CT file.
    all_ratios = []
    for (hk, hc), profile in table.iter_profiles(
        fold_table_region=fold_table_region, regions=regions
    ):
        # First, get the CT file for the profile, if any.
        ct_file = ct_files.get(profile.profile)
        if ct_file is None:
            # If no CT file was given for this profile via struct_file, then
            # determine the path of the CT file from the fold step.
            ct_file = profile.get_ct_file(table.top, branch)
            if not ct_file.is_file():
                # If no CT file exists from the fold step either, then stop.
                logger.warning(
                    "CT file for {!r} profile {!r} does not exist: {}",
                    table.ref,
                    profile.profile,
                    ct_file,
                )
                continue
        structures = iter(from_ct(ct_file))
        # Check if there is at least one structure in the CT file.
        try:
            struct = next(structures)
        except StopIteration:
            logger.warning("{} contains 0 structures: skipping", ct_file)
            continue
        # Check if there is more than one structure in the CT file.
        try:
            next(structures)
        except StopIteration:
            pass
        else:
            logger.warning("{} contains > 1 structure; using the first", ct_file)
        counts = table.fetch_count(k=hk, clust=hc)
        counts = counts.set_axis(counts.columns.get_level_values(-1), axis=1)
        counts = counts.reindex(struct.is_paired.index)
        if min_aucroc > 0.0:
            state = RNAState.from_struct_profile(struct, profile)
            auc = state.calc_auc()
            if not (auc >= min_aucroc):
                cluster = f" cluster {hk}-{hc}" if hk > 0 else ""
                logger.warning(
                    "Skipping {}{}: AUC-ROC = {} (< {})",
                    table,
                    cluster,
                    auc,
                    min_aucroc,
                )
                continue
        all_ratios.append(_calc_ratios(counts, struct.is_paired.values))
    return _accumulate_ratios(iter(all_ratios))


def _format_param_tokens(
    paired_means: dict, paired_fvar: dict, unpaired_means: dict, unpaired_fvar: dict
) -> list[str]:
    """Format the parameter tokens emitted by ``sim abstract``.

    Skips keys that are not parameters of ``make_pmut_means`` so the
    output can be fed straight back into ``sim total`` / ``sim muts``.
    """
    tokens: list[str] = []
    for key, mean in paired_means.items():
        tokens.append(f"{opt_pmut_paired.opts[-1]} {key} {mean}")
    for base, var in paired_fvar.items():
        tokens.append(f"{opt_vmut_paired.opts[-1]} {base} {var}")
    for key, mean in unpaired_means.items():
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
    fold_table_region: bool,
    fold_coords: Iterable[tuple[str, int, int]],
    fold_primers: Iterable[tuple[str, DNA, DNA]],
    fold_regions_file: str | None,
    struct_file: Iterable[str | Path],
    branch: str,
    min_aucroc: float,
    print_params: bool = True,
    verify_times: bool,
    num_cpus: int,
):
    """Abstract simulation parameters from tables and structures."""
    # Ensure input_path is a list to prevent exhausting an iterable.
    input_path = list(input_path)
    tables = list(
        chain(
            FilterPositionTableLoader.load_tables(
                input_path, verify_times=verify_times
            ),
            ClusterPositionTableLoader.load_tables(
                input_path, verify_times=verify_times
            ),
        )
    )
    # Compute the fold regions from coordinates/primers/file.
    ref_regions = RefRegions(
        {(table.ref, table.refseq) for table in tables},
        regs_file=optional_path(fold_regions_file),
        ends=fold_coords,
        primers=fold_primers,
        default_full=(not fold_table_region),
    )
    # Ensure struct_file is a list to prevent exhausting an iterable.
    struct_file = list(struct_file)
    # Build a mapping from reference and profile name to CT file when
    # explicit structure files are given.
    ref_profile_to_ct: dict[str, dict[str, Path]] = defaultdict(dict)
    for ct in path.find_files_chain(struct_file, path.CT_FILE_LAST_SEGS):
        fields = path.parse(ct, path.CT_FILE_LAST_SEGS)
        ref = fields[path.REF]
        profile = fields[path.PROFILE]
        try:
            ct_curr = ref_profile_to_ct[ref][profile]
        except KeyError:
            # This is a new profile for this reference.
            ref_profile_to_ct[ref][profile] = ct
        else:
            # The profile is repeated.
            raise ValueError(
                "More than one CT file in struct_file matched reference "
                f"{repr(ref)} and profile {repr(profile)}: {ct_curr}, {ct}"
            )
    # List the tables and their region(s).
    args = list()
    for table in tables:
        # If there are no explicit regions for this table, then pass
        # None to use the default region.
        regions = ref_regions.list(table.ref) or None
        # Use the CT files for the reference if there are any, otherwise
        # an empty dict.
        ct_files = ref_profile_to_ct[table.ref]
        if not ct_files and struct_file:
            logger.warning(
                "No CT files given with --struct-file matched tables with reference {!r}",
                table.ref,
            )
        args.append((table, regions, ct_files))
    table_ratios = dispatch(
        abstract_table,
        num_cpus=num_cpus,
        pass_num_cpus=False,
        as_list=False,
        ordered=False,
        raise_on_error=False,
        args=args,
        kwargs=dict(
            fold_table_region=fold_table_region, branch=branch, min_aucroc=min_aucroc
        ),
    )
    # Accumulate all ratios and calculate statistics.
    paired_ratios, unpaired_ratios = _accumulate_ratios(table_ratios)
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
    opt_fold_table_region,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_struct_file,
    opt_branch,
    opt_min_aucroc,
    opt_verify_times,
    opt_num_cpus,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Abstract simulation parameters from tables and structures."""
    run(*args, **kwargs)
