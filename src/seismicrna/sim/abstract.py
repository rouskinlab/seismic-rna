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
from ..core.arg import (arg_input_path,
                        opt_struct_file,
                        opt_verify_times,
                        opt_num_cpus,
                        opt_pmut_paired,
                        opt_pmut_unpaired,
                        opt_vmut_paired,
                        opt_vmut_unpaired)
from ..core.error import DuplicateValueError, NoDataError
from ..core.logs import logger
from ..core.rna import UNPAIRED_MARK, from_ct
from ..core.run import run_func
from ..core.seq import DNA, Region, BASE_NAME, BASEN
from ..core.table import (COVER_REL,
                          INFOR_REL,
                          MUTAT_REL,
                          SUBST_REL,
                          SUB_A_REL,
                          SUB_C_REL,
                          SUB_G_REL,
                          SUB_T_REL,
                          DELET_REL,
                          INSRT_REL)
from ..core.task import as_list_of_tuples, dispatch
from ..core.validate import require_equal, require_between, require_isinstance
from ..export.web import (META_SYMBOL,
                          REF_SEQ,
                          REG_END5,
                          REG_END3,
                          REG_POS,
                          STRUCTURE,
                          POS_DATA)
from ..mask.table import MaskPositionTableLoader

COMMAND = __name__.split(os.path.extsep)[-1]


@cache
def get_acgt_parameters():
    return {"m": (MUTAT_REL, INFOR_REL),
            "a": (SUB_A_REL, MUTAT_REL),
            "c": (SUB_C_REL, MUTAT_REL),
            "g": (SUB_G_REL, MUTAT_REL),
            "t": (SUB_T_REL, MUTAT_REL)}


@cache
def get_other_parameters():
    return {"loq": (INFOR_REL, COVER_REL),
            "nm": (MUTAT_REL, INFOR_REL),
            "nd": (DELET_REL, MUTAT_REL)}


def new_parameter_dict():
    parameters = {code: np.array([], dtype=float)
                  for code in get_other_parameters()}
    for base in DNA.four():
        lower_base = base.lower()
        for code in get_acgt_parameters():
            if lower_base != code:
                parameters[f"{lower_base}{code}"] = np.array([], dtype=float)
    return parameters


def _calc_ratios(counts: pd.DataFrame, is_paired: np.ndarray):
    message = "\n".join(map(str, ["calculating ratios from",
                                  "counts:",
                                  counts,
                                  "is_paired:",
                                  is_paired]))
    logger.routine(f"Began {message}")
    is_unpaired = ~is_paired

    def calc_ratio(num: str, den: str, rows):
        return (counts.loc[rows, num] / counts.loc[rows, den]).dropna().values

    # Get fraction of low-quality bases.
    key = "loq"
    numer, denom = get_other_parameters()[key]
    paired_ratios = {key: 1. - calc_ratio(numer, denom, is_paired)}
    logger.detail(f"paired[{key}] = {paired_ratios[key]}")
    unpaired_ratios = {key: 1. - calc_ratio(numer, denom, is_unpaired)}
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
                paired_ratios[key] = calc_ratio(numer,
                                                denom,
                                                is_base_paired)
                logger.detail(f"paired[{key}] = {paired_ratios[key]}")
                unpaired_ratios[key] = calc_ratio(numer,
                                                  denom,
                                                  is_base_unpaired)
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
            all_paired_ratios[key] = np.concatenate([all_paired_ratios[key],
                                                     ratios])
            logger.detail(f"accum_paired[{key}] = {all_paired_ratios[key]}")
        assert unpaired_ratios.keys() == all_unpaired_ratios.keys()
        for key, ratios in unpaired_ratios.items():
            all_unpaired_ratios[key] = np.concatenate([all_unpaired_ratios[key],
                                                       ratios])
            logger.detail(f"accum_unpaired[{key}] = {all_unpaired_ratios[key]}")
        count += 1
    logger.routine(f"Ended accumulating {count} group(s) of ratios")
    return all_paired_ratios, all_unpaired_ratios


def abstract_table(table: MaskPositionTableLoader,
                   struct_file: str | Path):
    path_fields = path.parse(struct_file, path.CT_FILE_LAST_SEGS)
    require_equal("table.ref",
                  table.ref,
                  path_fields[path.REF],
                  "path_fields[ref]",
                  classes=str)
    structures = iter(from_ct(struct_file))
    try:
        struct = next(structures)
    except StopIteration:
        raise NoDataError(f"{struct_file} contains 0 structures") from None
    try:
        next(structures)
    except StopIteration:
        pass
    else:
        logger.warning(f"{struct_file} contains > 1 structure; using the first")
    require_equal("struct.ref",
                  struct.ref,
                  table.ref,
                  "table.ref",
                  classes=str)
    require_equal("struct.region.name",
                  struct.region.name,
                  table.region.name,
                  "table.region.name",
                  classes=str)
    require_equal("struct.region.seq",
                  struct.region.seq,
                  table.region.seq,
                  "table.region.seq",
                  classes=DNA)
    return _calc_ratios(table.fetch_count(), struct.is_paired.values)


def _abstract_seismicgraph_file(seismicgraph_file: Path):
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
                    {rel: dict(zip(region.unmasked,
                                   profile_data[key],
                                   strict=True))
                     for key, rel in POS_DATA.items()}
                ).reindex(region.range)
                counts[MUTAT_REL] = (counts[SUBST_REL]
                                     + counts[DELET_REL]
                                     + counts[INSRT_REL])
                dot_bracket = profile_data[STRUCTURE]
                require_equal("len(dot_bracket)",
                              len(dot_bracket),
                              region.length,
                              "region.length",
                              classes=int)
                is_paired = np.array([x != UNPAIRED_MARK for x in dot_bracket])
                yield _calc_ratios(counts, is_paired)


def abstract_seismicgraph_file(seismicgraph_file: Path):
    return _accumulate_ratios(_abstract_seismicgraph_file(seismicgraph_file))


def _calc_ratio_stats(ratios: dict[str, np.ndarray], margin: float = 1.e-6):
    require_between("margin", margin, 0., 1., classes=float)
    max_value = 1. - margin
    means = {key: np.clip(np.nanmean(values) if values.size > 0 else 0.,
                          margin,
                          max_value)
             for key, values in ratios.items()}
    if ratios:
        ms = np.concatenate([values for key, values in ratios.items()
                             if key.endswith("m")])
    else:
        ms = np.array([], dtype=float)
    if ms.size >= 2:
        mean = np.clip(np.nanmean(ms), margin, max_value)
        fvar = np.nanvar(ms) / (mean * (1. - mean))
    else:
        fvar = 0.
    return means, fvar


@run_func(COMMAND, default=None)
def run(input_path: Iterable[str | Path], *,
        struct_file: Iterable[str | Path],
        print_params: bool = True,
        verify_times: bool,
        num_cpus: int):
    """ Abstract simulation parameters from existing datasets. """
    input_path = list(input_path)
    # Group structure files by reference.
    ref_struct_files = dict()
    for file in path.find_files_chain(struct_file, path.CT_FILE_LAST_SEGS):
        ref = path.parse(file, path.CT_FILE_LAST_SEGS)[path.REF]
        if ref in ref_struct_files:
            raise DuplicateValueError(
                f"More than one structure file for reference {repr(ref)}"
            )
        ref_struct_files[ref] = file
    # Accumulate ratios from table files.
    args = list()
    for table in MaskPositionTableLoader.load_tables(
            input_path, verify_times=verify_times
    ):
        try:
            args.append((table, ref_struct_files[table.ref]))
        except KeyError:
            logger.error(f"No structure for reference {repr(table.ref)}")
    table_ratios = dispatch(abstract_table,
                            num_cpus=num_cpus,
                            pass_num_cpus=False,
                            as_list=False,
                            ordered=False,
                            raise_on_error=True,
                            args=args)
    # Accumulate ratios from SEISMICgraph files.
    seismicgraph_files = path.find_files_chain(input_path, [path.WebAppFileSeg])
    seismicgraph_ratios = dispatch(abstract_seismicgraph_file,
                                   num_cpus=num_cpus,
                                   pass_num_cpus=False,
                                   as_list=False,
                                   ordered=False,
                                   raise_on_error=True,
                                   args=as_list_of_tuples(seismicgraph_files))
    # Accumulate all ratios and calculate statistics.
    paired_ratios, unpaired_ratios = _accumulate_ratios(
        chain(table_ratios, seismicgraph_ratios)
    )
    paired_means, paired_fvar = _calc_ratio_stats(paired_ratios)
    unpaired_means, unpaired_fvar = _calc_ratio_stats(unpaired_ratios)
    if print_params:
        params_items = [f"{opt_vmut_paired.opts[-1]} {paired_fvar}",
                        f"{opt_vmut_unpaired.opts[-1]} {unpaired_fvar}"]
        params_items.extend(f"{opt_pmut_paired.opts[-1]} {key} {mean}"
                            for key, mean in paired_means.items())
        params_items.extend(f"{opt_pmut_unpaired.opts[-1]} {key} {mean}"
                            for key, mean in unpaired_means.items())
        params_text = f"Abstracted parameters:\n{' '.join(params_items)}\n"
        sys.stdout.write(params_text)
    return (paired_means, paired_fvar), (unpaired_means, unpaired_fvar)


params = [arg_input_path,
          opt_struct_file,
          opt_verify_times,
          opt_num_cpus]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Abstract simulation parameters from existing datasets. """
    run(*args, **kwargs)
