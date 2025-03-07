import os
import sys
from functools import cache
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from ..core import path
from ..core.arg import (arg_input_path,
                        opt_struct_file,
                        opt_verify_times,
                        opt_max_procs,
                        opt_pmut_paired,
                        opt_pmut_unpaired,
                        opt_vmut_paired,
                        opt_vmut_unpaired)
from ..core.error import DuplicateValueError, NoDataError
from ..core.logs import logger
from ..core.rna import from_ct
from ..core.run import run_func
from ..core.seq import DNA, BASE_NAME, BASEN
from ..core.table import (COVER_REL,
                          INFOR_REL,
                          MUTAT_REL,
                          SUB_A_REL,
                          SUB_C_REL,
                          SUB_G_REL,
                          SUB_T_REL,
                          DELET_REL)
from ..core.task import dispatch
from ..core.validate import require_equal, require_between
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
    is_paired = struct.is_paired
    is_unpaired = ~is_paired
    counts = table.fetch_count()

    def calc_ratio(num: str, den: str, rows):
        return (counts.loc[rows, num] / counts.loc[rows, den]).dropna().values

    # Get fraction of low-quality bases.
    key = "loq"
    numer, denom = get_other_parameters()[key]
    ratios_paired = {key: 1. - calc_ratio(numer, denom, is_paired)}
    ratios_unpaired = {key: 1. - calc_ratio(numer, denom, is_unpaired)}
    # Get fraction of mutated Ns.
    is_n = counts.index.get_level_values(BASE_NAME) == BASEN
    is_n_paired = is_n & is_paired
    is_n_unpaired = is_n & is_unpaired
    for key, (numer, denom) in get_other_parameters().items():
        if key.startswith(BASEN.lower()):
            ratios_paired[key] = calc_ratio(numer, denom, is_n_paired)
            ratios_unpaired[key] = calc_ratio(numer, denom, is_n_unpaired)
    # Get parameters of the four DNA bases.
    for base in DNA.four():
        lower_base = base.lower()
        is_base = counts.index.get_level_values(BASE_NAME) == base
        is_base_paired = is_base & is_paired
        is_base_unpaired = is_base & is_unpaired
        for code, (numer, denom) in get_acgt_parameters().items():
            if lower_base != code:
                key = f"{lower_base}{code}"
                ratios_paired[key] = calc_ratio(numer,
                                                denom,
                                                is_base_paired)
                ratios_unpaired[key] = calc_ratio(numer,
                                                  denom,
                                                  is_base_unpaired)
    return ratios_paired, ratios_unpaired


def calc_ratio_stats(ratios: dict[str, np.ndarray], margin: float = 1.e-6):
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
        verify_times: bool,
        max_procs: int):
    """ Abstract simulation parameters from existing datasets. """
    # Group structure files by reference.
    ref_struct_files = dict()
    for file in path.find_files_chain(struct_file, path.CT_FILE_LAST_SEGS):
        ref = path.parse(file, path.CT_FILE_LAST_SEGS)[path.REF]
        if ref in ref_struct_files:
            raise DuplicateValueError(
                f"More than one structure file for reference {repr(ref)}"
            )
        ref_struct_files[ref] = file
    args = list()
    for table in MaskPositionTableLoader.load_tables(
            input_path, verify_times=verify_times
    ):
        try:
            args.append((table, ref_struct_files[table.ref]))
        except KeyError:
            logger.error(f"No structure for reference {repr(table.ref)}")
    paired_ratios = new_parameter_dict()
    unpaired_ratios = new_parameter_dict()
    for p, u in dispatch(abstract_table,
                         max_procs=max_procs,
                         pass_n_procs=False,
                         args=args):
        assert p.keys() == paired_ratios.keys()
        # Iterate over list(paired_ratios) to avoid iterating over an
        # object that is being modified.
        for key in list(paired_ratios):
            paired_ratios[key] = np.concatenate([paired_ratios[key],
                                                 p[key]])
        assert u.keys() == unpaired_ratios.keys()
        # Iterate over list(unpaired_ratios) to avoid iterating over an
        # object that is being modified.
        for key in list(unpaired_ratios):
            unpaired_ratios[key] = np.concatenate([unpaired_ratios[key],
                                                   u[key]])
    paired_means, paired_fvar = calc_ratio_stats(paired_ratios)
    unpaired_means, unpaired_fvar = calc_ratio_stats(unpaired_ratios)
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
          opt_max_procs]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Abstract simulation parameters from existing datasets. """
    run(*args, **kwargs)
