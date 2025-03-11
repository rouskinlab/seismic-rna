import json
from functools import cache, partial
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .meta import combine_metadata
from ..core import path
from ..core.header import format_clust_name
from ..core.logs import logger
from ..core.rna import parse_db_strings
from ..core.write import need_write, write_mode
from ..fold.rnastructure import parse_energy
from ..mask.dataset import MaskMutsDataset
from ..mask.report import MaskReport
from ..relate.dataset import load_relate_dataset
from ..relate.report import RelateReport
from ..core.table import (COVER_REL,
                          INFOR_REL,
                          MUTAT_REL,
                          SUBST_REL,
                          SUB_A_REL,
                          SUB_C_REL,
                          SUB_G_REL,
                          SUB_T_REL,
                          DELET_REL,
                          INSRT_REL,
                          Table,
                          PositionTable,
                          ReadTable,
                          AbundanceTable)

META_SYMBOL = '#'
SAMPLE = "sample"
REF_SEQ = "sequence"
REF_NUM_ALIGN = "num_aligned"
REG_END5 = "section_start"
REG_END3 = "section_end"
REG_POS = "positions"
COVER_COUNT = "cov"
INFOR_COUNT = "info"
SUBST_COUNT = "sub_N"
SUB_A_COUNT = "sub_A"
SUB_C_COUNT = "sub_C"
SUB_G_COUNT = "sub_G"
SUB_T_COUNT = "sub_T"
DELET_COUNT = "del"
INSRT_COUNT = "ins"
POS_DATA = {COVER_COUNT: COVER_REL,
            INFOR_COUNT: INFOR_REL,
            SUBST_COUNT: SUBST_REL,
            SUB_A_COUNT: SUB_A_REL,
            SUB_C_COUNT: SUB_C_REL,
            SUB_G_COUNT: SUB_G_REL,
            SUB_T_COUNT: SUB_T_REL,
            DELET_COUNT: DELET_REL,
            INSRT_COUNT: INSRT_REL}
SUBST_RATE = "sub_rate"
SUBST_HIST = "sub_hist"
CLUST_PROP = "proportion"
STRUCTURE = "structure"
FREE_ENERGY = "deltaG"
PRECISION = 6


def format_metadata(metadata: dict[str, Any]):
    """ Prefix each key with the metadata symbol. """
    return {f"{META_SYMBOL}{key}": value for key, value in metadata.items()}


def get_sample_metadata(sample: str,
                        samples_metadata: dict[str, dict]):
    sample_metadata = {SAMPLE: sample}
    return format_metadata(combine_metadata(sample_metadata,
                                            samples_metadata,
                                            sample,
                                            "sample"))


def get_ref_metadata(top: Path,
                     sample: str,
                     ref: str,
                     refs_metadata: dict[str, dict]):
    dataset = load_relate_dataset(RelateReport.build_path(
        {path.TOP: top,
         path.SAMPLE: sample,
         path.BRANCHES: dict(),
         path.REF: ref}
    ))
    ref_metadata = {REF_SEQ: str(dataset.refseq),
                    REF_NUM_ALIGN: dataset.num_reads}
    return format_metadata(combine_metadata(ref_metadata,
                                            refs_metadata,
                                            ref,
                                            "reference"))


def get_reg_metadata(top: Path,
                     sample: str,
                     ref: str,
                     reg: str,
                     all_pos: bool):
    dataset = MaskMutsDataset(MaskReport.build_path(
        {path.TOP: top,
         path.SAMPLE: sample,
         path.BRANCHES: dict(),
         path.REF: ref,
         path.REG: reg}
    ))
    positions = (dataset.region.range_int if all_pos
                 else dataset.region.unmasked_int)
    reg_metadata = {REG_END5: dataset.region.end5,
                    REG_END3: dataset.region.end3,
                    REG_POS: positions.tolist()}
    return format_metadata(reg_metadata)


def conform_series(series: pd.Series | pd.DataFrame):
    if isinstance(series, pd.DataFrame):
        if series.columns.size != 1:
            raise TypeError("If series is a DataFrame, then it must have "
                            f"exactly 1 column, but got {series.columns}")
        series = series[series.columns[0]]
    if not isinstance(series, pd.Series):
        raise TypeError(f"Expected Series, but got {type(series.__name__)}")
    return series


def get_db_structs(table: PositionTable,
                   k: int | None = None,
                   clust: int | None = None):
    structs = dict()
    energies = dict()
    for profile in table.iter_profiles(k=k, clust=clust):
        db_file = profile.get_db_file(table.top, branch="")
        if db_file.is_file():
            try:
                # Parse all structures in the dot-bracket file.
                seq, profile_structs = parse_db_strings(db_file)
                # Keep only the first (minimum free energy) structure.
                header, struct = list(profile_structs.items())[0]
                # Parse the minimum free energy of folding.
                energy = parse_energy(header)
            except Exception as error:
                logger.error(error)
            else:
                structs[profile.data_name] = struct
                energies[profile.data_name] = energy
        else:
            logger.warning(f"No structure model available for {profile}: "
                           f"{db_file} does not exist")
    return structs, energies


def iter_pos_table_struct(table: PositionTable, k: int, clust: int):
    structs, energies = get_db_structs(table, k, clust)
    keys = list(structs)
    if keys != list(energies):
        raise ValueError(f"Names of structures {keys} and energies "
                         f"{list(energies)} do not match")
    if keys:
        if len(keys) != 1:
            raise ValueError(f"Expected exactly one structure, but got {keys}")
        key = keys[0]
        yield STRUCTURE, structs[key]
        yield FREE_ENERGY, energies[key]


def iter_pos_table_series(table: PositionTable,
                          k: int,
                          clust: int,
                          all_pos: bool):
    exclude_masked = not all_pos
    for key, rel in POS_DATA.items():
        yield key, conform_series(
            table.fetch_count(rel=rel,
                              k=k,
                              clust=clust,
                              exclude_masked=exclude_masked)
        ).to_list()
    yield SUBST_RATE, conform_series(
        table.fetch_ratio(rel=SUBST_REL,
                          k=k,
                          clust=clust,
                          exclude_masked=exclude_masked,
                          precision=PRECISION)
    ).to_list()


def iter_pos_table_data(table: PositionTable, k: int, clust: int, all_pos: bool):
    yield from iter_pos_table_series(table, k, clust, all_pos)
    yield from iter_pos_table_struct(table, k, clust)


def iter_read_table_data(table: ReadTable, k: int, clust: int):
    read_counts = np.asarray(
        conform_series(table.fetch_count(rel=SUBST_REL,
                                         k=k,
                                         clust=clust)).values,
        dtype=int
    )
    yield SUBST_HIST, np.bincount(read_counts, minlength=1).tolist()


def iter_clust_table_data(table: AbundanceTable, k: int, clust: int):
    clust_count = table.data[table.header.select(k=k,
                                                 clust=clust)].squeeze()
    k_count = table.data[table.header.select(k=k)].sum().squeeze()
    proportion = (round(clust_count / k_count, PRECISION)
                  if k_count > 0
                  else np.nan)
    yield CLUST_PROP, proportion


def iter_table_data(table: Table, k: int, clust: int, all_pos: bool):
    if isinstance(table, PositionTable):
        yield from iter_pos_table_data(table, k, clust, all_pos)
    elif isinstance(table, ReadTable):
        yield from iter_read_table_data(table, k, clust)
    elif isinstance(table, AbundanceTable):
        yield from iter_clust_table_data(table, k, clust)
    else:
        raise TypeError(f"Invalid table type: {type(table).__name__}")


def get_table_data(table: Table, all_pos: bool):
    data = dict()
    for k, clust in table.header.clusts:
        name = format_clust_name(k, clust)
        data[name] = dict(iter_table_data(table, k, clust, all_pos))
    return data


def get_sample_data(top: Path,
                    sample: str,
                    tables: list[Table], *,
                    samples_metadata: dict[str, dict],
                    refs_metadata: dict[str, dict],
                    all_pos: bool):
    # Cache results from the metadata functions to improve speed.
    ref_metadata = cache(partial(get_ref_metadata,
                                 refs_metadata=refs_metadata))
    reg_metadata = cache(get_reg_metadata)
    # Add the metadata for the sample.
    data = get_sample_metadata(sample, samples_metadata)
    # Use a while loop with list.pop() until tables is empty rather
    # than a for loop because each table caches a DataFrame of its
    # data once the data are requested (but not before). Keeping a
    # table in the list after processing it would thus needlessly
    # bloat the memory usage.
    while tables:
        table = tables.pop()
        ref = table.ref
        reg = table.reg
        # Add the metadata for the reference and region.
        if ref not in data:
            data[ref] = ref_metadata(top, sample, ref)
        if reg not in data[ref]:
            data[ref][reg] = reg_metadata(top, sample, ref, reg, all_pos)
        for clust, clust_data in get_table_data(table, all_pos).items():
            if clust not in data[ref][reg]:
                data[ref][reg][clust] = dict()
            data[ref][reg][clust].update(clust_data)
    return data


def export_sample(top_sample: tuple[Path, str], *args, force: bool, **kwargs):
    top, sample = top_sample
    sample_file = path.buildpar([path.WebAppFileSeg],
                                {path.TOP: top,
                                 path.SAMPLE: sample,
                                 path.EXT: path.JSON_EXT})
    if need_write(sample_file, force):
        # Calculate the sample data before opening the file so that the
        # file will not be written if get_sample_data() fails.
        sample_data = get_sample_data(top, sample, *args, **kwargs)
        with open(sample_file, write_mode(force)) as f:
            json.dump(sample_data, f)
    return sample_file
