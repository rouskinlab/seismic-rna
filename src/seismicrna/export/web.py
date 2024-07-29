import json
from functools import cache, partial
from logging import getLogger
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .meta import combine_metadata
from ..core import path
from ..core.header import format_clust_name
from ..core.rna import parse_db_strings
from ..core.write import need_write, write_mode
from ..fold.rnastructure import parse_energy
from ..mask.data import MaskMutsDataset
from ..mask.report import MaskReport
from ..pool.data import load_relate_dataset
from ..relate.report import RelateReport
from ..table.base import (COVER_REL,
                          UNAMB_REL,
                          SUBST_REL,
                          SUB_A_REL,
                          SUB_C_REL,
                          SUB_G_REL,
                          SUB_T_REL,
                          DELET_REL,
                          INSRT_REL)
from ..table.base import (Table,
                          PosTable,
                          ReadTable,
                          ClustFreqTable)

logger = getLogger(__name__)

META_SYMBOL = '#'
SAMPLE = "sample"
REF_SEQ = "sequence"
REF_NUM_ALIGN = "num_aligned"
SECT_END5 = "section_start"
SECT_END3 = "section_end"
SECT_POS = "positions"
POS_DATA = {
    "cov": COVER_REL,
    "info": UNAMB_REL,
    "sub_N": SUBST_REL,
    "sub_A": SUB_A_REL,
    "sub_C": SUB_C_REL,
    "sub_G": SUB_G_REL,
    "sub_T": SUB_T_REL,
    "del": DELET_REL,
    "ins": INSRT_REL,
}
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
    dataset = load_relate_dataset(RelateReport.build_path(top=top,
                                                          sample=sample,
                                                          ref=ref))
    ref_metadata = {REF_SEQ: str(dataset.refseq),
                    REF_NUM_ALIGN: dataset.num_reads}
    return format_metadata(combine_metadata(ref_metadata,
                                            refs_metadata,
                                            ref,
                                            "reference"))


def get_sect_metadata(top: Path,
                      sample: str,
                      ref: str,
                      sect: str,
                      all_pos: bool):
    dataset = MaskMutsDataset.load(MaskReport.build_path(top=top,
                                                         sample=sample,
                                                         ref=ref,
                                                         sect=sect))
    positions = (dataset.section.range_int if all_pos
                 else dataset.section.unmasked_int)
    sect_metadata = {SECT_END5: dataset.end5,
                     SECT_END3: dataset.end3,
                     SECT_POS: positions.tolist()}
    return format_metadata(sect_metadata)


def conform_series(series: pd.Series | pd.DataFrame):
    if isinstance(series, pd.DataFrame):
        if series.columns.size != 1:
            raise TypeError("If series is a DataFrame, then it must have "
                            f"exactly 1 column, but got {series.columns}")
        series = series[series.columns[0]]
    if not isinstance(series, pd.Series):
        raise TypeError(f"Expected Series, but got {type(series.__name__)}")
    return series


def get_db_structs(table: PosTable,
                   k: int | None = None,
                   clust: int | None = None):
    structs = dict()
    energies = dict()
    for profile in table.iter_profiles(k=k, clust=clust):
        db_file = profile.get_db_file(table.top)
        if db_file.is_file():
            try:
                # Parse all structures in the dot-bracket file.
                seq, profile_structs = parse_db_strings(db_file)
                # Keep only the first (minimum free energy) structure.
                header, struct = list(profile_structs.items())[0]
                # Parse the minimum free energy of folding.
                energy = parse_energy(header)
            except Exception as error:
                logger.error("Failed to parse minimum free energy structure "
                             f"from dot-bracket file {db_file}: {error}")
            else:
                structs[profile.data_name] = struct
                energies[profile.data_name] = energy
        else:
            logger.warning(f"No structure model available for {profile} "
                           f"(file {db_file} does not exist)")
    return structs, energies


def iter_pos_table_struct(table: PosTable, k: int, clust: int):
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


def iter_pos_table_series(table: PosTable,
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


def iter_pos_table_data(table: PosTable, k: int, clust: int, all_pos: bool):
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


def iter_clust_table_data(table: ClustFreqTable, k: int, clust: int):
    clust_count = table.data[table.header.select(k=k,
                                                 clust=clust)].squeeze()
    k_count = table.data[table.header.select(k=k)].sum().squeeze()
    proportion = (round(clust_count / k_count, PRECISION)
                  if k_count > 0
                  else np.nan)
    yield CLUST_PROP, proportion


def iter_table_data(table: Table, k: int, clust: int, all_pos: bool):
    if isinstance(table, PosTable):
        yield from iter_pos_table_data(table, k, clust, all_pos)
    elif isinstance(table, ReadTable):
        yield from iter_read_table_data(table, k, clust)
    elif isinstance(table, ClustFreqTable):
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
    sect_metadata = cache(get_sect_metadata)
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
        sect = table.sect
        # Add the metadata for the reference and section.
        if ref not in data:
            data[ref] = ref_metadata(top, sample, ref)
        if sect not in data[ref]:
            data[ref][sect] = sect_metadata(top, sample, ref, sect, all_pos)
        for clust, clust_data in get_table_data(table, all_pos).items():
            if clust not in data[ref][sect]:
                data[ref][sect][clust] = dict()
            data[ref][sect][clust].update(clust_data)
    return data


def export_sample(top_sample: tuple[Path, str], *args, force: bool, **kwargs):
    top, sample = top_sample
    sample_file = path.buildpar(path.WebAppFileSeg,
                                top=top,
                                sample=sample,
                                ext=path.JSON_EXT)
    if need_write(sample_file, force):
        with open(sample_file, write_mode(force)) as f:
            json.dump(get_sample_data(top, sample, *args, **kwargs), f)
    return sample_file

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
