import json
import os
from collections import defaultdict
from functools import cache, partial
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Any

import numpy as np
from click import command

from ..core import path
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.parallel import dispatch
from ..mask.data import MaskMerger
from ..mask.report import MaskReport
from ..relate.data import RelateLoader
from ..relate.report import RelateReport, NumReadsRel
from ..table.base import (COVER_REL,
                          INFOR_REL,
                          SUBST_REL,
                          SUB_A_REL,
                          SUB_C_REL,
                          SUB_G_REL,
                          SUB_T_REL,
                          DELET_REL,
                          INSRT_REL)
from ..table.load import (TableLoader,
                          PosTableLoader,
                          ReadTableLoader,
                          load_tables)

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

params = [
    arg_input_path,
    opt_force,
    opt_max_procs,
    opt_parallel,
]

META_SYMBOL = '#'
SAMPLE = "sample"
REF_SEQ = "sequence"
REF_NUM_ALIGN = "num_aligned"
SECT_END5 = "section_start"
SECT_END3 = "section_end"
SECT_POS = "positions"
POS_DATA = {
    COVER_REL: "cov",
    INFOR_REL: "info",
    SUBST_REL: "sub_N",
    SUB_A_REL: "sub_A",
    SUB_C_REL: "sub_C",
    SUB_G_REL: "sub_G",
    SUB_T_REL: "sub_T",
    DELET_REL: "del",
    INSRT_REL: "ins",
}
SUBST_RATE = "sub_rate"
SUBST_HIST = "sub_hist"
COUNT_PRECISION = 1
RATIO_PRECISION = 6


def format_metadata(metadata: dict[str, Any]):
    """ Prefix each key with the metadata symbol. """
    return {f"{META_SYMBOL}{key}": value for key, value in metadata.items()}


def combine_metadata(special_metadata: dict[str, Any],
                     parsed_metadata: dict[Any, dict],
                     item: Any,
                     what: str = "item"):
    try:
        item_metadata = parsed_metadata[item]
    except KeyError:
        logger.warning(f"No metadata were given for {what} {repr(item)}")
        return special_metadata
    # Check for any keys in the parsed metadata that match those in the
    # special metadata.
    for field in set(special_metadata) & set(item_metadata):
        if (s := special_metadata[field]) != (p := item_metadata[field]):
            raise ValueError(f"Metadata field {repr(field)} of {what} "
                             f"{repr(item)} is {repr(s)}, but was also given "
                             f"as {repr(p)} in the metadata file")
    return format_metadata(special_metadata | item_metadata)


def get_sample_metadata(sample: str,
                        samples_metadata: dict[str, dict]):
    sample_metadata = {SAMPLE: sample}
    return combine_metadata(sample_metadata, samples_metadata, sample, "sample")


def get_ref_metadata(top: Path,
                     sample: str,
                     ref: str,
                     refs_metadata: dict[str, dict]):
    dataset = RelateLoader.load(RelateReport.build_path(top=top,
                                                        sample=sample,
                                                        ref=ref))
    ref_metadata = {REF_SEQ: str(dataset.refseq),
                    REF_NUM_ALIGN: dataset.report.get_field(NumReadsRel)}
    return combine_metadata(ref_metadata, refs_metadata, ref, "reference")


def get_sect_metadata(top: Path,
                      sample: str,
                      ref: str,
                      sect: str,
                      sects_metadata: dict[tuple[str, str], dict]):
    dataset = MaskMerger.load(MaskReport.build_path(top=top,
                                                    sample=sample,
                                                    ref=ref,
                                                    sect=sect))
    sect_metadata = {SECT_END5: dataset.end5,
                     SECT_END3: dataset.end3,
                     SECT_POS: dataset.section.unmasked_int.tolist()}
    return combine_metadata(sect_metadata, sects_metadata, sect, "section")


def iter_pos_table_data(table: PosTableLoader):
    for rel, key in POS_DATA.items():
        table.fetch(ratio=False, precision=COUNT_PRECISION, rels=[rel])
        yield key, table.data[col].values.tolist()
    yield SUBST_RATE, (umd[SUBST_REL] / umd[INFOR_REL]).values.tolist()


def iter_read_table_data(table: ReadTableLoader):
    read_counts = np.asarray(table.data[SUBST_REL].values.round(), dtype=int)
    yield SUBST_HIST, np.bincount(read_counts, minlength=1).tolist()


def iter_table_data(table: TableLoader):
    if isinstance(table, PosTableLoader):
        yield from iter_pos_table_data(table)
    elif isinstance(table, ReadTableLoader):
        yield from iter_read_table_data(table)
    else:
        yield from ()


def get_table_data(table: TableLoader):
    return dict(iter_table_data(table))


def get_sample_data(top: Path,
                    sample: str,
                    tables: list[TableLoader], *,
                    samples_metadata: dict[str, dict],
                    refs_metadata: dict[str, dict],
                    sects_metadata: dict[tuple[str, str], dict]):
    # Cache results from the metadata functions to improve speed.
    ref_metadata = cache(partial(get_ref_metadata,
                                 refs_metadata=refs_metadata))
    sect_metadata = cache(partial(get_sect_metadata,
                                  sects_metadata=sects_metadata))
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
            data[ref][sect] = sect_metadata(top, sample, ref, sect)
        data[ref][sect].update(get_table_data(table))
    return data


def export_sample(top_sample: tuple[Path, str], *args, force: bool, **kwargs):
    top, sample = top_sample
    sample_file = path.buildpar(path.WebAppFileSeg,
                                top=top,
                                sample=sample,
                                ext=path.JSON_EXT)
    if force or not sample_file.is_file():
        with open(sample_file, 'w') as f:
            json.dump(get_sample_data(top, sample, *args, **kwargs), f)
    else:
        logger.warning(f"File exists: {sample_file}")
    return sample_file


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        force: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Export a JSON file of each sample for the DREEM Web App. """
    tables = defaultdict(list)
    for table in load_tables(input_path):
        tables[(table.top, table.sample)].append(table)
    return list(dispatch(export_sample,
                         max_procs,
                         parallel,
                         pass_n_procs=False,
                         args=list(tables.items()),
                         kwargs=dict(samples_metadata=dict(),
                                     refs_metadata=dict(),
                                     sects_metadata=dict(),
                                     force=force)))


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Export a JSON file of each sample for the DREEM Web App. """
    return run(*args, **kwargs)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
