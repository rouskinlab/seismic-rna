import os
from abc import ABC, abstractmethod
from collections import defaultdict
from functools import cache, cached_property, partial
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

from click import command
from plotly import graph_objects as go

from ..relate.report import RelateReport, NumReadsRel
from ..relate.data import RelateLoader
from ..core import path
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_rels,
                        opt_y_ratio,
                        opt_quantile,
                        opt_arrange,
                        opt_csv,
                        opt_html,
                        opt_pdf,
opt_quantile,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.parallel import dispatch
from ..table.base import REL_CODES
from ..mask.data import MaskLoader
from ..mask.report import MaskReport
from ..table.load import (TableLoader,
                          MaskPosTableLoader,
                          MaskReadTableLoader,
                          ClustPosTableLoader,
                          load_tables,
                          get_clusters)

logger = getLogger(__name__)

META_SYMBOL = '#'

COMMAND = __name__.split(os.path.extsep)[-1]

params = [
    arg_input_path,
    opt_force,
    opt_max_procs,
    opt_parallel,
]


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
    return special_metadata | item_metadata


def get_sample_metadata(sample: str,
                        samples_metadata: dict[str, dict]):
    sample_metadata = {"sample": sample}
    return format_metadata(combine_metadata(sample_metadata,
                                            samples_metadata,
                                            sample,
                                            "sample"))


def get_ref_metadata(top: Path,
                     sample: str,
                     ref: str,
                     refs_metadata: dict[str, dict]):
    dataset = RelateLoader.load(RelateReport.build_path(top=top,
                                                        sample=sample,
                                                        ref=ref))
    ref_metadata = {"sequence": dataset.refseq,
                    "num_aligned": dataset.report.get_field(NumReadsRel)}
    return format_metadata(combine_metadata(ref_metadata,
                                            refs_metadata,
                                            ref,
                                            "reference"))


def get_section_metadata(top: Path,
                         sample: str,
                         ref: str,
                         sect: str,
                         sects_metadata: dict[tuple[str, str], dict]):
    dataset = MaskLoader.load(MaskReport.build_path(top=top,
                                                    sample=sample,
                                                    ref=ref,
                                                    sect=sect))
    sect_metadata = {"section_start": dataset.end5,
                     "section_end": dataset.end3,
                     "positions": dataset.section.unmasked_int.tolist()}
    return format_metadata(combine_metadata(sect_metadata,
                                            sects_metadata,
                                            sect,
                                            "section"))


def get_table_data(table: TableLoader):
    if isinstance(table, MaskReadTableLoader):
        pass  # FIXME
    elif isinstance(table, MaskPosTableLoader):
        pass


def export_sample(top_sample: tuple[Path, str],
                  tables: list[TableLoader], *,
                  samples_metadata: dict[str, dict],
                  refs_metadata: dict[str, dict],
                  sects_metadata: dict[tuple[str, str], dict],
                  force: bool):
    top, sample = top_sample
    sample_file = path.buildpar(path.WebAppFileSeg,
                                top=top,
                                sample=sample,
                                ext=path.JSON_EXT)
    if force or not sample_file.is_file():
        # Add the metadata for the sample.
        data = get_metadata(samples_metadata, sample, "sample")
        # Use a while loop with list.pop() until tables is empty rather
        # than a for loop because each table caches a DataFrame of its
        # data once the data are requested (but not before). Keeping a
        # table in the list after processing it would thus needlessly
        # bloat the memory usage.
        while tables:
            table = tables.pop()
            ref_sect = table.ref, table.sect
            ref, sect = ref_sect
            # Add the metadata for the reference and section.
            if ref not in data:
                data[ref] = get_metadata(refs_metadata, ref, "reference")
            if ref_sect not in data[ref]:
                data[ref][sect] = get_metadata(sects_metadata,
                                               ref_sect,
                                               "section")


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        force: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Export a JSON file of each sample for the DREEM Web App. """
    tables = defaultdict(list)
    for table in load_tables(input_path):
        tables[(table.top, table.sample)].append(table)
    return list(chain(*dispatch(export_sample,
                                max_procs,
                                parallel,
                                pass_n_procs=False,
                                args=list(tables.items()),
                                kwargs=dict(force=force))))


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
