from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .write import mask_section
from ..core import path
from ..core.arg import (CMD_MASK,
                        docdef,
                        arg_input_path,
                        opt_coords,
                        opt_primers,
                        opt_primer_gap,
                        opt_sections_file,
                        opt_count_del,
                        opt_count_ins,
                        opt_discount_mut,
                        opt_exclude_polya,
                        opt_exclude_gu,
                        opt_exclude_pos,
                        opt_min_finfo_read,
                        opt_max_fmut_read,
                        opt_min_mut_gap,
                        opt_min_ninfo_pos,
                        opt_max_fmut_pos,
                        opt_brotli_level,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.data import load_data
from ..core.parallel import dispatch
from ..core.seq import DNA, RefSections
from ..relate.data import RelateLoader

logger = getLogger(__name__)

params = [
    # Input/output paths
    arg_input_path,
    # Sections
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_sections_file,
    # Mutation counting
    opt_count_del,
    opt_count_ins,
    opt_discount_mut,
    # Filtering
    opt_exclude_polya,
    opt_exclude_gu,
    opt_exclude_pos,
    opt_min_finfo_read,
    opt_max_fmut_read,
    opt_min_mut_gap,
    opt_min_ninfo_pos,
    opt_max_fmut_pos,
    # Compression
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Effort
    opt_force,
]


@command(CMD_MASK, params=params)
def cli(*args, **kwargs):
    """ Select a section of the reference, define which relationships
    count as mutations, and filter out unusable positions and reads. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        # Sections
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        sections_file: str,
        # Mutation counting
        count_del: bool,
        count_ins: bool,
        discount_mut: tuple[str, ...],
        # Filtering
        exclude_polya: int,
        exclude_gu: bool,
        exclude_pos: tuple[tuple[str, int], ...],
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        # Compression
        brotli_level: int,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        force: bool) -> list[Path]:
    """ Run the mask command. """
    # Load all relation vector datasets and get the sections for each.
    datasets, sections = load_sections(map(Path, input_path),
                                       coords=coords,
                                       primers=primers,
                                       primer_gap=primer_gap,
                                       sections_file=(Path(sections_file)
                                                      if sections_file
                                                      else None))
    # List the relation datasets and their sections.
    args = [(dataset, section) for dataset in datasets
            for section in sections.list(dataset.ref)]
    # Define the keyword arguments.
    kwargs = dict(count_del=count_del,
                  count_ins=count_ins,
                  discount=discount_mut,
                  exclude_polya=exclude_polya,
                  exclude_gu=exclude_gu,
                  exclude_pos=exclude_pos,
                  min_finfo_read=min_finfo_read,
                  max_fmut_read=max_fmut_read,
                  min_mut_gap=min_mut_gap,
                  min_ninfo_pos=min_ninfo_pos,
                  max_fmut_pos=max_fmut_pos,
                  brotli_level=brotli_level,
                  force=force)
    # Call the mutations and filter the relation vectors.
    reports = dispatch(mask_section,
                       max_procs=max_procs,
                       parallel=parallel,
                       pass_n_procs=False,
                       args=args,
                       kwargs=kwargs)
    return list(map(Path, reports))


def load_sections(report_files: Iterable[Path],
                  coords: Iterable[tuple[str, int, int]],
                  primers: Iterable[tuple[str, DNA, DNA]],
                  primer_gap: int,
                  sections_file: Path | None = None):
    """ Open sections of relate reports. """
    datasets = list(load_data(path.find_files_chain(report_files,
                                                    [path.RelateRepSeg]),
                              RelateLoader))
    sections = RefSections({(loader.ref, loader.refseq) for loader in datasets},
                           coords=coords,
                           primers=primers,
                           primer_gap=primer_gap,
                           sects_file=sections_file)
    return datasets, sections

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
