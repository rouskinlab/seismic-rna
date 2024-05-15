from collections import defaultdict
from itertools import chain, product
from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .write import mask_section
from ..core.arg import (CMD_MASK,
                        docdef,
                        arg_input_path,
                        opt_coords,
                        opt_primers,
                        opt_primer_gap,
                        opt_sections_file,
                        opt_mask_del,
                        opt_mask_ins,
                        opt_mask_mut,
                        opt_mask_polya,
                        opt_mask_gu,
                        opt_mask_pos_file,
                        opt_mask_pos,
                        opt_mask_discontig,
                        opt_min_ncov_read,
                        opt_min_finfo_read,
                        opt_max_fmut_read,
                        opt_min_mut_gap,
                        opt_min_ninfo_pos,
                        opt_max_fmut_pos,
                        opt_quick_unbias,
                        opt_quick_unbias_thresh,
                        opt_brotli_level,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.data import load_datasets
from ..core.parallel import dispatch
from ..core.seq import DNA, RefSections
from ..pool.data import load_relate_dataset

logger = getLogger(__name__)


def load_sections(input_path: Iterable[str | Path],
                  coords: Iterable[tuple[str, int, int]],
                  primers: Iterable[tuple[str, DNA, DNA]],
                  primer_gap: int,
                  sections_file: Path | None = None):
    """ Open sections of relate reports. """
    # Load all datasets, grouped by their reference names.
    datasets = defaultdict(list)
    for dataset in load_datasets(input_path, load_relate_dataset):
        datasets[dataset.ref].append(dataset)
    # Determine the sections for each reference in the datasets.
    sections = RefSections({(loader.ref, loader.refseq)
                            for loader in chain(*datasets.values())},
                           coords=coords,
                           primers=primers,
                           primer_gap=primer_gap,
                           sects_file=sections_file)
    return datasets, sections


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        # Sections
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        sections_file: str,
        # Mutation counting
        mask_del: bool,
        mask_ins: bool,
        mask_mut: tuple[str, ...],
        # Filtering
        mask_polya: int,
        mask_gu: bool,
        mask_pos: tuple[tuple[str, int], ...],
        mask_pos_file: str | None,
        mask_discontig: bool,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        min_ncov_read: int,
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        # Observer bias correction
        quick_unbias: bool,
        quick_unbias_thresh: float,
        # Compression
        brotli_level: int,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        force: bool) -> list[Path]:
    """ Define mutations and sections to filter reads and positions. """
    # Load all Relate datasets and get the sections for each.
    datasets, sections = load_sections(input_path,
                                       coords=coords,
                                       primers=primers,
                                       primer_gap=primer_gap,
                                       sections_file=(Path(sections_file)
                                                      if sections_file
                                                      else None))
    # List the datasets and their sections.
    args = [(dataset, section)
            for ref, ref_datasets in datasets.items()
            for dataset, section in product(ref_datasets, sections.list(ref))]
    # Define the keyword arguments.
    kwargs = dict(mask_del=mask_del,
                  mask_ins=mask_ins,
                  mask_mut=mask_mut,
                  mask_polya=mask_polya,
                  mask_gu=mask_gu,
                  mask_pos=list(mask_pos),
                  mask_pos_file=(Path(mask_pos_file)
                                 if mask_pos_file
                                 else None),
                  mask_discontig=mask_discontig,
                  min_ncov_read=min_ncov_read,
                  min_finfo_read=min_finfo_read,
                  max_fmut_read=max_fmut_read,
                  min_mut_gap=min_mut_gap,
                  min_ninfo_pos=min_ninfo_pos,
                  max_fmut_pos=max_fmut_pos,
                  quick_unbias=quick_unbias,
                  quick_unbias_thresh=quick_unbias_thresh,
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


params = [
    # Input/output paths
    arg_input_path,
    # Sections
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_sections_file,
    # Mutation counting
    opt_mask_del,
    opt_mask_ins,
    opt_mask_mut,
    # Filtering
    opt_mask_polya,
    opt_mask_gu,
    opt_mask_pos,
    opt_mask_pos_file,
    opt_min_ninfo_pos,
    opt_max_fmut_pos,
    opt_mask_discontig,
    opt_min_ncov_read,
    opt_min_finfo_read,
    opt_max_fmut_read,
    opt_min_mut_gap,
    # Observer bias correction
    opt_quick_unbias,
    opt_quick_unbias_thresh,
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
    """ Define mutations and sections to filter reads and positions. """
    return run(*args, **kwargs)

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
