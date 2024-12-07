from collections import defaultdict
from itertools import chain, product
from pathlib import Path
from typing import Iterable

from click import command

from .write import mask_region
from ..core.arg import (CMD_MASK,
                        arg_input_path,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_mask_coords,
                        opt_mask_primers,
                        opt_primer_gap,
                        opt_mask_regions_file,
                        opt_max_mask_iter,
                        opt_mask_del,
                        opt_mask_ins,
                        opt_mask_mut,
                        opt_mask_polya,
                        opt_mask_gu,
                        opt_mask_pos,
                        opt_mask_pos_file,
                        opt_mask_read,
                        opt_mask_read_file,
                        opt_mask_discontig,
                        opt_min_ncov_read,
                        opt_min_finfo_read,
                        opt_max_fmut_read,
                        opt_min_mut_gap,
                        opt_min_ninfo_pos,
                        opt_max_fmut_pos,
                        opt_quick_unbias,
                        opt_quick_unbias_thresh,
                        opt_mask_pos_table,
                        opt_mask_read_table,
                        opt_brotli_level,
                        opt_max_procs,
                        opt_force,
                        optional_path,
                        extra_defaults)
from ..core.data import load_datasets
from ..core.logs import logger
from ..core.run import run_func
from ..core.seq import DNA, RefRegions
from ..core.task import dispatch
from ..relate.data import load_relate_dataset


def load_regions(input_path: Iterable[str | Path],
                 coords: Iterable[tuple[str, int, int]],
                 primers: Iterable[tuple[str, DNA, DNA]],
                 primer_gap: int,
                 regions_file: Path | None = None):
    """ Open regions of relate reports. """
    # Load all datasets, grouped by their reference names.
    datasets = defaultdict(list)
    for dataset in load_datasets(input_path, load_relate_dataset):
        try:
            datasets[dataset.ref].append(dataset)
        except Exception as error:
            logger.error(error)
    # Determine the regions for each reference in the datasets.
    regions = RefRegions({(loader.ref, loader.refseq)
                          for loader in chain(*datasets.values())},
                         regs_file=regions_file,
                         coords=coords,
                         primers=primers,
                         primer_gap=primer_gap,
                         exclude_primers=True)
    return datasets, regions


@run_func(CMD_MASK, with_tmp=True, extra_defaults=extra_defaults)
def run(input_path: tuple[str, ...], *,
        tmp_dir: Path,
        # Regions
        mask_coords: tuple[tuple[str, int, int], ...],
        mask_primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        mask_regions_file: str | None,
        # Mutation counting
        mask_del: bool,
        mask_ins: bool,
        mask_mut: tuple[str, ...],
        # Filtering
        mask_polya: int,
        mask_gu: bool,
        mask_pos: tuple[tuple[str, int], ...],
        mask_pos_file: str | None,
        mask_read: tuple[str, ...],
        mask_read_file: str | None,
        mask_discontig: bool,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        min_ncov_read: int,
        min_finfo_read: float,
        max_fmut_read: float,
        min_mut_gap: int,
        # Observer bias correction
        quick_unbias: bool,
        quick_unbias_thresh: float,
        # Iteration
        max_mask_iter: int,
        # Table options
        mask_pos_table: bool,
        mask_read_table: bool,
        # Compression
        brotli_level: int,
        # Parallelization
        max_procs: int,
        # Effort
        force: bool) -> list[Path]:
    """ Define mutations and regions to filter reads and positions. """
    # Load all Relate datasets and get the regions for each.
    datasets, regions = load_regions(
        input_path,
        coords=mask_coords,
        primers=mask_primers,
        primer_gap=primer_gap,
        regions_file=optional_path(mask_regions_file)
    )
    # List the datasets and their regions.
    args = [(dataset, region)
            for ref, ref_datasets in datasets.items()
            for dataset, region in product(ref_datasets, regions.list(ref))]
    # Define the keyword arguments.
    kwargs = dict(tmp_dir=tmp_dir,
                  mask_del=mask_del,
                  mask_ins=mask_ins,
                  mask_mut=mask_mut,
                  mask_polya=mask_polya,
                  mask_gu=mask_gu,
                  mask_pos=list(mask_pos),
                  mask_pos_file=optional_path(mask_pos_file),
                  mask_read=list(mask_read),
                  mask_read_file=optional_path(mask_read_file),
                  mask_discontig=mask_discontig,
                  min_ncov_read=min_ncov_read,
                  min_finfo_read=min_finfo_read,
                  max_fmut_read=max_fmut_read,
                  min_mut_gap=min_mut_gap,
                  min_ninfo_pos=min_ninfo_pos,
                  max_fmut_pos=max_fmut_pos,
                  quick_unbias=quick_unbias,
                  quick_unbias_thresh=quick_unbias_thresh,
                  max_mask_iter=max_mask_iter,
                  mask_pos_table=mask_pos_table,
                  mask_read_table=mask_read_table,
                  brotli_level=brotli_level,
                  force=force)
    # Call the mutations and filter the relation vectors.
    return dispatch(mask_region,
                    max_procs=max_procs,
                    args=args,
                    kwargs=kwargs)


params = [
    # Input/output paths
    arg_input_path,
    opt_tmp_pfx,
    opt_keep_tmp,
    # Regions
    opt_mask_coords,
    opt_mask_primers,
    opt_primer_gap,
    opt_mask_regions_file,
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
    opt_mask_read,
    opt_mask_read_file,
    opt_mask_discontig,
    opt_min_ncov_read,
    opt_min_finfo_read,
    opt_max_fmut_read,
    opt_min_mut_gap,
    # Observer bias correction
    opt_quick_unbias,
    opt_quick_unbias_thresh,
    # Iteration
    opt_max_mask_iter,
    # Table options
    opt_mask_pos_table,
    opt_mask_read_table,
    # Compression
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    # Effort
    opt_force,
]


@command(CMD_MASK, params=params)
def cli(*args, **kwargs):
    """ Define mutations and regions to filter reads and positions. """
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
