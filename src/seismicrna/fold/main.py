from datetime import datetime
from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command

from .report import FoldReport
from .rnastructure import fold, ct2dot
from ..core.arg import (CMD_FOLD,
                        docdef,
                        arg_input_path,
                        opt_temp_dir,
                        opt_keep_temp,
                        opt_sections_file,
                        opt_coords,
                        opt_primers,
                        opt_primer_gap,
                        opt_quantile,
                        opt_fold_temp,
                        opt_fold_constraint,
                        opt_fold_md,
                        opt_fold_mfe,
                        opt_fold_max,
                        opt_fold_percent,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.extern import (RNASTRUCTURE_CT2DOT_CMD,
                           RNASTRUCTURE_FOLD_CMD,
                           require_dependency)
from ..core.parallel import as_list_of_tuples, dispatch, lock_temp_dir
from ..core.rna import RNAProfile
from ..core.seq import DNA, RefSections, RefSeqs, Section
from ..pool.data import load_relate_pool_dataset
from ..relate.report import RelateReport
from ..table.base import MaskPosTable, ClustPosTable
from ..table.load import load_pos_tables

logger = getLogger(__name__)

DEFAULT_QUANTILE = 0.95

params = [
    arg_input_path,
    opt_sections_file,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_quantile,
    opt_fold_temp,
    opt_fold_constraint,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_percent,
    opt_temp_dir,
    opt_keep_temp,
    opt_max_procs,
    opt_parallel,
    opt_force,
]


@command(CMD_FOLD, params=params)
def cli(*args, **kwargs):
    """ Predict RNA secondary structures using mutation rates. """
    return run(*args, **kwargs)


@lock_temp_dir
@docdef.auto()
def run(input_path: tuple[str, ...],
        *,
        sections_file: str | None,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        quantile: float,
        fold_temp: float,
        fold_constraint: str | None,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        temp_dir: str,
        keep_temp: bool,
        max_procs: int,
        parallel: bool,
        force: bool):
    """ Predict RNA secondary structures using mutation rates. """
    require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
    require_dependency(RNASTRUCTURE_CT2DOT_CMD, __name__)
    # Reactivities must be normalized before using them to fold.
    if quantile <= 0.:
        logger.warning("Fold needs normalized mutation rates, but got quantile "
                       f"= {quantile}; setting quantile to {DEFAULT_QUANTILE}")
        quantile = DEFAULT_QUANTILE
    # Gather the tables and reference sequences.
    tables = list()
    ref_seqs = RefSeqs()
    for table in load_pos_tables(input_path):
        # Fold can use only positional tables from Mask and Cluster.
        if isinstance(table, (MaskPosTable, ClustPosTable)):
            tables.append(table)
            # Fetch the reference sequence from the Relate step.
            ref_seqs.add(
                table.ref,
                load_relate_pool_dataset(
                    RelateReport.build_path(top=table.top,
                                            sample=table.sample,
                                            ref=table.ref)
                ).refseq
            )
    # Get the sections for every reference sequence.
    ref_sections = RefSections(ref_seqs,
                               sects_file=(Path(sections_file)
                                           if sections_file
                                           else None),
                               coords=coords,
                               primers=primers,
                               primer_gap=primer_gap)
    # Fold the RNA profiles.
    return list(chain(dispatch(fold_profile,
                               max_procs,
                               parallel,
                               args=[(loader, ref_sections.list(loader.ref))
                                     for loader in tables],
                               kwargs=dict(temp_dir=Path(temp_dir),
                                           keep_temp=keep_temp,
                                           quantile=quantile,
                                           fold_temp=fold_temp,
                                           fold_constraint=(
                                               Path(fold_constraint)
                                               if fold_constraint is not None
                                               else None
                                           ),
                                           fold_md=fold_md,
                                           fold_mfe=fold_mfe,
                                           fold_max=fold_max,
                                           fold_percent=fold_percent,
                                           force=force),
                               pass_n_procs=True)))


def fold_profile(table: MaskPosTable | ClustPosTable,
                 sections: list[Section],
                 n_procs: int,
                 quantile: float,
                 **kwargs):
    """ Fold an RNA molecule from one table of reactivities. """
    return dispatch(fold_section,
                    n_procs,
                    parallel=True,
                    args=as_list_of_tuples(table.iter_profiles(
                        sections=sections, quantile=quantile)
                    ),
                    kwargs=dict(out_dir=table.top,
                                quantile=quantile,
                                **kwargs),
                    pass_n_procs=False)


def fold_section(rna: RNAProfile,
                 out_dir: Path,
                 quantile: float,
                 fold_temp: float,
                 fold_constraint: Path | None,
                 fold_md: int,
                 fold_mfe: bool,
                 fold_max: int,
                 fold_percent: float,
                 force: bool,
                 **kwargs):
    """ Fold a section of an RNA from one mutational profile. """
    began = datetime.now()
    rna.to_varna_color_file(out_dir)
    ct_file = fold(rna,
                   out_dir=out_dir,
                   fold_temp=fold_temp,
                   fold_constraint=fold_constraint,
                   fold_md=fold_md,
                   fold_mfe=fold_mfe,
                   fold_max=fold_max,
                   fold_percent=fold_percent,
                   force=force,
                   **kwargs)
    ct2dot(ct_file)
    ended = datetime.now()
    report = FoldReport(sample=rna.sample,
                        ref=rna.ref,
                        sect=rna.sect,
                        profile=rna.profile,
                        quantile=quantile,
                        fold_temp=fold_temp,
                        fold_md=fold_md,
                        fold_mfe=fold_mfe,
                        fold_max=fold_max,
                        fold_percent=fold_percent,
                        began=began,
                        ended=ended)
    return report.save(out_dir, force=force)

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
