from logging import getLogger
from pathlib import Path

from click import command

from .rnastructure import fold, ct2dot
from ..core import path
from ..core.arg import (CMD_FOLD,
                        docdef,
                        arg_input_path,
                        opt_temp_dir,
                        opt_keep_temp,
                        arg_fasta,
                        opt_sections_file,
                        opt_coords,
                        opt_primers,
                        opt_primer_gap,
                        opt_quantile,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.extern import (RNASTRUCTURE_CT2DOT_CMD,
                           RNASTRUCTURE_FOLD_CMD,
                           require_dependency)
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rna import RnaProfile
from ..core.seq import DNA, RefSections, Section, parse_fasta
from ..core.parallel import lock_temp_dir
from ..table.load import load, MaskPosTableLoader, ClustPosTableLoader

logger = getLogger(__name__)

params = [
    arg_fasta,
    arg_input_path,
    opt_sections_file,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_quantile,
    opt_temp_dir,
    opt_keep_temp,
    opt_max_procs,
    opt_parallel,
    opt_force,
]


@command(CMD_FOLD, params=params)
def cli(*args, **kwargs):
    """ Predict the structure(s) of an RNA using mutation rates from the
    individual clusters or the ensemble average ('mask' step). """
    return run(*args, **kwargs)


@lock_temp_dir
@docdef.auto()
def run(fasta: str,
        input_path: tuple[str, ...],
        *,
        sections_file: str | None,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        quantile: float,
        temp_dir: str,
        keep_temp: bool,
        max_procs: int,
        parallel: bool,
        force: bool):
    """
    Predict RNA structures using mutation rates as constraints.
    """

    require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
    require_dependency(RNASTRUCTURE_CT2DOT_CMD, __name__)

    # Reactivities must be normalized before using them to fold.
    if quantile <= 0.:
        logger.warning("Fold requires normalized mutation rates, but got "
                       f"quantile = {quantile}; setting quantile to 1.0")
        quantile = 1.
    # Get the sections for every reference sequence.
    ref_sections = RefSections(parse_fasta(Path(fasta), DNA),
                               sects_file=(Path(sections_file) if sections_file
                                           else None),
                               coords=coords,
                               primers=primers,
                               primer_gap=primer_gap)
    # Initialize the table loaders.
    tab_files = path.find_files_chain(map(Path, input_path), [path.TableSeg])
    loaders = [loader for loader in dispatch(load, max_procs, parallel,
                                             args=as_list_of_tuples(tab_files),
                                             pass_n_procs=False)
               if isinstance(loader, (MaskPosTableLoader, ClustPosTableLoader))]
    # Fold the RNA profiles.
    return dispatch(fold_rna, max_procs, parallel,
                    args=[(loader, ref_sections.list(loader.ref))
                          for loader in loaders],
                    kwargs=dict(temp_dir=Path(temp_dir), keep_temp=keep_temp,
                                quantile=quantile, force=force),
                    pass_n_procs=True)


def fold_rna(loader: MaskPosTableLoader | ClustPosTableLoader,
             sections: list[Section], n_procs: int, quantile: float, **kwargs):
    """ Fold an RNA molecule from one table of reactivities. """
    return dispatch(fold_profile, n_procs, parallel=True,
                    args=[(profile,)
                          for profile in loader.iter_profiles(sections,
                                                              quantile)],
                    kwargs=dict(out_dir=loader.top, **kwargs),
                    pass_n_procs=False)


def fold_profile(rna: RnaProfile, out_dir: Path, **kwargs):
    """ Fold a section of an RNA from one mutational profile. """
    ct_file = fold(rna, out_dir=out_dir, **kwargs)
    dot_file = ct2dot(ct_file)
    varna_color_file = rna.to_varna_color_file(out_dir)
    return ct_file, dot_file, varna_color_file

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
