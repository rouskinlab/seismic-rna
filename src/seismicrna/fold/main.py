from logging import getLogger
from pathlib import Path

from click import command

from .rnastructure import fold, ct2dot
from ..core import docdef, path
from ..core.cli import (arg_input_path, opt_temp_dir, opt_save_temp,
                        arg_fasta, opt_sections_file,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_quantile,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.cmd import CMD_FOLD
from ..core.depend import require_dependency
from ..core.fasta import parse_fasta
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rna import RnaProfile
from ..core.sect import RefSections, Section
from ..core.seq import DNA
from ..core.shell import RNASTRUCTURE_FOLD_CMD
from ..core.temp import lock_temp_dir
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
    opt_save_temp,
    opt_max_procs,
    opt_parallel,
    opt_rerun,
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
        save_temp: bool,
        max_procs: int,
        parallel: bool,
        rerun: bool):
    """
    Run the structure module.
    """

    require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)

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
                    kwargs=dict(temp_dir=Path(temp_dir), save_temp=save_temp,
                                quantile=quantile, rerun=rerun),
                    pass_n_procs=True)


def fold_rna(loader: MaskPosTableLoader | ClustPosTableLoader,
             sections: list[Section], n_procs: int, quantile: float, **kwargs):
    """ Fold an RNA molecule from one table of reactivities. """
    return dispatch(fold_profile, n_procs, parallel=True,
                    args=[(profile,)
                          for profile in loader.iter_profiles(sections,
                                                              quantile)],
                    kwargs=dict(out_dir=loader.out_dir, **kwargs),
                    pass_n_procs=False)


def fold_profile(rna: RnaProfile, out_dir: Path, **kwargs):
    """ Fold a section of an RNA from one mutational profile. """
    ct_file = fold(rna, out_dir=out_dir, **kwargs)
    dot_file = ct2dot(ct_file)
    varna_color_file = rna.to_varna_color_file(out_dir)
    return ct_file, dot_file, varna_color_file

########################################################################
#                                                                      #
# ©2023, the Rouskin Lab.                                              #
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
