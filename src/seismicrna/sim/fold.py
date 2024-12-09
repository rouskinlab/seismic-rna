import os
from pathlib import Path

from click import command

from .ref import get_fasta_path
from ..core import path
from ..core.arg import (arg_fasta,
                        opt_sim_dir,
                        opt_tmp_pfx,
                        opt_profile_name,
                        opt_fold_regions_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_constraint,
                        opt_fold_temp,
                        opt_fold_md,
                        opt_fold_mfe,
                        opt_fold_max,
                        opt_fold_percent,
                        opt_keep_tmp,
                        opt_force,
                        opt_max_procs,
                        optional_path,
                        extra_defaults)
from ..core.extern import (RNASTRUCTURE_FOLD_CMD,
                           require_dependency,
                           args_to_cmd,
                           run_cmd)
from ..core.rna import renumber_ct
from ..core.run import run_func
from ..core.seq import DNA, RefRegions, Region, parse_fasta, write_fasta
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write
from ..fold.rnastructure import make_fold_cmd, retitle_ct, require_data_path

COMMAND = __name__.split(os.path.extsep)[-1]


def get_ct_path(top: Path, region: Region, profile: str):
    """ Get the path of a connectivity table (CT) file. """
    return path.buildpar(path.RefSeg,
                         path.RegSeg,
                         path.ConnectTableSeg,
                         top=top.joinpath(path.SIM_PARAM_DIR),
                         ref=region.ref,
                         reg=region.name,
                         profile=profile,
                         ext=path.CT_EXT)


def fold_region(region: Region, *,
                sim_dir: Path,
                tmp_dir: Path,
                profile_name: str,
                fold_constraint: Path | None,
                fold_temp: float,
                fold_md: int,
                fold_mfe: bool,
                fold_max: int,
                fold_percent: float,
                keep_tmp: bool,
                force: bool,
                n_procs: int):
    ct_sim = get_ct_path(sim_dir, region, profile_name)
    if need_write(ct_sim, force):
        fasta_tmp = get_fasta_path(tmp_dir, region.ref)
        ct_tmp = get_ct_path(tmp_dir, region, profile_name)
        try:
            # Write a temporary FASTA file for this region only.
            write_fasta(fasta_tmp,
                        [(region.ref, region.seq.tr())],
                        force=force)
            # Predict the RNA structure.
            run_cmd(args_to_cmd(make_fold_cmd(fasta_tmp,
                                              ct_tmp,
                                              dms_file=None,
                                              fold_constraint=fold_constraint,
                                              fold_temp=fold_temp,
                                              fold_md=fold_md,
                                              fold_mfe=fold_mfe,
                                              fold_max=fold_max,
                                              fold_percent=fold_percent,
                                              n_procs=n_procs)))
            # Reformat the CT file title lines so that each is unique.
            retitle_ct(ct_tmp, ct_tmp, force=True)
            # Renumber the CT file so that it has the same numbering
            # scheme as the region, rather than always starting at 1,
            # the latter of which is always output by the Fold program.
            renumber_ct(ct_tmp, ct_sim, region.end5, force=force)
        finally:
            if not keep_tmp:
                fasta_tmp.unlink(missing_ok=True)
                if ct_tmp != ct_sim:
                    ct_tmp.unlink(missing_ok=True)
    return ct_sim


@run_func(COMMAND,
          with_tmp=True,
          pass_keep_tmp=True,
          extra_defaults=extra_defaults)
def run(fasta: str, *,
        sim_dir: str,
        profile_name: str,
        fold_coords: tuple[tuple[str, int, int], ...],
        fold_primers: tuple[tuple[str, DNA, DNA], ...],
        fold_regions_file: str | None,
        fold_constraint: str | None,
        fold_temp: float,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        keep_tmp: bool,
        tmp_dir: Path,
        force: bool,
        max_procs: int):
    # Check for the dependencies and the DATAPATH environment variable.
    require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
    require_data_path()
    # List the regions.
    regions = RefRegions(parse_fasta(Path(fasta), DNA),
                         regs_file=(Path(fold_regions_file)
                                    if fold_regions_file
                                    else None),
                         coords=fold_coords,
                         primers=fold_primers)
    return dispatch(fold_region,
                    max_procs=max_procs,
                    pass_n_procs=True,
                    args=as_list_of_tuples(regions.regions),
                    kwargs=dict(sim_dir=Path(sim_dir),
                                tmp_dir=tmp_dir,
                                profile_name=profile_name,
                                fold_constraint=optional_path(fold_constraint),
                                fold_temp=fold_temp,
                                fold_md=fold_md,
                                fold_mfe=fold_mfe,
                                fold_max=fold_max,
                                fold_percent=fold_percent,
                                keep_tmp=keep_tmp,
                                force=force))


params = [arg_fasta,
          opt_sim_dir,
          opt_tmp_pfx,
          opt_profile_name,
          opt_fold_regions_file,
          opt_fold_coords,
          opt_fold_primers,
          opt_fold_constraint,
          opt_fold_temp,
          opt_fold_md,
          opt_fold_mfe,
          opt_fold_max,
          opt_fold_percent,
          opt_keep_tmp,
          opt_force,
          opt_max_procs]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate secondary structure(s) a reference sequence. """
    run(*args, **kwargs)

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
