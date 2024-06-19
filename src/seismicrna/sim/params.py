import os
from logging import getLogger

from click import command

from . import (clusts as clusts_mod,
               ends as ends_mod,
               muts as muts_mod)
from ..core.arg import merge_params
from ..core.run import run_func

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


@run_func(logger.critical)
def run(*,
        ct_file: tuple[str, ...],
        pmut_paired: tuple[tuple[str, float], ...],
        pmut_unpaired: tuple[tuple[str, float], ...],
        vmut_paired: float,
        vmut_unpaired: float,
        end3_fmean: float,
        insert_fmean: float,
        ends_var: float,
        clust_conc: float,
        force: bool,
        max_procs: int,
        parallel: bool):
    """ Simulate parameter files. """
    muts_mod.run(ct_file=ct_file,
                 pmut_paired=pmut_paired,
                 pmut_unpaired=pmut_unpaired,
                 vmut_paired=vmut_paired,
                 vmut_unpaired=vmut_unpaired,
                 force=force,
                 parallel=parallel,
                 max_procs=max_procs)
    ends_mod.run(ct_file=ct_file,
                 end3_fmean=end3_fmean,
                 insert_fmean=insert_fmean,
                 ends_var=ends_var,
                 force=force,
                 parallel=parallel,
                 max_procs=max_procs)
    clusts_mod.run(ct_file=ct_file,
                   clust_conc=clust_conc,
                   force=force,
                   parallel=parallel,
                   max_procs=max_procs)


params = merge_params(clusts_mod.params,
                      ends_mod.params,
                      muts_mod.params)


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate parameter files. """
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
