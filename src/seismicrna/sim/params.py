import os
from pathlib import Path
from typing import Iterable

from click import command

from . import (clusts as clusts_mod,
               ends as ends_mod,
               muts as muts_mod)
from ..core.arg import merge_params
from ..core.run import run_func

COMMAND = __name__.split(os.path.extsep)[-1]


@run_func(COMMAND)
def run(*,
        ct_file: Iterable[str | Path],
        pmut_paired: Iterable[tuple[str, float]],
        pmut_unpaired: Iterable[tuple[str, float]],
        vmut_paired: float,
        vmut_unpaired: float,
        center_fmean: float,
        center_fvar: float,
        length_fmean: float,
        length_fvar: float,
        clust_conc: float,
        force: bool,
        num_cpus: int):
    """ Simulate parameter files. """
    # Since ct_file is used three times, ensure it is not an exhaustible
    # generator.
    ct_file = list(ct_file)
    muts_mod.run(ct_file=ct_file,
                 pmut_paired=pmut_paired,
                 pmut_unpaired=pmut_unpaired,
                 vmut_paired=vmut_paired,
                 vmut_unpaired=vmut_unpaired,
                 force=force,
                 num_cpus=num_cpus)
    ends_mod.run(ct_file=ct_file,
                 center_fmean=center_fmean,
                 center_fvar=center_fvar,
                 length_fmean=length_fmean,
                 length_fvar=length_fvar,
                 force=force,
                 num_cpus=num_cpus)
    clusts_mod.run(ct_file=ct_file,
                   clust_conc=clust_conc,
                   force=force,
                   num_cpus=num_cpus)


params = merge_params(clusts_mod.params,
                      ends_mod.params,
                      muts_mod.params)


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate parameter files. """
    return run(*args, **kwargs)
