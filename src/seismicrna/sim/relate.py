import os
from pathlib import Path
from typing import Iterable

import numpy as np
from click import command

from .clusts import load_pclust
from .ends import load_pends
from .muts import load_pmut
from ..core import path
from ..core.arg import (opt_param_dir,
                        opt_profile_name,
                        opt_sample_sim,
                        opt_branch,
                        opt_paired_end,
                        opt_read_length,
                        opt_reverse_fraction,
                        opt_probe,
                        opt_min_mut_gap,
                        opt_mut_collisions,
                        opt_mut_probs,
                        opt_num_reads,
                        opt_batch_size,
                        opt_write_read_names,
                        opt_brotli_level,
                        opt_force,
                        opt_num_cpus,
                        opt_seed)
from ..core.rna import find_ct_region
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..mask.main import set_mut_gap_params
from ..relate.sim import simulate_relate

COMMAND = __name__.split(os.path.extsep)[-1]


def _get_param_dir_fields(param_dir: Path):
    fields = path.parse(param_dir, [path.RefSeg, path.RegSeg])
    params_dir = fields[path.TOP]
    if params_dir.name != path.SIM_PARAM_DIR:
        raise ValueError(
            f"Expected parameter directory named {repr(path.SIM_PARAM_DIR)}, "
            f"but got {repr(params_dir.name)}"
        )
    return params_dir.parent, fields[path.REF], fields[path.REG]


def _load_param_dir(param_dir: Path, profile: str):
    """ Load all parameters for a profile in a directory. """
    prefix = param_dir.joinpath(profile)
    region = find_ct_region(prefix.with_suffix(path.CT_EXT))
    pmut = load_pmut(prefix.with_suffix(path.PARAM_MUTS_EXT))
    u5s, u3s, pends = load_pends(prefix.with_suffix(path.PARAM_ENDS_EXT))
    pclust = load_pclust(prefix.with_suffix(path.PARAM_CLUSTS_EXT))
    return region, pmut, u5s, u3s, pends, pclust


def _from_param_dir(param_dir: Path,
                   profile: str,
                   **kwargs):
    """ Simulate a Relate dataset given parameter files. """
    sim_dir, _, _ = _get_param_dir_fields(param_dir)
    region, pmut, u5s, u3s, pends, pclust = _load_param_dir(param_dir, profile)
    return simulate_relate(out_dir=sim_dir.joinpath(path.SIM_SAMPLES_DIR),
                           ref=region.ref,
                           refseq=region.seq,
                           pmut=pmut,
                           uniq_end5s=u5s,
                           uniq_end3s=u3s,
                           pends=pends,
                           pclust=pclust,
                           **kwargs)


@run_func(COMMAND, with_tmp=True)
def run(*,
        param_dir: Iterable[str | Path],
        profile_name: str,
        sample: str,
        branch: str,
        paired_end: bool,
        read_length: int,
        reverse_fraction: float,
        probe: str,
        min_mut_gap: int | None,
        mut_collisions: str,
        mut_probs: str | None,
        num_reads: int,
        batch_size: int,
        write_read_names: bool,
        brotli_level: int,
        tmp_dir: Path,
        force: bool,
        num_cpus: int,
        seed: int | None):
    """ Simulate a Relate dataset. """
    min_mut_gap, mut_collisions = set_mut_gap_params(probe,
                                                     min_mut_gap,
                                                     mut_collisions)
    mut_probs_arr = (np.array(list(map(float, mut_probs.split(","))), dtype=float)
                     if mut_probs is not None else None)
    return dispatch(_from_param_dir,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(map(Path, param_dir)),
                    kwargs=dict(sample=sample,
                                branch=branch,
                                profile=profile_name,
                                paired=paired_end,
                                read_length=read_length,
                                p_rev=reverse_fraction,
                                min_mut_gap=min_mut_gap,
                                mut_probs=mut_probs_arr,
                                mut_collisions=mut_collisions,
                                num_reads=num_reads,
                                batch_size=batch_size,
                                write_read_names=write_read_names,
                                brotli_level=brotli_level,
                                tmp_dir=tmp_dir,
                                force=force,
                                seed=seed))


params = [
    opt_param_dir,
    opt_profile_name,
    opt_sample_sim,
    opt_branch,
    opt_paired_end,
    opt_read_length,
    opt_reverse_fraction,
    opt_probe,
    opt_min_mut_gap,
    opt_mut_collisions,
    opt_mut_probs,
    opt_num_reads,
    opt_batch_size,
    opt_write_read_names,
    opt_brotli_level,
    opt_force,
    opt_num_cpus,
    opt_seed,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate a Relate dataset. """
    run(*args, **kwargs)
