import os
from logging import getLogger
from pathlib import Path

from click import command

from .clusts import load_pclust
from .ends import load_pends
from .muts import load_pmut
from ..core import path
from ..core.arg import (opt_param_dir,
                        opt_profile_name,
                        opt_sample,
                        opt_paired_end,
                        opt_read_length,
                        opt_reverse_fraction,
                        opt_min_mut_gap,
                        opt_num_reads,
                        opt_batch_size,
                        opt_brotli_level,
                        opt_force,
                        opt_parallel,
                        opt_max_procs)
from ..core.rna import find_ct_section
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..relate.sim import simulate_relate

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


def get_param_dir_fields(param_dir: Path):
    fields = path.parse(param_dir, path.RefSeg, path.SectSeg)
    params_dir = fields[path.TOP]
    if params_dir.name != path.SIM_PARAM_DIR:
        raise ValueError(
            f"Expected parameter directory named {repr(path.SIM_PARAM_DIR)}, "
            f"but got {repr(params_dir.name)}"
        )
    return params_dir.parent, fields[path.REF], fields[path.SECT]


def load_param_dir(param_dir: Path, profile: str):
    """ Load all parameters for a profile in a directory. """
    prefix = param_dir.joinpath(profile)
    section = find_ct_section(prefix.with_suffix(path.CT_EXT))
    pmut = load_pmut(prefix.with_suffix(path.PARAM_MUTS_EXT))
    u5s, u3s, pends = load_pends(prefix.with_suffix(path.PARAM_ENDS_EXT))
    pclust = load_pclust(prefix.with_suffix(path.PARAM_CLUSTS_EXT))
    return section, pmut, u5s, u3s, pends, pclust


def from_param_dir(param_dir: Path,
                   profile: str,
                   min_mut_gap: int,
                   **kwargs):
    """ Simulate a Relate dataset given parameter files. """
    sim_dir, _, _ = get_param_dir_fields(param_dir)
    section, pmut, u5s, u3s, pends, pclust = load_param_dir(param_dir, profile)
    return simulate_relate(out_dir=sim_dir.joinpath(path.SIM_SAMPLES_DIR),
                           ref=section.ref,
                           refseq=section.seq,
                           pmut=pmut,
                           uniq_end5s=u5s,
                           uniq_end3s=u3s,
                           pends=pends,
                           pclust=pclust,
                           min_mut_gap=min_mut_gap,
                           **kwargs)


@run_func(logger.critical, with_tmp=True)
def run(*,
        param_dir: tuple[str, ...],
        profile_name: str,
        sample: str,
        paired_end: bool,
        read_length: int,
        reverse_fraction: float,
        min_mut_gap: int,
        num_reads: int,
        batch_size: int,
        brotli_level: int,
        tmp_dir: Path,
        force: bool,
        parallel: bool,
        max_procs: int):
    """ Simulate a Relate dataset. """
    return dispatch(from_param_dir,
                    max_procs=max_procs,
                    parallel=parallel,
                    pass_n_procs=False,
                    args=as_list_of_tuples(map(Path, param_dir)),
                    kwargs=dict(sample=sample,
                                profile=profile_name,
                                paired=paired_end,
                                read_length=read_length,
                                p_rev=reverse_fraction,
                                min_mut_gap=min_mut_gap,
                                num_reads=num_reads,
                                batch_size=batch_size,
                                brotli_level=brotli_level,
                                tmp_dir=tmp_dir,
                                force=force))


params = [
    opt_param_dir,
    opt_profile_name,
    opt_sample,
    opt_paired_end,
    opt_read_length,
    opt_reverse_fraction,
    opt_min_mut_gap,
    opt_num_reads,
    opt_batch_size,
    opt_brotli_level,
    opt_force,
    opt_parallel,
    opt_max_procs
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate a Relate dataset. """
    run(*args, **kwargs)
