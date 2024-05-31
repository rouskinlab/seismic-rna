import os
from pathlib import Path

from click import command

from .muts import (sim_pmut,
                   get_paired,
                   make_pmut_means_paired,
                   make_pmut_means_unpaired)
from ..core.seq import DNA
from ..relate.sim import simulate_relate

import os
from pathlib import Path

import numpy as np
import pandas as pd
from click import command

from .clusts import load_pclust
from .ends import load_pends
from .muts import load_pmut
from ..core import path
from ..core.arg import (docdef,
                        opt_param_dir,
                        opt_profile_name,
                        opt_sample,
                        opt_num_reads,
                        opt_batch_size,
                        opt_brotli_level,
                        opt_force,
                        opt_parallel,
                        opt_max_procs)
from ..core.header import index_order_clusts
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rna import find_ct_section
from ..core.write import need_write

COMMAND = __name__.split(os.path.extsep)[-1]


def test_simulate():
    sample = "mysample"
    num_reads = 2 ** 16
    batch_size = 2 ** 16
    ref = "myref"
    seq = DNA.random(2000)
    nclust = 2
    paired = get_paired(seq, nclust)
    pclust = sim_pclust(nclust)

    u5s, u3s, pends = sim_pends(1, len(seq), len(seq) * 0.8, 250, 0.1)
    pm = make_pmut_means_paired()
    um = make_pmut_means_unpaired()

    pmut = [sim_pmut(paired[cluster], pm, um, 0.001, 0.04) for cluster in paired.columns]

    out_dir = Path.cwd().joinpath("out")
    return simulate_relate(out_dir, sample, ref, seq, batch_size, num_reads, pmut, u5s, u3s, pends, pclust.values,
                           brotli_level=10, force=True)


def get_param_dir_fields(param_dir: Path):
    fields = path.parse(param_dir, path.RefSeg, path.SectSeg)
    return fields[path.TOP], fields[path.REF], fields[path.SECT]


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
                   **kwargs):
    """ Simulate a Relate dataset given parameter files. """
    sim_dir, _, _ = get_param_dir_fields(param_dir)
    section, pmut, u5s, u3s, pends, pclust = load_param_dir(param_dir, profile)
    return simulate_relate(out_dir=sim_dir,
                           ref=section.ref,
                           refseq=section.seq,
                           pmut=pmut,
                           uniq_end5s=u5s,
                           uniq_end3s=u3s,
                           pends=pends,
                           pclust=pclust,
                           **kwargs)


@docdef.auto()
def run(param_dir: tuple[str, ...],
        profile_name: str,
        sample: str,
        num_reads: int,
        batch_size: int,
        brotli_level: int,
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
                                num_reads=num_reads,
                                batch_size=batch_size,
                                brotli_level=brotli_level,
                                force=force))


params = [
    opt_param_dir,
    opt_profile_name,
    opt_sample,
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
