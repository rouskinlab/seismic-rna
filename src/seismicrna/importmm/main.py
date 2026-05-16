from pathlib import Path
from typing import Iterable

from click import command

from .write import import_mm
from ..core.arg import (CMD_IMPORTMM,
                        arg_input_path,
                        opt_out_dir,
                        opt_tmp_pfx,
                        opt_branch,
                        opt_importmm_sample,
                        opt_min_reads,
                        opt_batch_size,
                        opt_insert3,
                        opt_brotli_level,
                        opt_write_read_names,
                        opt_force,
                        opt_num_cpus)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch


@run_func(CMD_IMPORTMM, with_tmp=True, pass_keep_tmp=False)
def run(input_path: Iterable[str | Path], *,
        sample: str,
        out_dir: str | Path,
        tmp_dir: Path,
        branch: str,
        min_reads: int,
        batch_size: int,
        insert3: bool,
        write_read_names: bool,
        brotli_level: int,
        num_cpus: int,
        force: bool):
    """ Import RNA Framework Mutation Map (MM) files as relate outputs. """
    return dispatch(import_mm,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(list(input_path)),
                    kwargs=dict(sample=sample,
                                out_dir=Path(out_dir),
                                tmp_dir=tmp_dir,
                                branch=branch,
                                min_reads=min_reads,
                                batch_size=batch_size,
                                insert3=insert3,
                                write_read_names=write_read_names,
                                brotli_level=brotli_level,
                                force=force))


params = [
    # Input files
    arg_input_path,
    # Output directories
    opt_out_dir,
    opt_tmp_pfx,
    opt_branch,
    # Sample identification
    opt_importmm_sample,
    # Import options
    opt_min_reads,
    opt_batch_size,
    opt_insert3,
    opt_brotli_level,
    # Output options
    opt_write_read_names,
    # Parallelization
    opt_num_cpus,
    # File generation
    opt_force,
]


@command(CMD_IMPORTMM, params=params)
def cli(**kwargs):
    """ Import RNA Framework Mutation Map (MM) files as relate outputs. """
    return run(**kwargs)
