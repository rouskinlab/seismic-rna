from pathlib import Path
from typing import Iterable

from click import command

from .cleanfa import clean_fasta
from ..core import path
from ..core.arg import (CMD_CLEANFA,
                        arg_input_path,
                        opt_inplace,
                        opt_out_dir,
                        opt_force,
                        opt_num_cpus)
from ..core.run import run_func
from ..core.task import dispatch


@run_func(CMD_CLEANFA)
def run(input_path: Iterable[str | Path], *,
        inplace: bool,
        out_dir: str | Path,
        force: bool,
        num_cpus: int):
    """ Clean the names and sequences in FASTA files. """
    # List all the files to clean.
    input_files = list(path.find_files_chain(input_path, [path.FastaSeg]))
    # Determine the output path of each cleaned file.
    if inplace:
        # If modifying in-place, then the output and input paths match.
        output_files = input_files
    else:
        # Otherwise, determine the output paths using out_dir.
        out_dir = path.sanitize(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        output_files = path.transpaths(out_dir, input_files, strict=True)
    # Generate the positional arguments for clean_fasta.
    args = list(zip(input_files, output_files, strict=True))
    # Clean the files; if modifying in-place, force must be True.
    return dispatch(clean_fasta,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=args,
                    kwargs=dict(force=force or inplace))


params = [
    arg_input_path,
    opt_inplace,
    opt_out_dir,
    opt_force,
    opt_num_cpus,
]


@command(CMD_CLEANFA, params=params)
def cli(*args, **kwargs):
    """ Clean the names and sequences in FASTA files. """
    return run(*args, **kwargs)
