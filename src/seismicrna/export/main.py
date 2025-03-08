from collections import defaultdict
from pathlib import Path
from typing import Iterable

from click import command

from .meta import parse_refs_metadata, parse_samples_metadata
from .web import export_sample
from ..cluster.data import ClusterPositionTableLoader, ClusterAbundanceTableLoader
from ..core.arg import (CMD_EXPORT,
                        arg_input_path,
                        opt_samples_meta,
                        opt_refs_meta,
                        opt_all_pos,
                        opt_force,
                        opt_num_cpus)
from ..core.run import run_func
from ..core.task import dispatch
from ..mask.table import MaskPositionTableLoader, MaskReadTableLoader


@run_func(CMD_EXPORT)
def run(input_path: Iterable[str | Path], *,
        samples_meta: str,
        refs_meta: str,
        all_pos: bool,
        force: bool,
        num_cpus: int) -> list[Path]:
    """ Export a file of each sample for the seismic-graph web app. """
    tables = defaultdict(list)
    samples_metadata = (parse_samples_metadata(Path(samples_meta))
                        if samples_meta
                        else dict())
    refs_metadata = (parse_refs_metadata(Path(refs_meta))
                     if refs_meta
                     else dict())
    for table_type in [MaskPositionTableLoader,
                       MaskReadTableLoader,
                       ClusterPositionTableLoader,
                       ClusterAbundanceTableLoader]:
        for table in table_type.load_tables(input_path):
            tables[(table.top, table.sample)].append(table)
    return dispatch(export_sample,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=list(tables.items()),
                    kwargs=dict(samples_metadata=samples_metadata,
                                refs_metadata=refs_metadata,
                                all_pos=all_pos,
                                force=force))


params = [
    arg_input_path,
    opt_samples_meta,
    opt_refs_meta,
    opt_all_pos,
    opt_force,
    opt_num_cpus,
]


@command(CMD_EXPORT, params=params)
def cli(*args, **kwargs):
    """ Export each sample to SEISMICgraph (https://seismicrna.org). """
    return run(*args, **kwargs)
