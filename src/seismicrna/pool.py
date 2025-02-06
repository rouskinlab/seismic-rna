from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

from click import command

from .core import path
from .core.arg import (CMD_POOL,
                       arg_input_path,
                       opt_pooled,
                       opt_relate_pos_table,
                       opt_relate_read_table,
                       opt_verify_times,
                       opt_tmp_pfx,
                       opt_keep_tmp,
                       opt_max_procs,
                       opt_force)
from .core.error import InconsistentValueError
from .core.logs import logger
from .core.report import BranchF, AncestorsF
from .core.run import run_func
from .core.task import dispatch
from .core.tmp import release_to_out, with_tmp_dir
from .core.write import need_write
from .relate.dataset import PoolDataset, load_relate_dataset
from .relate.report import PoolReport, RelateReport
from .relate.table import RelateDatasetTabulator
from .table import tabulate


def write_report(out_dir: Path, **kwargs):
    report = PoolReport(ended=datetime.now(), **kwargs)
    return report.save(out_dir, force=True)


@with_tmp_dir(pass_keep_tmp=False)
def pool_samples(out_dir: Path,
                 name: str,
                 branches: Iterable[str],
                 ref: str,
                 samples: Iterable[str], *,
                 tmp_dir: Path,
                 relate_pos_table: bool,
                 relate_read_table: bool,
                 verify_times: bool,
                 n_procs: int,
                 force: bool):
    """ Pool one or more samples (vertically).

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    name: str
        Name of the pool.
    branches: Iterable[str]
        Branches of the datasets being pooled.
    ref: str
        Name of the reference
    samples: Iterable[str]
        Names of the samples in the pool.
    tmp_dir: Path
        Temporary directory.
    relate_pos_table: bool
        Tabulate relationships per position for relate data.
    relate_read_table: bool
        Tabulate relationships per read for relate data
    verify_times: bool
        Verify that report files from later steps have later timestamps.
    n_procs: bool
        Number of processors to use.
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Pool report file.
    """
    began = datetime.now()
    # Ensure that branches is not an exhaustible generator.
    branches = list(branches)
    # Deduplicate and sort the samples.
    sample_counts = Counter(samples)
    if max(sample_counts.values()) > 1:
        logger.warning(f"Pool {repr(name)} with reference {repr(ref)} "
                       f"in {out_dir} got duplicate samples: {sample_counts}")
    samples = sorted(sample_counts)
    # Determine the output report file.
    pool_report_file = PoolReport.build_path({path.TOP: out_dir,
                                              path.SAMPLE: name,
                                              path.BRANCHES: branches,
                                              path.REF: ref})
    if need_write(pool_report_file, force):
        # Because Relate and Pool report files have the same name, it
        # would be possible to overwrite a Relate report with a Pool
        # report, rendering the Relate dataset unusable; prevent this.
        if pool_report_file.is_file():
            # Check whether the report file contains a Pool report.
            try:
                PoolReport.load(pool_report_file)
            except ValueError:
                # The report file does not contain a Pool report.
                raise TypeError(f"Cannot overwrite {pool_report_file} with "
                                f"{PoolReport.__name__}: would cause data loss")
        # To be able to load, the pooled dataset must have access to the
        # original relate dataset(s) in the temporary directory.
        report_kwargs = dict(sample=name,
                             ref=ref,
                             pooled_samples=samples,
                             began=began)
        for sample in samples:
            relate_dataset = load_relate_dataset(
                RelateReport.build_path({path.TOP: out_dir,
                                         path.SAMPLE: sample,
                                         path.BRANCHES: branches,
                                         path.REF: ref}),
                verify_times=verify_times
            )
            # Ensure all datasets have the same branch and ancestors;
            # it's possible for them to have the same branches but
            # different individual branch or ancestors attributes.
            if BranchF.key in report_kwargs:
                assert AncestorsF.key in report_kwargs
                if relate_dataset.branch != report_kwargs[BranchF.key]:
                    raise InconsistentValueError(
                        "Datasets have different values for branch: "
                        f"{report_kwargs[BranchF.key]} ≠ "
                        f"{repr(relate_dataset.branch)}"
                    )
                assert relate_dataset.ancestors == report_kwargs[AncestorsF.key]
            else:
                assert AncestorsF.key not in report_kwargs
                report_kwargs[BranchF.key] = relate_dataset.branch
                report_kwargs[AncestorsF.key] = relate_dataset.ancestors
            # Make links to the files for this dataset in the temporary
            # directory.
            relate_dataset.link_data_dirs_to_tmp(tmp_dir)
        if BranchF.key not in report_kwargs:
            assert AncestorsF.key not in report_kwargs
            # Default to these values if there are no samples.
            report_kwargs[BranchF.key] = ""
            report_kwargs[AncestorsF.key] = branches
        # Tabulate the pooled dataset.
        pool_dataset = load_relate_dataset(write_report(tmp_dir,
                                                        **report_kwargs),
                                           verify_times=verify_times)
        tabulate(pool_dataset,
                 RelateDatasetTabulator,
                 pos_table=relate_pos_table,
                 read_table=relate_read_table,
                 clust_table=False,
                 n_procs=n_procs,
                 force=True)
        # Rewrite the report file with the updated time.
        release_to_out(out_dir,
                       tmp_dir,
                       write_report(tmp_dir, **report_kwargs).parent)
    return pool_report_file.parent


@run_func(CMD_POOL)
def run(input_path: Iterable[str | Path], *,
        pooled: str,
        relate_pos_table: bool,
        relate_read_table: bool,
        verify_times: bool,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        max_procs: int,
        force: bool) -> list[Path]:
    """ Merge samples (vertically) from the Relate step. """
    if not pooled:
        raise ValueError("No name for the pooled sample was given via --pooled")
    # Group the datasets by output directory, branches, and reference.
    pools = defaultdict(list)
    for dataset in load_relate_dataset.iterate(input_path,
                                               verify_times=verify_times):
        # Check whether the dataset was pooled.
        if isinstance(dataset, PoolDataset):
            # If so, then use all samples in the pool.
            samples = dataset.samples
        else:
            # Otherwise, use just the sample of the dataset.
            samples = [dataset.sample]
        key = dataset.top, tuple(dataset.branches), dataset.ref
        pools[key].extend(samples)
    # Make each pool of samples.
    return dispatch(pool_samples,
                    max_procs=max_procs,
                    pass_n_procs=True,
                    args=[(out_dir, pooled, branches, ref, samples)
                          for (out_dir, branches, ref), samples
                          in pools.items()],
                    kwargs=dict(relate_pos_table=relate_pos_table,
                                relate_read_table=relate_read_table,
                                verify_times=verify_times,
                                tmp_pfx=tmp_pfx,
                                keep_tmp=keep_tmp,
                                force=force))


params = [
    arg_input_path,
    opt_pooled,
    opt_relate_pos_table,
    opt_relate_read_table,
    opt_verify_times,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_max_procs,
    opt_force,
]


@command(CMD_POOL, params=params)
def cli(*args, **kwargs):
    """ Merge samples (vertically) from the Relate step. """
    return run(*args, **kwargs)
