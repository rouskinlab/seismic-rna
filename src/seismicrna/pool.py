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
                       opt_num_cpus,
                       opt_force)
from .core.error import InconsistentValueError, NoDataError
from .core.logs import logger
from .core.report import BranchesF
from .core.run import run_func
from .core.task import dispatch
from .core.tmp import release_to_out, with_tmp_dir
from .core.write import need_write
from .relate.dataset import PoolMutsDataset, load_relate_dataset
from .relate.report import PoolReport, RelateReport
from .relate.table import RelateDatasetTabulator
from .table import tabulate


def write_report(out_dir: Path, **kwargs):
    report = PoolReport(ended=datetime.now(), **kwargs)
    return report.save(out_dir, force=True)


@with_tmp_dir(pass_keep_tmp=False)
def pool_samples(out_dir: Path,
                 name: str,
                 branches_flat: Iterable[str],
                 ref: str,
                 samples: Iterable[str], *,
                 tmp_dir: Path,
                 relate_pos_table: bool,
                 relate_read_table: bool,
                 verify_times: bool,
                 num_cpus: int,
                 force: bool):
    """ Pool one or more samples (vertically).

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    name: str
        Name of the pooled sample.
    branches_flat: Iterable[str]
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
    num_cpus: bool
        Number of processors to use.
    force: bool
        Force the report to be written, even if it exists.

    Returns
    -------
    pathlib.Path
        Path of the Pool report file.
    """
    began = datetime.now()
    branches_flat = list(branches_flat)
    # Deduplicate and sort the samples.
    sample_counts = Counter(samples)
    if max(sample_counts.values()) > 1:
        logger.warning(f"Pool {repr(name)} with reference {repr(ref)} "
                       f"in {out_dir} got duplicate samples: {sample_counts}")
    samples = sorted(sample_counts)
    # Determine the output report file.
    pool_report_file = PoolReport.build_path({path.TOP: out_dir,
                                              path.SAMPLE: name,
                                              path.BRANCHES: branches_flat,
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
                                         path.BRANCHES: branches_flat,
                                         path.REF: ref}),
                verify_times=verify_times
            )
            # Ensure all datasets have the same branches; it's possible
            # for two datasets to have different branches despite having
            # the same flattened branches.
            if BranchesF.key in report_kwargs:
                if relate_dataset.branches != report_kwargs[BranchesF.key]:
                    raise InconsistentValueError(
                        "Cannot pool datasets with different branches: "
                        f"{report_kwargs[BranchesF.key]} ≠ "
                        f"{relate_dataset.branches}"
                    )
            else:
                report_kwargs[BranchesF.key] = relate_dataset.branches
            # Make links to the files for this dataset in the temporary
            # directory.
            relate_dataset.link_data_dirs_to_tmp(tmp_dir)
        if BranchesF.key not in report_kwargs:
            raise NoDataError(
                f"No samples were given to make pooled sample {repr(name)} "
                f"with branches {branches_flat} and reference {repr(ref)} "
                f"in {out_dir}"
            )
        # Tabulate the pooled dataset.
        pool_dataset = load_relate_dataset(write_report(tmp_dir,
                                                        **report_kwargs),
                                           verify_times=verify_times)
        tabulate(pool_dataset,
                 RelateDatasetTabulator,
                 pos_table=relate_pos_table,
                 read_table=relate_read_table,
                 clust_table=False,
                 num_cpus=num_cpus,
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
        num_cpus: int,
        force: bool) -> list[Path]:
    """ Merge samples (vertically) from the Relate step. """
    if not pooled:
        raise ValueError("No name for the pooled sample was given via --pooled")
    # Group the datasets by output directory, branches, and reference.
    pools = defaultdict(list)
    for dataset in load_relate_dataset.iterate(input_path,
                                               verify_times=verify_times):
        # Check whether the dataset was pooled.
        if isinstance(dataset, PoolMutsDataset):
            # If so, then use all samples in the pool.
            samples = dataset.samples
        else:
            # Otherwise, use just the sample of the dataset.
            samples = [dataset.sample]
        # Flatten the dict of branches because the path preserves
        # only the flat structure, and this step must group datasets
        # that will have the same pooled path.
        branches_flat = tuple(path.flatten_branches(dataset.branches))
        key = (dataset.top,
               branches_flat,
               dataset.ref)
        pools[key].extend(samples)
        logger.detail(f"Added samples {samples} for {dataset}")
    # Make each pool of samples.
    return dispatch(pool_samples,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=[(out_dir, pooled, branches_flat, ref, samples)
                          for (out_dir, branches_flat, ref), samples
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
    opt_num_cpus,
    opt_force,
]


@command(CMD_POOL, params=params)
def cli(*args, **kwargs):
    """ Merge samples (vertically) from the Relate step. """
    return run(*args, **kwargs)
