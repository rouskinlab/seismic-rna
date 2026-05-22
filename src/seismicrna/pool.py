from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Iterable

from click import command

from .core import path
from .core.arg import (CMD_POOL,
                       arg_input_path,
                       arg_pooled_sample,
                       opt_min_pearson_pool,
                       opt_max_marcd_pool,
                       opt_idmut_pos_table,
                       opt_idmut_read_table,
                       opt_verify_times,
                       opt_tmp_pfx,
                       opt_keep_tmp,
                       opt_num_cpus,
                       opt_force)
from .core.error import InconsistentValueError, NoDataError
from .core.logs import logger
from .core.mu.compare import calc_mean_arcsine_distance, calc_pearson
from .core.report import BranchesF
from .core.run import run_func
from .core.table import INFOR_REL, MUTAT_REL
from .core.task import dispatch
from .core.tmp import release_to_out, with_tmp_dir
from .core.write import need_write
from .idmut.dataset import PoolMutsDataset, load_idmut_dataset
from .idmut.report import PoolReport, IDmutReport
from .idmut.table import IDmutDatasetTabulator, IDmutPositionTableLoader
from .table import tabulate


def _calc_sample_mus(idmut_dataset, num_cpus: int):
    """ Get per-position mutation rates for an IDmut dataset.

    Uses an existing position table CSV if available (fast), otherwise
    computes from batch files via IDmutDatasetTabulator (slow).
    """
    table_files = list(IDmutPositionTableLoader.find_tables([idmut_dataset.dir]))
    if len(table_files) > 1:
        raise ValueError(f"Expected at most 1 position table in {idmut_dataset.dir}, "
                         f"but found {len(table_files)}: {table_files}")
    if table_files:
        table = IDmutPositionTableLoader(table_files[0])
        return table.fetch_ratio(rel=MUTAT_REL, squeeze=True)
    tabulator = IDmutDatasetTabulator(
        dataset=idmut_dataset, count_pos=True, count_read=False, num_cpus=num_cpus
    )
    data = tabulator.data_per_pos
    return data[MUTAT_REL] / data[INFOR_REL]


def write_report(out_dir: Path, **kwargs):
    report = PoolReport(ended=datetime.now(), **kwargs)
    return report.save(out_dir, force=True)


@with_tmp_dir(pass_keep_tmp=False)
def pool_samples(out_dir: Path,
                 pooled_sample: str,
                 branches_flat: Iterable[str],
                 ref: str,
                 samples: Iterable[str], *,
                 tmp_dir: Path,
                 min_pearson: float,
                 max_marcd: float,
                 idmut_pos_table: bool,
                 idmut_read_table: bool,
                 verify_times: bool,
                 num_cpus: int,
                 force: bool):
    """ Pool one or more samples (vertically).

    Parameters
    ----------
    out_dir: pathlib.Path
        Output directory.
    pooled_sample: str
        Name of the pooled sample.
    branches_flat: Iterable[str]
        Branches of the datasets being pooled.
    ref: str
        Name of the reference
    samples: Iterable[str]
        Names of the samples in the pool.
    tmp_dir: Path
        Temporary directory.
    min_pearson: float
        Skip pooling if any pair of samples has Pearson r below this.
    max_marcd: float
        Skip pooling if any pair of samples has MARCD above this.
    idmut_pos_table: bool
        Tabulate relationships per position for IDmut data.
    idmut_read_table: bool
        Tabulate relationships per read for IDmut data.
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
        logger.warning(f"Pool {repr(pooled_sample)} with reference {repr(ref)} "
                       f"in {out_dir} got duplicate samples: {sample_counts}")
    samples = sorted(sample_counts)
    # Determine the output report file.
    pool_report_file = PoolReport.build_path({path.TOP: out_dir,
                                              path.SAMPLE: pooled_sample,
                                              path.BRANCHES: branches_flat,
                                              path.REF: ref})
    if need_write(pool_report_file, force):
        # Because IDmut and Pool report files have the same name, it
        # would be possible to overwrite an IDmut report with a Pool
        # report, rendering the IDmut dataset unusable; prevent this.
        if pool_report_file.is_file():
            # Check whether the report file contains a Pool report.
            try:
                PoolReport.load(pool_report_file)
            except ValueError:
                # The report file does not contain a Pool report.
                raise TypeError(f"Cannot overwrite {pool_report_file} with "
                                f"{PoolReport.__name__}: would cause data loss")
        # To be able to load, the pooled dataset must have access to the
        # original IDmut dataset(s) in the temporary directory.
        report_kwargs = dict(sample=pooled_sample,
                             ref=ref,
                             pooled_samples=samples,
                             min_pearson=min_pearson,
                             max_marcd=max_marcd,
                             began=began)
        sample_mus = {}
        for sample in samples:
            idmut_dataset = load_idmut_dataset(
                IDmutReport.build_path({path.TOP: out_dir,
                                         path.SAMPLE: sample,
                                         path.BRANCHES: branches_flat,
                                         path.REF: ref}),
                verify_times=verify_times
            )
            # Ensure all datasets have the same branches; it's possible
            # for two datasets to have different branches despite having
            # the same flattened branches.
            if BranchesF.key in report_kwargs:
                if idmut_dataset.branches != report_kwargs[BranchesF.key]:
                    raise InconsistentValueError(
                        "Cannot pool datasets with different branches: "
                        f"{report_kwargs[BranchesF.key]} ≠ "
                        f"{idmut_dataset.branches}"
                    )
            else:
                report_kwargs[BranchesF.key] = idmut_dataset.branches
            # Make links to the files for this dataset in the temporary
            # directory.
            idmut_dataset.link_data_dirs_to_tmp(tmp_dir)
            # Calculate the mutation rates for the sample.
            sample_mus[sample] = _calc_sample_mus(idmut_dataset, num_cpus)
        if BranchesF.key not in report_kwargs:
            raise NoDataError(
                f"No samples were given to make pooled sample {repr(pooled_sample)} "
                f"with branches {branches_flat} and reference {repr(ref)} "
                f"in {out_dir}"
            )
        # Pairwise similarity filter.
        for s1, s2 in combinations(samples, 2):
            mu1 = sample_mus[s1]
            mu2 = sample_mus[s2]
            pearson = float(calc_pearson(mu1, mu2))
            if pearson < min_pearson:
                logger.warning(
                    f"Skipping pool {repr(pooled_sample)} with reference {repr(ref)} "
                    f"in {out_dir}: Pearson r = {pearson:.4f} between "
                    f"{repr(s1)} and {repr(s2)} is less than "
                    f"min_pearson = {min_pearson}"
                )
                return None
            marcd = float(calc_mean_arcsine_distance(mu1, mu2))
            if marcd > max_marcd:
                logger.warning(
                    f"Skipping pool {repr(pooled_sample)} with reference {repr(ref)} "
                    f"in {out_dir}: MARCD = {marcd:.4f} between "
                    f"{repr(s1)} and {repr(s2)} exceeds "
                    f"max_marcd = {max_marcd}"
                )
                return None
        # Tabulate the pooled dataset.
        pool_dataset = load_idmut_dataset(write_report(tmp_dir,
                                                        **report_kwargs),
                                           verify_times=verify_times)
        tabulate(pool_dataset,
                 IDmutDatasetTabulator,
                 pos_table=idmut_pos_table,
                 read_table=idmut_read_table,
                 clust_table=False,
                 num_cpus=num_cpus,
                 force=True)
        # Rewrite the report file with the updated time.
        release_to_out(out_dir,
                       tmp_dir,
                       write_report(tmp_dir, **report_kwargs).parent)
    return pool_report_file.parent


@run_func(CMD_POOL)
def run(pooled_sample: str,
        input_path: Iterable[str | Path], *,
        idmut_pos_table: bool,
        idmut_read_table: bool,
        min_pearson: float,
        max_marcd: float,
        verify_times: bool,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        num_cpus: int,
        force: bool) -> list[Path]:
    """ Merge samples (vertically) from the IDmut step. """
    if not pooled_sample:
        raise ValueError("No name for the pooled sample was given")
    # Group the datasets by output directory, branches, and reference.
    pools = defaultdict(list)
    for dataset in load_idmut_dataset.iterate(input_path,
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
    # Make each pool of samples, dropping any that were skipped by the
    # pairwise similarity filter (pool_samples returns None for those).
    return [result
            for result in dispatch(pool_samples,
                                   num_cpus=num_cpus,
                                   pass_num_cpus=True,
                                   as_list=False,
                                   ordered=False,
                                   raise_on_error=False,
                                   args=[(out_dir, pooled_sample, branches_flat, ref, samples)
                                         for (out_dir, branches_flat, ref), samples
                                         in pools.items()],
                                   kwargs=dict(min_pearson=min_pearson,
                                               max_marcd=max_marcd,
                                               idmut_pos_table=idmut_pos_table,
                                               idmut_read_table=idmut_read_table,
                                               verify_times=verify_times,
                                               tmp_pfx=tmp_pfx,
                                               keep_tmp=keep_tmp,
                                               force=force))
            if result is not None]


params = [
    arg_pooled_sample,
    arg_input_path,
    opt_idmut_pos_table,
    opt_idmut_read_table,
    opt_min_pearson_pool,
    opt_max_marcd_pool,
    opt_verify_times,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_num_cpus,
    opt_force,
]


@command(CMD_POOL, params=params)
def cli(*args, **kwargs):
    """ Merge samples (vertically) from the IDmut step. """
    return run(*args, **kwargs)
