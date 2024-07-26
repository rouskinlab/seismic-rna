from datetime import datetime
from logging import getLogger
from math import inf
from pathlib import Path
from typing import Iterable

import numpy as np

from .compare import EMRunsK, find_best_k, sort_runs
from .csv import write_log_counts, write_mus, write_props
from .em import EMRun
from .io import write_batches
from .report import ClusterReport
from .uniq import UniqReads
from ..core.header import validate_ks
from ..core.io import recast_file_path
from ..core.logs import exc_info
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import release_to_out
from ..core.types import get_max_uint
from ..core.write import need_write
from ..mask.data import load_mask_dataset
from ..mask.report import MaskReport

logger = getLogger(__name__)

SEED_DTYPE = np.uint32


def run_k(uniq_reads: UniqReads,
          k: int,
          em_runs: int, *,
          n_procs: int,
          **kwargs) -> list[EMRun]:
    """ Run EM with a specific number of clusters. """
    if k < 1:
        raise ValueError(f"k must be ≥ 1, but got {k}")
    if em_runs < 1:
        logger.warning(
            f"Expected em_runs to be ≥ 1, but got {em_runs}: setting to 1"
        )
        em_runs = 1
    rng = np.random.default_rng()
    logger.info(f"Began {em_runs} run(s) of EM with {k} cluster(s)")
    # On some but not all platforms, using this central source of seeds
    # is necessary because otherwise the runs would all take identical
    # trajectories, defeating the purpose of replicates.
    seeds = as_list_of_tuples(rng.integers(get_max_uint(SEED_DTYPE),
                                           size=em_runs,
                                           dtype=SEED_DTYPE))
    runs = list(dispatch([EMRun(uniq_reads, k, **kwargs).run
                          for _ in range(em_runs)],
                         n_procs,
                         parallel=True,
                         pass_n_procs=False,
                         args=seeds))
    if len(runs) < em_runs:
        logger.warning(f"Obtained only {len(runs)} (of {em_runs}) "
                       f"run(s) of {uniq_reads} with {k} cluster(s)")
    return sort_runs(runs)


def run_ks(uniq_reads: UniqReads,
           ks: Iterable[int],
           em_runs: int, *,
           try_all_ks: bool,
           min_nrmsd_run: float,
           max_pearson_run: float,
           max_loglike_vs_best: float,
           min_pearson_vs_best: float,
           max_nrmsd_vs_best: float,
           min_iter: int,
           max_iter: int,
           em_thresh: float,
           n_procs: int,
           top: Path,
           **kwargs):
    """ Run EM with multiple numbers of clusters. """
    path_kwargs = dict(top=top,
                       sample=uniq_reads.sample,
                       ref=uniq_reads.ref,
                       sect=uniq_reads.section.name)
    runs_ks = dict()
    ks = validate_ks(ks)
    for k in ks:
        try:
            # Cluster em_runs times with different starting points.
            runs = run_k(uniq_reads,
                         k,
                         em_runs=(em_runs if k > 1 else 1),
                         em_thresh=(em_thresh if k > 1 else inf),
                         min_iter=(min_iter * k if k > 1 else 2),
                         max_iter=(max_iter * k if k > 1 else 2),
                         n_procs=n_procs,
                         **kwargs)
            logger.info(f"Ended {em_runs} run(s) of EM with {k} cluster(s)")
            # Output each run's mutation rates and cluster proportions.
            for rank, run in enumerate(runs):
                write_mus(run, rank=rank, **path_kwargs)
                write_props(run, rank=rank, **path_kwargs)
            # Collect all runs for this number of clusters.
            runs_ks[k] = EMRunsK(runs,
                                 max_pearson_run=max_pearson_run,
                                 min_nrmsd_run=min_nrmsd_run)
            # Check whether this number of clusters passes the filters.
            runs_ks[k].set_passing(max_loglike_vs_best=max_loglike_vs_best,
                                   min_pearson_vs_best=min_pearson_vs_best,
                                   max_nrmsd_vs_best=max_nrmsd_vs_best)
            if not (try_all_ks or k == find_best_k(runs_ks.values())):
                # The current k is not the best so far.
                break
        except Exception as error:
            logger.error(
                f"Failed to split {uniq_reads} into {k} cluster(s): {error}",
                exc_info=exc_info()
            )
    return runs_ks


def cluster(mask_report_file: Path, *,
            min_clusters: int,
            max_clusters: int,
            try_all_ks: bool,
            write_all_ks: bool,
            em_runs: int,
            n_procs: int,
            brotli_level: int,
            force: bool,
            tmp_dir: Path,
            **kwargs):
    """ Cluster unique reads from one mask dataset. """
    # Check if the cluster report file already exists.
    cluster_report_file = recast_file_path(mask_report_file,
                                           MaskReport,
                                           ClusterReport)
    if need_write(cluster_report_file, force):
        began = datetime.now()
        logger.info(f"Began clustering {mask_report_file}")
        # Load the unique reads.
        dataset = load_mask_dataset(mask_report_file)
        if dataset.min_mut_gap != 3:
            logger.warning("For clustering, it is highly recommended to use "
                           "the observer bias correction with min_mut_gap=3, "
                           f"but got min_mut_gap={dataset.min_mut_gap}")
        uniq_reads = UniqReads.from_dataset_contig(dataset)
        # Run clustering for every number of clusters.
        if max_clusters >= 1:
            max_clusters_use = max_clusters
        elif max_clusters == 0:
            if try_all_ks:
                # Prevent accidentally forcing the use of a huge number
                # of clusters.
                raise ValueError("If using --try-all-ks, "
                                 "then must specify --max-clusters ≥ 1")
            # Set the maximum number of clusters to more than the number
            # of reads: effectively no limit.
            max_clusters_use = dataset.num_reads + 1
        else:
            raise ValueError(
                f"max_clusters must be ≥ 0, but got {max_clusters}"
            )
        runs_ks = run_ks(uniq_reads,
                         ks=range(min_clusters,
                                  max_clusters_use + 1),
                         try_all_ks=try_all_ks,
                         em_runs=em_runs,
                         n_procs=n_procs,
                         top=tmp_dir,
                         **kwargs)
        # Choose which numbers of clusters to write.
        if write_all_ks:
            write_ks = list(runs_ks.values())
        elif (best_k := find_best_k(runs_ks.values())) >= 1:
            write_ks = [runs_ks[best_k]]
        else:
            write_ks = []
        # Write the observed and expected counts for every best run.
        write_log_counts(uniq_reads,
                         write_ks,
                         top=tmp_dir,
                         sample=dataset.sample,
                         ref=dataset.ref,
                         sect=dataset.sect)
        # Output the cluster memberships in batches of reads.
        checksums, ks_written = write_batches(dataset,
                                              write_ks,
                                              brotli_level,
                                              tmp_dir)
        # Write the cluster report.
        ended = datetime.now()
        report = ClusterReport.from_clusters(list(runs_ks.values()),
                                             uniq_reads,
                                             min_clusters=min_clusters,
                                             max_clusters=max_clusters,
                                             try_all_ks=try_all_ks,
                                             write_all_ks=write_all_ks,
                                             ks_written=ks_written,
                                             em_runs=em_runs,
                                             checksums=checksums,
                                             began=began,
                                             ended=ended,
                                             **kwargs)
        report_saved = report.save(tmp_dir)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
        logger.info(f"Ended clustering {mask_report_file}")
    return cluster_report_file

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
