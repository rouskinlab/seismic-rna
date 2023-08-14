from logging import getLogger
from math import inf
from pathlib import Path

from .em import EmClustering
from .compare import RunOrderResults, find_best_order, sort_replicate_runs
from .report import ClustReport
from .write import write_batches, write_log_counts, write_mus, write_props
from ..core.bitvect import BitMonolith, UniqMutBits
from ..core.parallel import dispatch
from ..mask.load import MaskLoader

logger = getLogger(__name__)


def cluster(call_report: Path, max_order: int, n_runs: int, *,
            min_iter: int, max_iter: int, conv_thresh: float, min_muts: int,
            n_procs: int, rerun: bool):
    """ Run all processes of clustering reads from one filter. """
    # Load the vector calling report.
    loader = MaskLoader.open(Path(call_report))
    # Check if the clustering report file already exists.
    report_file = ClustReport.build_path(loader.out_dir,
                                         sample=loader.sample,
                                         ref=loader.ref,
                                         sect=loader.sect)
    if rerun or not report_file.is_file():
        logger.info(f"Began clustering {loader} up to order {max_order} "
                    f"cluster(s) and {n_runs} independent run(s) per order")
        # Get the unique bit vectors.
        uniq_muts = BitMonolith(loader.section,
                                loader.iter_batches_processed()).get_uniq_muts()
        # Run EM clustering for every number of clusters.
        results = run_max_order(loader, uniq_muts, max_order, n_runs,
                                min_iter=min_iter, max_iter=max_iter,
                                conv_thresh=conv_thresh, n_procs=n_procs)
        logger.info(
            f"Ended clustering {loader}: {find_best_order(results)} clusters")
        # Output the observed and expected counts for every best run.
        write_log_counts(results)
        # Output the cluster memberships in batches of reads.
        checksums = write_batches(results)
        report = ClustReport.from_clusters(results, loader, uniq_muts,
                                           max_order, n_runs,
                                           min_iter=min_iter,
                                           max_iter=max_iter,
                                           conv_thresh=conv_thresh,
                                           checksums=checksums)
        report.save()
    else:
        logger.warning(f"File exists: {report_file}")
    # Return the path of the clustering report file.
    return report_file


def run_max_order(loader: MaskLoader,
                  uniq_muts: UniqMutBits,
                  max_order: int,
                  n_runs: int, *,
                  min_iter: int,
                  max_iter: int,
                  conv_thresh: float,
                  n_procs: int):
    """
    Find the optimal order (i.e. number of clusters) for EM, up to max_order.
    """
    if n_runs < 1:
        logger.warning(f"Number of EM runs must be ≥ 1: setting to 1")
        n_runs = 1
    logger.info(f"Began clustering {loader} up to order {max_order} "
                f"with {n_runs} runs per order")
    # For each order, keep the best run and a summary of all runs.
    results: dict[int, RunOrderResults] = dict()
    # Run clustering for each number of clusters, up to max_clusters.
    while len(results) < max_order:
        # Determine the current order for clustering.
        order = len(results) + 1
        # Run EM clustering n_runs times with different starting points.
        runs = run_order(loader, uniq_muts, order,
                         n_runs=(n_runs if order > 1 else 1),
                         conv_thresh=(conv_thresh if order > 1 else inf),
                         min_iter=(min_iter * order if order > 1 else 2),
                         max_iter=(max_iter * order if order > 1 else 2),
                         n_procs=n_procs)
        # Output tables of the mutation rates and cluster proportions
        # for every run.
        for rank, run in enumerate(runs):
            write_mus(run, rank)
            write_props(run, rank)
        # Compute a summary of all runs.
        results[order] = RunOrderResults(runs)
        # Compare the BIC for this order to lower orders, if order > 1.
        best_bic = results[order].best.bic
        logger.debug(f"The best BIC with {order} cluster(s) is {best_bic}")
        if order > 1:
            prev_bic = results[order - 1].best.bic
            # Check if the best BIC is better (smaller) than the best
            # BIC from clustering with one cluster fewer.
            if best_bic < prev_bic:
                # If the best BIC is < the previous best BIC, then the
                # current model is the best.
                logger.debug(f"The BIC decreased from {prev_bic} "
                             f"(order {order - 1}) to {best_bic} "
                             f"(order {order})")
            else:
                # If the best BIC is ≥ than the previous best BIC, then
                # this model is worse than the previous model. Stop.
                logger.info(f"The BIC failed to decrease from {prev_bic} "
                            f"(order {order - 1}) to {best_bic} "
                            f"(order {order}): stopping")
                break
    return results


def run_order(loader: MaskLoader,
              uniq_muts: UniqMutBits,
              order: int,
              n_runs: int, *,
              min_iter: int,
              max_iter: int,
              conv_thresh: float,
              n_procs: int) -> list[EmClustering]:
    """ Run EM with a specific number of clusters. """
    logger.info(f"Began {n_runs} run(s) of EM with {order} cluster(s)")
    # Initialize one EmClustering object for each replicate run.
    runs = [EmClustering(loader, uniq_muts, order, conv_thresh=conv_thresh,
                         min_iter=min_iter, max_iter=max_iter)
            for _ in range(n_runs)]
    # Run independent replicates of the clustering algorithm.
    runs = dispatch([rep.run for rep in runs], n_procs,
                    parallel=True, pass_n_procs=False)
    logger.info(f"Ended {n_runs} run(s) of EM with {order} cluster(s)")
    return sort_replicate_runs(runs)
