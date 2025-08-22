from typing import Iterable

import numpy as np

from .em import EMRun
from ..core.error import NoDataError, InconsistentValueError
from ..core.logs import logger
from ..core.mu import (calc_sum_arcsine_distance,
                       calc_mean_arcsine_distance,
                       calc_pearson)
from ..core.validate import require_atleast, require_equal, require_array_equal

NOCONV = 0


def get_common_k(runs: list[EMRun]):
    """ Find the number of clusters (k) from among EM clustering runs.
    If there are multiple ks, then raise a ValueError. """
    ks: list[int] = sorted({run.k for run in runs})
    if not ks:
        raise NoDataError("Got 0 numbers of clusters")
    if len(ks) != 1:
        raise InconsistentValueError(f"Got > 1 number of clusters: {ks}")
    return ks[0]


def sort_runs(runs: list[EMRun]):
    """ Sort the runs of EM clustering by decreasing likelihood so that
    the run with the best (largest) likelihood comes first. """
    if runs:
        # Verify that every run has the same k; if not, the likelihood
        # is not directly comparable between runs.
        get_common_k(runs)
    return sorted(runs, key=lambda run: run.log_like, reverse=True)


class EMRunsK(object):
    """ One or more EM runs with the same number of clusters. """

    def __init__(self,
                 runs: list[EMRun],
                 max_pearson_run: float,
                 min_marcd_run: float,
                 max_arcd_vs_ens_avg: float,
                 max_gini_run: float,
                 max_jackpot_quotient: float,
                 max_loglike_vs_best: float,
                 min_pearson_vs_best: float,
                 max_marcd_vs_best: float):
        if not runs:
            raise NoDataError("Got no EM runs")
        # Sort the runs from largest to smallest likelihood.
        runs = sort_runs(runs)
        # Flag runs that fail to meet the filters.
        self.max_pearson_run = max_pearson_run
        self.min_marcd_run = min_marcd_run
        self.max_arcd_vs_ens_avg = max_arcd_vs_ens_avg
        self.max_gini_run = max_gini_run
        self.max_jackpot_quotient = max_jackpot_quotient
        # Set the criteria for whether this number of clusters passes.
        self.max_loglike_vs_best = max_loglike_vs_best
        self.min_pearson_vs_best = min_pearson_vs_best
        self.max_marcd_vs_best = max_marcd_vs_best
        # Check whether each run shows signs of being overclustered.
        # To select only the valid runs, use "not" with the opposite of
        # the desired inequality because runs with just one cluster will
        # produce NaN values, which should always compare as True here.
        self.run_not_overclustered = np.array(
            [not ((run.max_pearson > max_pearson_run
                   and run.min_marcd < min_marcd_run)
                  or
                  (run.max_arcd_vs_ens_avg > max_arcd_vs_ens_avg
                   and run.max_gini > max_gini_run))
             for run in runs]
        )
        # Check whether each run shows signs of being underclustered.
        # To select only the valid runs, use "not" with the opposite of
        # the desired inequality because runs for which the jackpotting
        # quotient was not calculated will produce NaN values, which
        # should always compare as True here.
        self.run_not_underclustered = np.array(
            [not (run.jackpot_quotient > max_jackpot_quotient)
             for run in runs]
        )
        # Select the best run.
        self.best = runs[self.best_index()]
        # Number of runs.
        self.n_runs_total = len(runs)
        # Number of clusters (K).
        self.k = get_common_k(runs)
        # Number of iterations until convergenge for each run.
        self.converged = np.array([run.iter if run.converged else NOCONV
                                   for run in runs])
        # Log-likelihood of each run.
        self.log_likes = np.array([run.log_like for run in runs])
        # BIC of each run.
        self.bics = np.array([run.bic for run in runs])
        # Jackpotting quotient of each run.
        self.jackpot_quotients = np.array([run.jackpot_quotient
                                           for run in runs])
        # Maximum ARCD between the clusters and average.
        self.max_arcds_vs_ens_avg = np.array([run.max_arcd_vs_ens_avg
                                              for run in runs])
        # Maximum Gini difference between the clusters and average.
        self.max_ginis = np.array([run.max_gini for run in runs])
        # Minimum MARCD between any two clusters in each run.
        self.min_marcds = np.array([run.min_marcd for run in runs])
        # Maximum correlation between any two clusters in each run.
        self.max_pearsons = np.array([run.max_pearson for run in runs])
        # MARCD between each run and the best run.
        self.marcds_vs_best = np.array(
            [calc_mean_arcsine_distance_clusters(run.mus.values,
                                                 self.best.mus.values)
             for run in runs]
        )
        # Correlation between each run and the best run.
        self.pearsons_vs_best = np.array(
            [calc_mean_pearson_clusters(run.mus.values,
                                        self.best.mus.values)
             for run in runs]
        )

    def run_passing(self, allow_underclustered: bool = False):
        """ Whether each run passed the filters. """
        if allow_underclustered:
            return self.run_not_overclustered
        return self.run_not_overclustered & self.run_not_underclustered

    def n_runs_passing(self, **kwargs):
        """ Number of runs passing the filters. """
        return int(np.count_nonzero(self.run_passing(**kwargs)))

    def get_valid_index(self, i: int | list[int] | np.ndarray, **kwargs):
        """ Index(es) of valid run number(s) `i`. """
        return np.flatnonzero(self.run_passing(**kwargs))[i]

    def best_index(self, **kwargs) -> int:
        """ Index of the best valid run. """
        try:
            # The best run is the valid run with the largest likelihood.
            return self.get_valid_index(0, **kwargs)
        except IndexError:
            # If no runs are valid, then use the best invalid run.
            return 0

    def subopt_indexes(self, **kwargs):
        """ Indexes of the valid suboptimal runs. """
        return self.get_valid_index(np.arange(1, self.n_runs_passing(**kwargs)),
                                    **kwargs)

    def loglike_vs_best(self, **kwargs):
        """ Log likelihood difference between the best and second-best
        runs. """
        try:
            index1, index2 = self.get_valid_index([0, 1], **kwargs)
            return float(self.log_likes[index1] - self.log_likes[index2])
        except IndexError:
            return np.nan

    def pearson_vs_best(self, **kwargs):
        """ Maximum Pearson correlation between the best run and any
        other run. """
        try:
            return float(np.max(
                self.pearsons_vs_best[self.subopt_indexes(**kwargs)]
            ))
        except ValueError:
            return np.nan

    def marcd_vs_best(self, **kwargs):
        """ Minimum MARCD between the best run and any other run. """
        try:
            return float(np.min(
                self.marcds_vs_best[self.subopt_indexes(**kwargs)]
            ))
        except ValueError:
            return np.nan

    def _n_min_runs_passing(self, **kwargs):
        n_runs_passing = self.n_runs_passing(**kwargs)
        min_runs_passing = min(self.n_runs_total, 2)
        return n_runs_passing, min_runs_passing

    def enough_runs_passing(self, **kwargs):
        """ Whether enough runs passed. """
        n_runs_passing, min_runs_passing = self._n_min_runs_passing(**kwargs)
        return n_runs_passing >= min_runs_passing

    def passing(self, **kwargs):
        """ Whether this number of clusters passes the filters. """
        n_runs_passing, min_runs_passing = self._n_min_runs_passing(**kwargs)
        if n_runs_passing < min_runs_passing:
            logger.detail(f"{self} did not pass: {n_runs_passing} runs passed, "
                          f"but needed {min_runs_passing}")
            return False
        # Make sure that if any attribute is NaN, the run will still be
        # able to pass; this can be done by requiring each inequality
        # to be True in order to not pass (since < and > will be False
        # if one side is NaN).
        loglike_vs_best = self.loglike_vs_best(**kwargs)
        if loglike_vs_best > self.max_loglike_vs_best > 0.:
            logger.detail(f"{self} did not pass: difference between 1st/2nd "
                          f"log likelihoods is {loglike_vs_best}, but needed "
                          f"to be ≤ {self.max_loglike_vs_best}")
            return False
        pearson_vs_best = self.pearson_vs_best(**kwargs)
        if pearson_vs_best < self.min_pearson_vs_best:
            logger.detail(f"{self} did not pass: Pearson correlation between "
                          f"best run and any other run is {pearson_vs_best}, "
                          f"but needed to be ≥ {self.min_pearson_vs_best}")
            return False
        marcd_vs_best = self.marcd_vs_best(**kwargs)
        if marcd_vs_best > self.max_marcd_vs_best:
            logger.detail(f"{self} did not pass: MARCD between best run and "
                          f"any other run is {marcd_vs_best}, but needed to "
                          f"be ≤ {self.max_marcd_vs_best}")
            return False
        logger.detail(f"{self} passed all filters using {kwargs}")
        return True

    def summarize(self, **kwargs):
        """ Summarize the results of the runs. """
        lines = [f"EM runs for K={self.k}",
                 "\nPARAMETERS\n"]
        for attr in ["max_pearson_run",
                     "min_marcd_run",
                     "max_jackpot_quotient",
                     "max_loglike_vs_best",
                     "max_arcd_vs_ens_avg",
                     "max_gini_run",
                     "min_pearson_vs_best",
                     "max_marcd_vs_best"]:
            lines.append(f"{attr} = {getattr(self, attr)}")
        lines.append("\nRUNS\n")
        for attr in ["n_runs_total",
                     "converged",
                     "log_likes",
                     "bics",
                     "jackpot_quotients",
                     "max_arcds_vs_ens_avg",
                     "max_ginis",
                     "min_marcds",
                     "max_pearsons",
                     "marcds_vs_best",
                     "pearsons_vs_best",
                     "run_not_overclustered",
                     "run_not_underclustered"]:
            lines.append(f"{attr} = {getattr(self, attr)}")
        lines.append("\nPASSING\n")
        for attr in ["run_passing",
                     "n_runs_passing",
                     "best_index",
                     "loglike_vs_best",
                     "pearson_vs_best",
                     "marcd_vs_best",
                     "enough_runs_passing",
                     "passing"]:
            func = getattr(self, attr)
            lines.append(f"{attr} = {func(**kwargs)}")
        return "\n".join(lines)

    def __str__(self):
        return (f"{type(self).__name__} with {self.k} clusters "
                f"and {self.n_runs_total} runs")


def find_best_k(ks: Iterable[EMRunsK], **kwargs):
    """ Find the best number of clusters. """
    # Sort the runs by increasing numbers of clusters.
    ks = sorted(ks, key=lambda runs: runs.k)
    if not ks:
        logger.warning("No numbers of clusters exist")
        return 0
    # Select only the numbers of clusters that pass the filters.
    # For the largest number of clusters, underclustering can be allowed
    # to permit this number to be identified as the best number so far;
    # for all other numbers, use all filters.
    ks_passing = [runs for runs in ks[:-1] if runs.passing()]
    if ks[-1].passing(**kwargs):
        ks_passing.append(ks[-1])
    if not ks_passing:
        logger.warning("No numbers of clusters pass the filters")
        return 0
    # Of the remaining numbers of clusters, find the number that gives
    # the smallest BIC.
    ks_passing = sorted(ks_passing, key=lambda runs: runs.best.bic)
    # Return that number of clusters.
    return ks_passing[0].k


def assign_clusterings(mus1: np.ndarray, mus2: np.ndarray):
    """ Optimally assign clusters from two groups to each other. """
    n1, k1 = mus1.shape
    n2, k2 = mus2.shape
    require_equal("mus1.shape[0]", n1, n2, "mus2.shape[0]")
    require_equal("mus1.shape[1]", k1, k2, "mus2.shape[1]")
    if n1 >= 1 and k1 >= 1:
        # Match the clusters using linear_sum_assignment.
        costs = np.array([[calc_sum_arcsine_distance(mus1[:, cluster1],
                                                     mus2[:, cluster2])
                           for cluster2 in range(k2)]
                          for cluster1 in range(k1)]).reshape((k1, k2))
        from scipy.optimize import linear_sum_assignment
        rows, cols = linear_sum_assignment(costs)
    else:
        # If n1 == 0, then the costs matrix will contain NaN, which will
        # cause linear_sum_assignment to raise an error.
        rows = np.arange(k1)
        cols = np.arange(k1)
    assert np.array_equal(rows, np.arange(k1))
    assert rows.shape == cols.shape
    return rows, cols


def _concat_clusters(mus: np.ndarray, clusters: np.ndarray):
    n, k = mus.shape
    require_equal("clusters.shape", clusters.shape, (k,), "(mus.shape[1],)")
    require_atleast("mus.shape[1]", k, 1)
    require_array_equal(
        "sorted(clusters)", np.sort(clusters), np.arange(k), "np.arange(k)"
    )
    return np.concatenate([mus[:, c] for c in clusters])


def _assign_concat_clusters(mus1: np.ndarray, mus2: np.ndarray):
    rows, cols = assign_clusterings(mus1, mus2)
    return _concat_clusters(mus1, rows), _concat_clusters(mus2, cols)


def calc_mean_arcsine_distance_clusters(mus1: np.ndarray, mus2: np.ndarray):
    """ MARCD between all clusters after matching. """
    return float(calc_mean_arcsine_distance(*_assign_concat_clusters(mus1,
                                                                     mus2)))


def calc_mean_pearson_clusters(mus1: np.ndarray, mus2: np.ndarray):
    """ Pearson correlation between all clusters after matching. """
    return float(calc_pearson(*_assign_concat_clusters(mus1, mus2)))
