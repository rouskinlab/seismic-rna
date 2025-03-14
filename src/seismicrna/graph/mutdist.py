import os
from functools import cached_property

import numpy as np
import pandas as pd
from click import command

from .cgroup import make_tracks
from .color import ColorMapGraph, RelColorMap
from .dataset import DatasetGraph, DatasetWriter, DatasetRunner
from .trace import HIST_COUNT_NAME, get_hist_trace
from ..core.arg import opt_mutdist_null
from ..core.dataset import MutsDataset
from ..core.header import REL_NAME, make_header
from ..core.logs import logger
from ..core.rel import RelPattern
from ..core.run import log_command
from ..core.seq import FIELD_END5, FIELD_END3
from ..core.table import PositionTable, all_patterns
from ..core.task import as_list_of_tuples, dispatch
from ..core.unbias import (calc_p_noclose_given_ends_auto,
                           calc_p_ends_observed,
                           triu_dot)
from ..table import get_tabulator_type

COMMAND = __name__.split(os.path.extsep)[-1]
NULL_SUFFIX = "-NULL"


def get_null_name(name: str):
    return f"{name}{NULL_SUFFIX}"


def _get_num_bins(dataset: MutsDataset):
    return dataset.region.length


def _init_hists(dataset: MutsDataset, rel_name: str):
    header = make_header(rels=[rel_name], ks=dataset.ks)
    return pd.DataFrame(0. if header.get_is_clustered() else 0,
                        pd.RangeIndex(_get_num_bins(dataset),
                                      name=HIST_COUNT_NAME),
                        header.index)


def _calc_hists(batch_num: int, *,
                dataset: MutsDataset,
                pattern: RelPattern,
                rel_name: str,
                count_pos_ends: bool):
    """ Calculate the histogram of the smallest distances between two
    mutations in a read. """
    batch = dataset.get_batch(batch_num)
    hists = _init_hists(dataset, rel_name)
    min_mut_dist = batch.calc_min_mut_dist(pattern)
    num_bins = _get_num_bins(dataset)
    if dataset.is_clustered:
        if not isinstance(batch.read_weights, pd.DataFrame):
            raise TypeError(batch.read_weights)
        for (k, clust), weights in batch.read_weights.items():
            col = (rel_name, k, clust)
            hists.loc[:, col] += np.bincount(min_mut_dist,
                                             weights=weights,
                                             minlength=num_bins)
    else:
        if batch.read_weights is not None:
            raise TypeError(batch.read_weights)
        hists.loc[:, rel_name] += np.bincount(min_mut_dist,
                                              minlength=num_bins)
    if count_pos_ends:
        # Also count the relationships per position and 5'/3' ends.
        batch_counts = batch.count_all(patterns=all_patterns(dataset.pattern),
                                       ks=dataset.ks,
                                       count_ends=True,
                                       count_pos=True,
                                       count_read=False)
    else:
        batch_counts = None
    return batch.read_lengths.max(initial=0), hists, batch_counts


class MutationDistanceGraph(DatasetGraph, ColorMapGraph):
    """ Distance between the closest two mutations in each read. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Distance between closest"

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def __init__(self, *, mutdist_null: bool, **kwargs):
        super().__init__(**kwargs)
        self.calc_null = mutdist_null

    @property
    def x_title(self):
        return "Distance (nt)"

    @property
    def y_title(self):
        return "Number of reads"

    @cached_property
    def row_tracks(self):
        return make_tracks(self.dataset, self.k, self.clust)

    @cached_property
    def _data(self):
        logger.routine(f"Began calculating real histogram for {self}")
        results = dispatch(
            _calc_hists,
            num_cpus=self.num_cpus,
            pass_num_cpus=False,
            as_list=False,
            ordered=False,
            raise_on_error=True,
            args=as_list_of_tuples(range(self.dataset.num_batches)),
            kwargs=dict(dataset=self.dataset,
                        pattern=self.pattern,
                        rel_name=self.rel_name,
                        count_pos_ends=self.calc_null)
        )
        max_read_length = 0
        hists = _init_hists(self.dataset, self.rel_name)
        counts = list()
        for (batch_max_read_length, batch_hists, batch_counts) in results:
            if batch_max_read_length > max_read_length:
                max_read_length = batch_max_read_length
            hists += batch_hists
            if self.calc_null:
                assert batch_counts is not None
                counts.append(batch_counts)
        if self.calc_null:
            # Calculate the mutation rates and 5'/3' ends.
            dataset_type = type(self.dataset)
            tabulator_type = get_tabulator_type(dataset_type, count=True)
            init_keywords = (set(get_tabulator_type(dataset_type).init_kws())
                             - {"get_batch_count_all", "num_batches"})
            kwargs = {kw: getattr(self.dataset, kw) for kw in init_keywords}
            tabulator = tabulator_type(batch_counts=counts,
                                       count_ends=True,
                                       count_pos=True,
                                       count_read=False,
                                       **kwargs)
        else:
            tabulator = None
        logger.routine(f"Ended calculating real histogram for {self}")
        return max_read_length, hists, tabulator

    @property
    def max_read_length(self):
        max_read_length, hists, tabulator = self._data
        return max_read_length

    @property
    def hists(self):
        max_read_length, hists, tabulator = self._data
        return hists

    @property
    def tabulator(self):
        max_read_length, hists, tabulator = self._data
        if tabulator is None:
            raise TypeError("Cannot tabulate if --mutdist-null is False")
        return tabulator

    @cached_property
    def table(self):
        tables = list(self.tabulator.generate_tables(pos=True,
                                                     read=False,
                                                     clust=False))
        assert len(tables) == 1
        table = tables[0]
        assert isinstance(table, PositionTable)
        return table

    @cached_property
    def loc_clusters(self):
        return (self.k if self.k is not None else slice(None),
                self.clust if self.clust is not None else slice(self.clust))

    @cached_property
    def _real_hist(self):
        if self.dataset.is_clustered:
            cols = (self.rel_name,) + self.loc_clusters
            return self.hists.loc[:, cols]
        return self.hists

    @cached_property
    def _null_hist(self):
        logger.routine(f"Began calculating null histogram for {self}")
        if self.dataset.is_clustered:
            end_counts = self.tabulator.end_counts.loc[:, self.loc_clusters]
            num_reads = self.tabulator.num_reads.loc[self.loc_clusters].values
        else:
            end_counts = self.tabulator.end_counts
            num_reads = self.tabulator.num_reads
        end5 = self.dataset.region.end5
        end5s = np.asarray(end_counts.index.get_level_values(FIELD_END5) - end5)
        end3s = np.asarray(end_counts.index.get_level_values(FIELD_END3) - end5)
        p_ends = np.atleast_3d(calc_p_ends_observed(self.dataset.region.length,
                                                    end5s,
                                                    end3s,
                                                    end_counts.values))
        p_mut = self.table.fetch_ratio(k=self.k, rel=self.rel_name)
        # For every possible gap, calculate the fraction of reads that
        # have no two mutations closer than that gap.
        p_noclose_gap = np.empty((self.max_read_length,
                                  self._real_hist.columns.size),
                                 dtype=float)
        if p_noclose_gap.size > 0:
            p_noclose_gap[0] = 1.
        logger.detail("Calculating null fraction of reads in which every pair "
                      "of mutations would have at least N bases between them, "
                      f"from N = 1 to {self.max_read_length - 1}")
        for gap in range(1, self.max_read_length):
            p_noclose_ends = calc_p_noclose_given_ends_auto(p_mut.values, gap)
            p_noclose_gap[gap] = triu_dot(p_noclose_ends, p_ends)
            logger.detail(f"N = {gap}: {p_noclose_gap[gap]}")
        # For every possible distance, calculate the fraction of reads
        # where the closest two mutations have exactly that distance,
        # and for 0 the fraction of reads with fewer than two mutations.
        p_dist = np.zeros_like(self._real_hist, dtype=float)
        if p_dist.size > 0:
            p_dist[0] = p_noclose_gap[self.max_read_length - 1]
        p_dist[1: self.max_read_length] = -np.diff(p_noclose_gap, axis=0)
        assert np.all(p_dist >= 0.)
        assert np.allclose(p_dist.sum(axis=0), 1.)
        # Multiply by the number of reads to obtain the histogram.
        logger.routine(f"Ended calculating null histogram for {self}")
        return pd.DataFrame(
            p_dist * num_reads,
            index=self._real_hist.index,
            columns=self._real_hist.columns
        ).rename(columns=get_null_name, level=REL_NAME, copy=False)

    @cached_property
    def g_test(self):
        """ G-test statistic and P-value. """
        if self.calc_null:
            observed = self._real_hist.values
            expected = self._null_hist.values
            assert observed.ndim == 2
            assert observed.shape == expected.shape
            assert np.allclose(observed.sum(axis=0), expected.sum(axis=0))
            n, k = observed.shape
            dof = n - 1
            if dof >= 1:
                with np.errstate(divide="ignore", invalid="ignore"):
                    g_stat = 2. * np.where(
                        np.logical_and(observed > 0, expected > 0),
                        observed * np.log(observed / expected),
                        0.
                    ).sum(axis=0)
                from scipy.stats import chi2
                p_value = 1. - chi2.cdf(g_stat, dof)
            else:
                g_stat = np.zeros(k, dtype=float)
                p_value = np.ones_like(g_stat)
        else:
            g_stat = np.nan
            p_value = np.nan
        return (pd.Series(g_stat, index=self._real_hist.columns),
                pd.Series(p_value, index=self._real_hist.columns))

    @cached_property
    def data(self):
        if self.calc_null:
            return pd.concat([self._real_hist, self._null_hist], axis=1)
        return self._real_hist

    def get_traces(self):
        if self.calc_null:
            for row, ((_, real_values), (_, null_values)) in enumerate(
                    zip(self._real_hist.items(),
                        self._null_hist.items(),
                        strict=True),
                    start=1
            ):
                yield (row, 1), get_hist_trace(real_values,
                                               self.rel_name,
                                               self.cmap)
                yield (row, 1), get_hist_trace(null_values,
                                               get_null_name(self.rel_name),
                                               self.cmap)
        else:
            for row, (_, real_values) in enumerate(
                    self._real_hist.items(),
                    start=1
            ):
                yield (row, 1), get_hist_trace(real_values,
                                               self.rel_name,
                                               self.cmap)


class MutationDistanceWriter(DatasetWriter):

    def get_graph(self, rel, **kwargs):
        return MutationDistanceGraph(dataset=self.dataset,
                                     rel=rel,
                                     **kwargs)


class MutationDistanceRunner(DatasetRunner):

    @classmethod
    def get_writer_type(cls):
        return MutationDistanceWriter

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_mutdist_null]

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=MutationDistanceRunner.params())
def cli(*args, **kwargs):
    """ Distance between the closest two mutations in each read. """
    return MutationDistanceRunner.run(*args, **kwargs)
