from functools import cached_property
from logging import getLogger
from pathlib import Path

import pandas as pd

from .names import ORD_CLS_NAME
from .report import ClustReport
from ..core import path
from ..core.bitvect import BitBatch, ClusterBitBatch, ClustBitCounter
from ..core.load import BatchChainLoader, no_kwargs
from ..core.mu import calc_mu_adj_df, calc_f_obs_df
from ..mask.load import MaskLoader

logger = getLogger(__name__)


def get_cluster_index(max_order: int):
    """ Get the index for clusters up to a given order. """
    return pd.MultiIndex.from_tuples((order, cluster)
                                     for order in range(1, max_order + 1)
                                     for cluster in range(1, order + 1))


class ClustLoader(BatchChainLoader):
    """ Load clustering results. """

    def __init__(self, report: ClustReport):
        super().__init__(report)
        self.best_order = report.best_order

    @classmethod
    def get_report_type(cls):
        return ClustReport

    @classmethod
    def get_import_type(cls):
        return MaskLoader

    @property
    def import_path_fields(self):
        return {**super().import_path_fields, path.SECT: self.sect}

    def get_refseq(self):
        return self.import_loader.get_refseq()

    @property
    def section(self):
        return self.import_loader.section

    @property
    def min_mut_gap(self) -> int:
        return self.import_loader.min_mut_gap

    @cached_property
    def clusters(self):
        """ Order and number of each cluster. """
        return get_cluster_index(self.best_order)

    @no_kwargs
    def load_data_personal(self, batch_file: Path):
        # Load the cluster memberships of the reads in the batch.
        data = pd.read_csv(batch_file,
                           index_col=[0],
                           header=list(range(len(ORD_CLS_NAME))))
        # Cast the cluster numbers from strings to integers.
        data.columns = pd.MultiIndex.from_tuples([(int(order), int(k))
                                                  for order, k in data.columns],
                                                 names=ORD_CLS_NAME)
        return data

    def iter_batches_personal(self):
        yield from super().iter_batches_personal()

    @no_kwargs
    def process_batch(self, imported_batch: BitBatch,
                      personal_batch: pd.DataFrame):
        return ClusterBitBatch(imported_batch.section,
                               imported_batch.info,
                               imported_batch.affi,
                               personal_batch)

    def iter_batches_processed(self, **kwargs):
        yield from super().iter_batches_processed(**kwargs)

    @cached_property
    def resps(self):
        """ Cluster memberships. """
        try:
            # If any batches exist, then stack them into one DataFrame.
            return pd.concat(self.iter_batches_personal(), axis=0)
        except ValueError:
            # Otherwise, return an empty DataFrame.
            return pd.DataFrame(index=[], columns=self.clusters, dtype=float)

    @cached_property
    def mus(self):
        """ Mutation rates, adjusted for observer bias. """
        # Count the informative and affirmative bits at each position
        # in each cluster.
        counts = ClustBitCounter(self.section,
                                 self.clusters,
                                 self.iter_batches_processed())
        # Return the adjusted mutation rates per position.
        return calc_mu_adj_df(counts.f_affi_per_pos,
                              self.section,
                              self.min_mut_gap)

    @cached_property
    def f_obs(self):
        """ Fraction of all reads in each cluster that would actually be
        observed, as a Series with dimension (clusters). """
        return calc_f_obs_df(self.mus, self.section, self.min_mut_gap)

    @cached_property
    def n_reads_obs(self):
        """ Observed number of reads in each cluster. """
        return self.resps.sum(axis=0)

    @cached_property
    def n_reads_adj(self):
        """ Adjusted number of reads in each cluster. """
        return self.n_reads_obs / self.f_obs
