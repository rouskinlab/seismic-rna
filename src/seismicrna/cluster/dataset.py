from abc import ABC
from functools import cached_property
from itertools import combinations, product
from typing import Any

import numpy as np
import pandas as pd

from ._table import (ClusterBatchTabulator,
                     ClusterPositionTableLoader,
                     ClusterPositionTableWriter,
                     ClusterAbundanceTableLoader,
                     ClusterAbundanceTableWriter)
from .batch import ClusterMutsBatch
from .io import ClusterBatchIO
from .report import ClusterReport, JoinClusterReport
from ..core import path
from ..core.batch import MutsBatch
from ..core.dataset import (LoadFunction,
                            Dataset,
                            LoadedDataset,
                            MultistepDataset,
                            UnbiasDataset,
                            MergedUnbiasDataset)
from ..core.header import (NUM_CLUSTS_NAME,
                           ClustHeader,
                           list_clusts,
                           list_ks_clusts,
                           validate_ks)
from ..core.join import (BATCH_NUM,
                         READ_NUMS,
                         SEG_END5S,
                         SEG_END3S,
                         MUTS,
                         RESPS,
                         JoinMutsDataset)
from ..core.report import JoinedClustersF, KsWrittenF, BestKF
from ..core.seq import POS_NAME, BASE_NAME
from ..core.table import MUTAT_REL
from ..mask.batch import MaskMutsBatch
from ..mask.dataset import load_mask_dataset


class ClusterDataset(Dataset, ABC):
    """ Dataset of clustered data. """


class ClusterReadDataset(ClusterDataset, LoadedDataset):
    """ Load clustering results. """

    @classmethod
    def get_report_type(cls):
        return ClusterReport

    @classmethod
    def get_batch_type(cls):
        return ClusterBatchIO

    @cached_property
    def ks(self):
        return validate_ks(self.report.get_field(KsWrittenF))

    @cached_property
    def best_k(self):
        return self.report.get_field(BestKF)

    @property
    def pattern(self):
        return None


class ClusterMutsDataset(ClusterDataset, MultistepDataset, UnbiasDataset):
    """ Merge cluster responsibilities with mutation data. """

    @classmethod
    def get_dataset1_load_func(cls):
        return load_mask_dataset

    @classmethod
    def get_dataset2_type(cls):
        return ClusterReadDataset

    @property
    def pattern(self):
        return self.data1.pattern

    @pattern.setter
    def pattern(self, pattern):
        self.data1.pattern = pattern

    @property
    def region(self):
        return self.data1.region

    @property
    def min_mut_gap(self):
        return getattr(self.data1, "min_mut_gap")

    @min_mut_gap.setter
    def min_mut_gap(self, min_mut_gap):
        self.data1.min_mut_gap = min_mut_gap

    @property
    def quick_unbias(self):
        return getattr(self.data1, "quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return getattr(self.data1, "quick_unbias_thresh")

    @property
    def ks(self):
        return getattr(self.data2, "ks")

    @property
    def best_k(self):
        return getattr(self.data2, "best_k")

    def _integrate(self, batch1: MaskMutsBatch, batch2: ClusterBatchIO):
        return ClusterMutsBatch(batch=batch1.batch,
                                region=batch1.region,
                                seg_end5s=batch1.seg_end5s,
                                seg_end3s=batch1.seg_end3s,
                                muts=batch1.muts,
                                resps=batch2.resps,
                                sanitize=False)


def _get_clust_params(dataset: ClusterMutsDataset, max_procs: int = 1):
    """ Get the mutation rates and proportion for each cluster. If table
    files already exist, then use them to get the parameters; otherwise,
    calculate the parameters from the dataset. """
    # Try to load the tables from files.
    path_fields = {path.TOP: dataset.top,
                   path.SAMP: dataset.sample,
                   path.REF: dataset.ref,
                   path.REG: dataset.region.name}
    pos_table_file = ClusterPositionTableLoader.build_path(**path_fields)
    try:
        pos_table = ClusterPositionTableLoader(pos_table_file)
    except FileNotFoundError:
        pos_table = None
    abundance_table_file = ClusterAbundanceTableLoader.build_path(**path_fields)
    try:
        abundance_table = ClusterAbundanceTableLoader(abundance_table_file)
    except FileNotFoundError:
        abundance_table = None
    # If either table file does not exist, then calculate the tables.
    if pos_table is None or abundance_table is None:
        def get_batch_count_all(batch_num: int, **kwargs):
            batch = dataset.get_batch(batch_num)
            return batch.count_all(**kwargs)

        tabulator = ClusterBatchTabulator(
            top=dataset.top,
            sample=dataset.sample,
            region=dataset.region,
            refseq=dataset.refseq,
            pattern=dataset.pattern,
            min_mut_gap=dataset.min_mut_gap,
            quick_unbias=dataset.quick_unbias,
            quick_unbias_thresh=dataset.quick_unbias_thresh,
            ks=dataset.ks,
            num_batches=dataset.num_batches,
            get_batch_count_all=get_batch_count_all,
            count_ends=True,
            count_pos=(pos_table is None),
            count_read=False,
            max_procs=max_procs
        )
        pos_table = ClusterPositionTableWriter(tabulator)
        abundance_table = ClusterAbundanceTableWriter(tabulator)
    # Calculate the parameters from the tables.
    mus = pos_table.fetch_ratio(rel=MUTAT_REL).loc[:, MUTAT_REL]
    pis = abundance_table.proportions
    # Merge the parameters into one DataFrame with the proportions as
    # the first row with the index (0, "p") and the mutation rates on
    # subsequent rows.
    pis_df = pis.to_frame().T
    assert mus.index.names == [POS_NAME, BASE_NAME]
    pis_df.index = pd.MultiIndex.from_tuples([(0, "p")], names=mus.index.names)
    assert pis.index.equals(mus.columns)
    params = pd.concat([pis_df, mus])
    return params


def _calc_diff_log_odds(a: np.ndarray | pd.Series | pd.DataFrame,
                        b: np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the log odds difference between probabilities a and b:
    log(a / (1-a)) - log(b / (1-b))

    Wherever a == b, the log odds difference will be 0 (even if a and b
    are 0 or 1). Otherwise, if a or b == 0 or 1, the return value will
    be inf or -inf.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(a != b,
                        np.log((a * (1. - b)) / (b * (1. - a))),
                        0.)


def _join_regions_k(region_params: dict[str, pd.DataFrame]):
    """ Determine the optimal way to join regions . """
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import csr_matrix
    # Validate the arguments.
    n = len(region_params)
    assert n >= 1
    dfs = iter(region_params.values())
    df = next(dfs)
    assert isinstance(df, pd.DataFrame)
    clusters = df.columns
    assert isinstance(clusters, pd.Index)
    assert not isinstance(clusters, pd.MultiIndex)
    k = clusters.size
    assert k >= 1
    assert clusters.tolist() == list_clusts(k)
    for df in dfs:
        assert isinstance(df, pd.DataFrame)
        assert df.columns.equals(clusters)
    # Calculate matrices of the cost of joining each pair of clusters
    # from each pair of regions.
    cost_matrices = dict()
    for (reg1, df1), (reg2, df2) in combinations(region_params.items(), 2):
        # Even if the regions share no positions, the intersection of
        # their indexes will still include the proportion (0, "p").
        overlap = df1.index.intersection(df2.index)
        assert overlap.size > 0
        # Collect the cost of joining each cluster from region 1 with
        # each cluster from region 2.
        cost_matrix = pd.DataFrame(np.nan, clusters, clusters)
        for cluster1, cluster2 in product(cost_matrix.index,
                                          cost_matrix.columns):
            # Use log odds differences as the costs.
            cost = np.sum(np.abs(
                _calc_diff_log_odds(df1.loc[overlap, cluster1],
                                    df2.loc[overlap, cluster2])
            ))
            cost_matrix.loc[cluster1, cluster2] = cost
        assert not np.any(np.isnan(cost_matrix))
        cost_matrices[reg1, reg2] = cost_matrix
    # Build a hypergraph where every hyperedge connects n nodes, one
    # from each region.
    region_names = list(region_params)
    nodes = {node: i for i, node in enumerate(product(region_names,
                                                      clusters))}
    hyperedges = [tuple(zip(region_names, cluster_nums, strict=True))
                  for cluster_nums in product(clusters, repeat=n)]
    hyperedge_costs = np.array(
        [sum(cost_matrices[reg1, reg2].loc[clust1, clust2]
             for ((reg1, clust1), (reg2, clust2))
             in combinations(hyperedge, 2))
         for hyperedge in hyperedges]
    )
    # Build a sparse boolean matrix where rows are nodes and columns are
    # edges, with a 1 if the edge contains the node and 0 otherwise.
    matrix_rows = list()
    matrix_cols = list()
    for col, hyperedge in enumerate(hyperedges):
        for node in hyperedge:
            row = nodes[node]
            matrix_rows.append(row)
            matrix_cols.append(col)
    num_indices = n * len(hyperedges)
    assert len(matrix_rows) == num_indices
    assert len(matrix_cols) == num_indices
    incidence_matrix = csr_matrix((np.ones(num_indices, dtype=int),
                                   (np.array(matrix_rows),
                                    np.array(matrix_cols))),
                                  shape=(len(nodes), len(hyperedges)))
    # Require every node to appear in exactly one hyperedge.
    node_bounds = np.ones(len(nodes), dtype=int)
    constraints = LinearConstraint(incidence_matrix, node_bounds, node_bounds)
    # Require every possible edge to occur zero or one time.
    integrality = np.ones(len(hyperedges), dtype=bool)
    edge_bounds = Bounds(0, 1)
    # Find the edges that give the smallest cost.
    result = milp(hyperedge_costs,
                  integrality=integrality,
                  bounds=edge_bounds,
                  constraints=constraints)
    if result.status != 0 or not result.success:
        raise RuntimeError(
            f"Failed to determine optimal way to join regions {region_names} "
            f"with {k} cluster(s)"
        )
    # Return a list of the selected hyperedges.
    selected_hyperedges = [hyperedge for hyperedge, is_selected
                           in zip(hyperedges, result.x, strict=True)
                           if is_selected]
    assert len(selected_hyperedges) == k
    return selected_hyperedges


class JoinClusterMutsDataset(ClusterDataset,
                             JoinMutsDataset,
                             MergedUnbiasDataset):

    @classmethod
    def get_report_type(cls):
        return JoinClusterReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_cluster_dataset

    @classmethod
    def get_batch_type(cls):
        return ClusterMutsBatch

    @classmethod
    def name_batch_attrs(cls):
        return [BATCH_NUM, READ_NUMS, SEG_END5S, SEG_END3S, MUTS, RESPS]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # This field must be loaded from the report within __init__()
        # because LoadFunctions find the type of dataset for a report
        # file by checking whether __init__() raises an error, so all
        # fields from the report must be loaded here.
        self._report_joined_clusts = self.report.get_field(JoinedClustersF)

    @cached_property
    def ks(self):
        return self._get_common_attr("ks")

    @cached_property
    def best_k(self):
        return self._get_common_attr("best_k")

    @cached_property
    def clusts(self):
        """ Index of k and cluster numbers. """
        return list_ks_clusts(self.ks)

    @cached_property
    def joined_clusts(self):
        assert self.region_names
        if self._report_joined_clusts:
            if sorted(self._report_joined_clusts) != sorted(self.region_names):
                raise ValueError(
                    f"{self} expected clusters for {self.region_names}, "
                    f"but got {self._report_joined_clusts}"
                )
            return self._report_joined_clusts
        regions_params = {dataset.region.name: _get_clust_params(dataset)
                          for dataset in self.datasets}
        # Determine the best way to join each number of clusters.
        joined_clusts = {reg: {k: dict() for k in self.ks}
                         for reg in regions_params}
        for k in self.ks:
            hyperedges = _join_regions_k(
                {reg: params.loc[:, k]
                 for reg, params in regions_params.items()}
            )
            for hyperedge in hyperedges:
                assert len(hyperedge) >= 1
                _, joined_clust = hyperedge[0]
                for reg, clust in hyperedge:
                    joined_clusts[reg][k][joined_clust] = clust
        return joined_clusts

    def _reg_cols(self, reg: str):
        """ Get the columns for a region's responsibilities. """
        clusts = self.joined_clusts[reg]
        return pd.MultiIndex.from_tuples(
            [(k, clusts[k][clust]) for k, clust in self.clusts],
            names=ClustHeader.level_names()
        )

    def _reg_resps(self, reg: str, resps: pd.DataFrame):
        """ Get the cluster responsibilities for a region. """
        # Reorder the columns.
        reordered = resps.loc[:, self._reg_cols(reg)]
        # Rename the columns by increasing k and cluster.
        reordered.columns = pd.MultiIndex.from_tuples(
            self.clusts,
            names=ClustHeader.level_names()
        )
        return reordered

    def _get_batch_attrs(self, batch: MutsBatch, reg: str):
        attrs = super()._get_batch_attrs(batch, reg)
        # Adjust the cluster labels based on the region.
        attrs[RESPS] = self._reg_resps(reg, attrs[RESPS])
        return attrs

    def _join_attrs(self, attrs: dict[str, Any], add_attrs: dict[str, Any]):
        super()._join_attrs(attrs, add_attrs)
        # Join the cluster memberships taking the mean over all regions.
        # Because there can be more than two regions, accumulate the sum
        # here rather than taking the mean; then in _finalize_attrs(),
        # normalize to ensure the sum is 1 for each number of clusters.
        attrs[RESPS] = attrs[RESPS].add(add_attrs[RESPS], fill_value=0.)

    def _finalize_attrs(self, attrs: dict[str, Any]):
        # Ensure that cluster memberships for each read sum to 1.
        attrs[RESPS] /= attrs[RESPS].T.groupby(level=NUM_CLUSTS_NAME).sum().T
        # Fill any missing values with 0 and sort the read numbers.
        attrs[RESPS] = attrs[RESPS].fillna(0.).sort_index()
        # Delete read_nums (which is the index of resps).
        attrs.pop(READ_NUMS)


load_cluster_dataset = LoadFunction(ClusterMutsDataset, JoinClusterMutsDataset)

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
