from abc import ABC
from functools import cached_property
from itertools import combinations, product
from typing import Any

import numpy as np
import pandas as pd

from .batch import ClusterMutsBatch
from .io import ClusterFile, ClusterBatchIO
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
                           RelClustHeader,
                           list_clusts,
                           list_ks_clusts,
                           validate_ks,
                           make_header,
                           parse_header)
from ..core.join import (BATCH_NUM,
                         READ_NUMS,
                         SEG_END5S,
                         SEG_END3S,
                         MUTS,
                         RESPS,
                         JoinMutsDataset)
from ..core.logs import logger
from ..core.mu import calc_sum_arcsine_distance
from ..core.report import JoinedClustersF, KsWrittenF, BestKF
from ..core.seq import POS_NAME, BASE_NAME
from ..core.table import (MUTAT_REL,
                          TableLoader,
                          PositionTableLoader,
                          BatchTabulator,
                          CountTabulator,
                          AbundanceTable,
                          RelTypeTable,
                          PositionTableWriter,
                          AbundanceTableWriter)
from ..mask.batch import MaskMutsBatch
from ..mask.dataset import load_mask_dataset
from ..mask.table import (PartialTable,
                          PartialPositionTable,
                          PartialTabulator,
                          PartialDatasetTabulator)


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
        return self.dataset1.pattern

    @pattern.setter
    def pattern(self, pattern):
        self.dataset1.pattern = pattern

    @property
    def region(self):
        return self.dataset1.region

    @property
    def min_mut_gap(self):
        return getattr(self.dataset1, "min_mut_gap")

    @min_mut_gap.setter
    def min_mut_gap(self, min_mut_gap):
        self.dataset1.min_mut_gap = min_mut_gap

    @property
    def quick_unbias(self):
        return getattr(self.dataset1, "quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return getattr(self.dataset1, "quick_unbias_thresh")

    @property
    def ks(self):
        return getattr(self.dataset2, "ks")

    @property
    def best_k(self):
        return getattr(self.dataset2, "best_k")

    def _integrate(self, batch1: MaskMutsBatch, batch2: ClusterBatchIO):
        return ClusterMutsBatch(batch=batch1.batch,
                                region=batch1.region,
                                seg_end5s=batch1.seg_end5s,
                                seg_end3s=batch1.seg_end3s,
                                muts=batch1.muts,
                                resps=batch2.resps,
                                sanitize=False)


def get_clust_params(dataset: ClusterMutsDataset, num_cpus: int = 1):
    """ Get the mutation rates and proportion for each cluster. If table
    files already exist, then use them to get the parameters; otherwise,
    calculate the parameters from the dataset. """
    logger.routine(f"Began obtaining cluster parameters from {dataset}")
    # Try to load the tables from files.
    path_fields = {path.TOP: dataset.top,
                   path.SAMPLE: dataset.sample,
                   path.BRANCHES: dataset.branches,
                   path.REF: dataset.ref,
                   path.REG: dataset.region.name}
    pos_table_file = ClusterPositionTableLoader.build_path(path_fields)
    if pos_table_file.is_file():
        pos_table = ClusterPositionTableLoader(pos_table_file)
        logger.detail(f"Position table {pos_table_file} exists")
    else:
        pos_table = None
        logger.detail(f"Position table {pos_table_file} does not exist")
    abundance_table_file = ClusterAbundanceTableLoader.build_path(path_fields)
    if abundance_table_file.is_file():
        abundance_table = ClusterAbundanceTableLoader(abundance_table_file)
        logger.detail(f"Abundance table {abundance_table_file} exists")
    else:
        abundance_table = None
        logger.detail(f"Abundance table {abundance_table_file} does not exist")
    # If either table file does not exist, then calculate the tables.
    if pos_table is None or abundance_table is None:
        logger.detail(
            "Tabulating is needed because at least one table does not exist"
        )
        tabulator = ClusterBatchTabulator(
            top=dataset.top,
            sample=dataset.sample,
            branches=dataset.branches,
            region=dataset.region,
            refseq=dataset.refseq,
            pattern=dataset.pattern,
            min_mut_gap=dataset.min_mut_gap,
            quick_unbias=dataset.quick_unbias,
            quick_unbias_thresh=dataset.quick_unbias_thresh,
            ks=dataset.ks,
            num_batches=dataset.num_batches,
            get_batch_count_all=dataset.get_batch_count_all,
            count_ends=True,
            count_pos=(pos_table is None),
            count_read=False,
            num_cpus=num_cpus
        )
        if pos_table is None:
            pos_table = ClusterPositionTableWriter(tabulator)
        if abundance_table is None:
            abundance_table = ClusterAbundanceTableWriter(tabulator)
    else:
        logger.detail("Tabulating is not needed because all tables exist")
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
    logger.routine(f"Ended obtaining cluster parameters from {dataset}")
    return params


def _join_regions_k(region_params: dict[str, pd.DataFrame]):
    """ Determine the optimal way to join regions . """
    logger.routine("Began determining the optimal way to join clusters")
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
    logger.detail(f"There are {n} regions and {k} clusters")
    # Calculate matrices of the cost of joining each pair of clusters
    # from each pair of regions.
    cost_matrices = dict()
    for (reg1, df1), (reg2, df2) in combinations(region_params.items(), 2):
        # Even if the regions share no positions, their indexes will
        # both include the proportion (0, "p").
        overlap = df1.index.intersection(df2.index)
        assert overlap.size > 0
        logger.detail(f"Regions {repr(reg1)} and {repr(reg2)} "
                      f"share {overlap.size} parameter(s)")
        # Collect the cost of joining each cluster from region 1 with
        # each cluster from region 2.
        cost_matrix = pd.DataFrame(np.nan, clusters, clusters)
        for cluster1, cluster2 in product(clusters, repeat=2):
            # Use total arcsine distances as the costs.
            cost = calc_sum_arcsine_distance(df1.loc[overlap, cluster1],
                                             df2.loc[overlap, cluster2])
            cost_matrix.at[cluster1, cluster2] = cost
        logger.detail(f"Regions {repr(reg1)} and {repr(reg2)} "
                      f"have a cost matrix of\n{cost_matrix}")
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
        [sum(cost_matrices[reg1, reg2].at[clust1, clust2]
             for ((reg1, clust1), (reg2, clust2))
             in combinations(hyperedge, 2))
         for hyperedge in hyperedges]
    )
    logger.detail(f"Built a hypergraph with {len(nodes)} node(s) "
                  f"and {len(hyperedges)} hyperedge(s)")
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
    logger.detail("Created mixed-integer linear program: "
                  "min_x(cx), subject to Ax = 1, x ∈ {0, 1}; "
                  f"c and x are length {len(hyperedges)}, "
                  f"and A has dimensions {incidence_matrix.shape}")
    # Find the edges that give the smallest cost.
    logger.detail("Began solving mixed-integer linear program "
                  "(this could take a while)")
    result = milp(hyperedge_costs,
                  integrality=integrality,
                  bounds=edge_bounds,
                  constraints=constraints)
    if result.status != 0 or not result.success:
        raise RuntimeError(
            f"Failed to determine optimal way to join regions {region_names} "
            f"with {k} clusters"
        )
    logger.detail("Ended solving mixed-integer linear program: "
                  f"minimum total cost is {result.fun}")
    # Return a list of the selected hyperedges.
    selected_hyperedges = [hyperedge for hyperedge, is_selected
                           in zip(hyperedges, result.x, strict=True)
                           if is_selected]
    assert len(selected_hyperedges) == k
    logger.detail(f"Selected {k} hyperedges:\n"
                  + "\n".join(map(str, selected_hyperedges)))
    logger.routine("Ended determining the optimal way to join clusters")
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
        report_joined_clusts = self.report.get_field(JoinedClustersF)
        assert self.region_names
        if report_joined_clusts:
            if sorted(report_joined_clusts) != sorted(self.region_names):
                raise ValueError(
                    f"{self} expected clusters for {self.region_names}, "
                    f"but got {report_joined_clusts}"
                )
            self.joined_clusts = report_joined_clusts
        else:
            regions_params = {dataset.region.name: get_clust_params(dataset)
                              for dataset in self.datasets}
            # Determine the best way to join each number of clusters.
            optimal_joined_clusts = {reg: {k: dict() for k in self.ks}
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
                        optimal_joined_clusts[reg][k][joined_clust] = clust
            self.joined_clusts = optimal_joined_clusts

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

    def _reg_cols(self, reg: str):
        """ Get the columns for a region's responsibilities. """
        clusts = self.joined_clusts[reg]
        return pd.MultiIndex.from_tuples(
            [(k, clusts[k][clust]) for k, clust in self.clusts],
            names=ClustHeader.get_level_names()
        )

    def _reg_resps(self, reg: str, resps: pd.DataFrame):
        """ Get the cluster responsibilities for a region. """
        # Reorder the columns.
        reordered = resps.loc[:, self._reg_cols(reg)]
        # Rename the columns by increasing k and cluster.
        reordered.columns = pd.MultiIndex.from_tuples(
            self.clusts,
            names=ClustHeader.get_level_names()
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
        # divide by the sum for each number of clusters.
        attrs[RESPS] = attrs[RESPS].add(add_attrs[RESPS], fill_value=0.)

    def _finalize_attrs(self, attrs: dict[str, Any]):
        # Ensure that cluster memberships for each read sum to 1.
        attrs[RESPS] /= attrs[RESPS].T.groupby(level=NUM_CLUSTS_NAME).sum().T
        # Fill any missing values with 0 and sort the read numbers.
        attrs[RESPS] = attrs[RESPS].fillna(0.).sort_index()
        # Delete read_nums (which is the index of resps).
        attrs.pop(READ_NUMS)


load_cluster_dataset = LoadFunction(ClusterMutsDataset, JoinClusterMutsDataset)


class ClusterTable(RelTypeTable, ClusterFile, ABC):

    @classmethod
    def get_load_function(cls):
        return load_cluster_dataset

    @classmethod
    def get_header_type(cls):
        return RelClustHeader


class ClusterPositionTable(ClusterTable, PartialPositionTable, ABC):
    pass


class ClusterAbundanceTable(AbundanceTable, PartialTable, ClusterFile, ABC):

    @classmethod
    def get_load_function(cls):
        return load_cluster_dataset

    @classmethod
    def get_header_type(cls):
        return ClustHeader

    @classmethod
    def get_index_depth(cls):
        return cls.get_header_depth()

    def _get_header(self):
        return parse_header(self.data.index)

    @cached_property
    def proportions(self):
        return self.data / self.data.groupby(level=NUM_CLUSTS_NAME).sum()


class ClusterPositionTableWriter(PositionTableWriter, ClusterPositionTable):
    pass


class ClusterAbundanceTableWriter(AbundanceTableWriter, ClusterAbundanceTable):
    pass


class ClusterPositionTableLoader(PositionTableLoader, ClusterPositionTable):
    """ Load cluster data indexed by position. """


class ClusterAbundanceTableLoader(TableLoader, ClusterAbundanceTable):
    """ Load cluster data indexed by cluster. """

    @cached_property
    def data(self) -> pd.Series:
        data = pd.read_csv(self.path,
                           index_col=self.get_index_cols()).squeeze(axis=1)
        if not isinstance(data, pd.Series):
            raise ValueError(f"{self} must have one column, but got\n{data}")
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.index)
        # The index must be replaced with the header index for the
        # type casting to take effect.
        data.index = header.index
        return data


class ClusterTabulator(PartialTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [ClusterPositionTableWriter, ClusterAbundanceTableWriter]

    def __init__(self, *, ks: list[int], **kwargs):
        super().__init__(**kwargs)
        if ks is None:
            raise ValueError(
                f"{type(self).__name__} requires clusters, but got ks={ks}"
            )
        self.ks = ks

    @cached_property
    def clust_header(self):
        """ Header of the per-cluster data. """
        return make_header(ks=self.ks)

    @cached_property
    def data_per_clust(self):
        """ Number of reads in each cluster. """
        n_rels, n_clust = self._adjusted
        n_clust.name = "Number of Reads"
        return n_clust


class ClusterBatchTabulator(BatchTabulator, ClusterTabulator):
    pass


class ClusterCountTabulator(CountTabulator, ClusterTabulator):
    pass


class ClusterDatasetTabulator(PartialDatasetTabulator, ClusterTabulator):

    @classmethod
    def init_kws(cls):
        return super().init_kws() + ["ks"]
