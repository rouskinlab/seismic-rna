from abc import ABC
from functools import cache, cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .dataset import load_mask_dataset
from .io import MaskFile
from ..core import path
from ..core.batch import END5_COORD, END3_COORD
from ..core.header import NUM_CLUSTS_NAME, format_clust_name, validate_ks
from ..core.logs import logger
from ..core.rel import RelPattern
from ..core.rna import RNAProfile
from ..core.seq import DNA, Region
from ..core.table import (COVER_REL,
                          MATCH_REL,
                          MUTAT_REL,
                          INFOR_REL,
                          SUBMUTS,
                          Tabulator,
                          BatchTabulator,
                          CountTabulator,
                          DatasetTabulator,
                          Table,
                          PositionTable,
                          ReadTable,
                          PositionTableLoader,
                          ReadTableLoader,
                          PositionTableWriter,
                          ReadTableWriter)
from ..core.unbias import (calc_p_ends_observed,
                           calc_p_noclose_given_clust,
                           calc_p_noclose_given_ends_auto,
                           calc_params)
from ..relate.table import (AverageTable,
                            AverageTabulator)


class PartialTable(Table, path.HasRegFilePath, ABC):
    """ Table of filtered reads over a region of the sequence. """


class PartialPositionTable(PartialTable, PositionTable, ABC):

    def _iter_profiles(self, *,
                       regions: Iterable[Region] | None,
                       quantile: float,
                       rel: str,
                       k: int | None,
                       clust: int | None):
        """ Yield RNA mutational profiles from a table. """
        if regions is not None:
            regions = list(regions)
        else:
            regions = [self.region]
        for hk, hc in self.header.clusts:
            if (k is None or k == hk) and (clust is None or clust == hc):
                data_name = path.fill_whitespace(format_clust_name(hk, hc),
                                                 fill="-")
                for region in regions:
                    yield RNAProfile(region=region,
                                     sample=self.sample,
                                     branches=self.branches,
                                     data_reg=self.reg,
                                     data_name=data_name,
                                     data=self.fetch_ratio(quantile=quantile,
                                                           rel=rel,
                                                           k=hk,
                                                           clust=hc,
                                                           squeeze=True))


class PartialReadTable(PartialTable, ReadTable, ABC):
    pass


class MaskTable(AverageTable, MaskFile, ABC):

    @classmethod
    def get_load_function(cls):
        return load_mask_dataset


class MaskPositionTable(MaskTable, PartialPositionTable, ABC):
    pass


class MaskReadTable(MaskTable, PartialReadTable, ABC):
    pass


class MaskPositionTableWriter(PositionTableWriter, MaskPositionTable):
    pass


class MaskReadTableWriter(ReadTableWriter, MaskReadTable):
    pass


class MaskPositionTableLoader(PositionTableLoader, MaskPositionTable):
    pass


class MaskReadTableLoader(ReadTableLoader, MaskReadTable):
    pass


class PartialTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return np.nan

    def __init__(self, *,
                 refseq: DNA,
                 region: Region,
                 pattern: RelPattern,
                 min_mut_gap: int,
                 quick_unbias: bool,
                 quick_unbias_thresh: float,
                 count_ends: bool = True,
                 **kwargs):
        # Partial tabulators must count 5'/3' ends or else calculating
        # self.p_ends_given_clust_noclose will fail.
        if not count_ends:
            logger.warning(f"count_ends must be True for {type(self.__name__)}")
        super().__init__(region=region, count_ends=True, **kwargs)
        self.refseq = refseq
        self.pattern = pattern
        self.min_mut_gap = min_mut_gap
        self.quick_unbias = quick_unbias
        self.quick_unbias_thresh = quick_unbias_thresh

    @cached_property
    def p_ends_given_clust_noclose(self):
        """ Probability of each end coordinate. """
        # Ensure end_counts has 2 dimensions.
        if self.end_counts.ndim == 1:
            end_counts = self.end_counts.values[:, np.newaxis]
        else:
            end_counts = self.end_counts.values
        end5s = (self.end_counts.index.get_level_values(END5_COORD).values
                 - self.region.end5)
        end3s = (self.end_counts.index.get_level_values(END3_COORD).values
                 - self.region.end5)
        return calc_p_ends_observed(self.region.length,
                                    end5s,
                                    end3s,
                                    end_counts)

    @cached_property
    def _adjusted(self):
        table_per_pos = super().data_per_pos
        if self.min_mut_gap > 0:
            if self.region.length > np.sqrt(1_000_000_000):
                logger.warning("Using bias correction on a region with "
                               f"{self.region.length} positions requires "
                               ">1 GB of memory. If this is impractical, you "
                               "can (at the cost of lower accuracy) disable "
                               "bias correction using --min-mut-gap 0.")
            try:
                return adjust_counts(table_per_pos,
                                     self.p_ends_given_clust_noclose,
                                     self.num_reads,
                                     self.region,
                                     self.min_mut_gap,
                                     self.quick_unbias,
                                     self.quick_unbias_thresh)
            except Exception as error:
                logger.warning(error)
        return table_per_pos, self.num_reads

    @cached_property
    def data_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        n_rels, n_clust = self._adjusted
        return n_rels


class PartialDatasetTabulator(DatasetTabulator, PartialTabulator, ABC):

    @classmethod
    @cache
    def init_kws(cls):
        return super().init_kws() + cls._list_args(PartialTabulator.__init__)


class MaskTabulator(PartialTabulator, AverageTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [MaskPositionTableWriter, MaskReadTableWriter]


class MaskBatchTabulator(BatchTabulator, MaskTabulator):
    pass


class MaskCountTabulator(CountTabulator, MaskTabulator):
    pass


class MaskDatasetTabulator(PartialDatasetTabulator, MaskTabulator):
    pass


def _insert_masked(p_mut: pd.Series | pd.DataFrame,
                   region: Region):
    """ 2D array where masked positions are filled with 0. """
    # Fill masked positions with 0.
    p_mut = p_mut.reindex(index=region.range, fill_value=0.)
    # Convert to a 2D NumPy array.
    return p_mut.values.reshape((region.length, -1))


def adjust_counts(table_per_pos: pd.DataFrame,
                  p_ends_given_clust_noclose: np.ndarray,
                  n_reads_clust: pd.Series | int,
                  region: Region,
                  min_mut_gap: int,
                  quick_unbias: bool,
                  quick_unbias_thresh: float):
    """ Adjust the given table of masked/clustered counts per position
    to correct for observer bias. """
    if not isinstance(table_per_pos, pd.DataFrame):
        raise TypeError(table_per_pos)
    action = (f"unbiasing counts (min_mut_gap={min_mut_gap}) "
              f"of table with {table_per_pos.index.size} positions")
    if isinstance(n_reads_clust, pd.Series):
        k = n_reads_clust.size
        action += f", {k} clusters,"
        n_rels, rem = divmod(table_per_pos.columns.size, k)
        if rem:
            raise ValueError(
                f"Number of columns in table ({table_per_pos.columns.size}) "
                f"is not a multiple of number of clusters ({k})"
            )
    else:
        n_rels = table_per_pos.columns.size
    action += f" and {n_rels} relationships"
    logger.routine(f"Began {action}")
    # Determine which positions are unmasked.
    unmask = region.unmasked_bool
    # Calculate the fraction of mutations at each position among reads
    # with no two mutations too close.
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # resulting NaN values are zeroed by nan_to_num.
        p_mut_given_noclose = np.nan_to_num(_insert_masked(
            table_per_pos[MUTAT_REL] / table_per_pos[INFOR_REL],
            region
        ))
    if isinstance(n_reads_clust, int):
        # There is only one cluster, so the probability that each read
        # belongs to that cluster is 1.
        p_clust_given_noclose = np.array([1.])
        # Calculate the parameters.
        p_mut, p_ends, p_clust = calc_params(
            p_mut_given_noclose,
            p_ends_given_clust_noclose,
            p_clust_given_noclose,
            min_mut_gap,
            quick_unbias=quick_unbias,
            quick_unbias_thresh=quick_unbias_thresh
        )
        # Compute the probability that reads would have no two mutations
        # too close.
        p_noclose_given_clust = calc_p_noclose_given_clust(
            p_ends,
            calc_p_noclose_given_ends_auto(p_mut, min_mut_gap)
        )
        # Drop the cluster dimension from the parameters.
        if p_mut.shape != (region.length, 1):
            raise ValueError(f"p_mut must have shape {(region.length, 1)}, "
                             f"but got {p_mut.shape}")
        p_mut = p_mut.reshape((-1,))
        if p_clust.shape != (1,):
            raise ValueError(
                f"p_clust must have shape {1,}, but got {p_clust.shape}"
            )
        p_clust = float(p_clust[0])
        if not np.isclose(p_clust, 1.):
            raise ValueError(f"p_clust must equal 1., but got {p_clust}")
        if p_noclose_given_clust.shape != (1,):
            raise ValueError(f"p_noclose_given_clust must have shape {1,}, "
                             f"but got {p_noclose_given_clust.shape}")
        p_noclose = p_noclose_given_clust = float(p_noclose_given_clust[0])
        # Compute the number of reads.
        n_clust = n_reads_clust / p_noclose
    elif isinstance(n_reads_clust, pd.Series):
        # Calculate the number of reads with no two mutations too close
        # for each k.
        n_reads_noclose_ks = n_reads_clust.groupby(level=NUM_CLUSTS_NAME).sum()
        # Determine the numbers of clusters.
        ks = validate_ks(n_reads_noclose_ks.index.values)
        # Calculate the parameters for each k separately.
        p_mut = np.empty_like(p_mut_given_noclose)
        p_clust = np.empty_like(n_reads_clust.values)
        p_noclose_given_clust = np.empty_like(n_reads_clust.values)
        n_clust = pd.Series(index=n_reads_clust.index)
        for k in ks:
            ki = n_reads_clust.index.get_level_values(NUM_CLUSTS_NAME) == k
            # Calculate the fraction of reads with no two mutations too
            # close in each cluster.
            n_reads_noclose = float(n_reads_noclose_ks.at[k])
            p_clust_given_noclose = (n_reads_clust.loc[k].values
                                     / n_reads_noclose)
            # Calculate the parameters for each cluster.
            p_mut[:, ki], p_ends, p_clust[ki] = calc_params(
                p_mut_given_noclose[:, ki],
                p_ends_given_clust_noclose[:, :, ki],
                p_clust_given_noclose,
                min_mut_gap
            )
            # Compute the probability that reads from each cluster would
            # have no two mutations too close.
            p_noclose_given_clust[ki] = calc_p_noclose_given_clust(
                p_ends,
                calc_p_noclose_given_ends_auto(p_mut[:, ki], min_mut_gap)
            )
            # Compute the probability that reads from any cluster would
            # have no two mutations too close.
            p_noclose = float(np.vdot(p_noclose_given_clust[ki], p_clust[ki]))
            # Compute the number of reads in each cluster.
            n_clust.at[k] = (n_reads_noclose / p_noclose) * p_clust[ki]
    else:
        raise TypeError("n_reads_clust must be an int or Series, "
                        f"but got {type(n_reads_clust).__name__}")
    # Remove masked positions from the mutation rates.
    p_mut = p_mut[unmask]
    # Create the table of adjusted counts.
    # Initialize an empty DataFrame of the adjusted counts with the same
    # index and columns as the observed counts.
    n_rels = pd.DataFrame(np.nan, table_per_pos.index, table_per_pos.columns)
    # Assume that the observance bias affects the counts of covered and
    # informative bases equally, so p_noclose_given_ends is defined as:
    # p_noclose_given_ends := ninfo_obs / n_info
    # from which we can estimate the informative bases after adjustment:
    # n_info = ninfo_obs / p_noclose_given_ends
    n_cov = table_per_pos.loc[unmask, COVER_REL].values / p_noclose_given_clust
    n_info = table_per_pos.loc[unmask, INFOR_REL].values / p_noclose_given_clust
    n_rels.loc[unmask, COVER_REL] = n_cov
    n_rels.loc[unmask, INFOR_REL] = n_info
    # From the definition of the adjusted fraction of mutations:
    # p_mut := n_mut / n_info
    # we can also estimate the mutated bases after adjustment:
    n_mut = p_mut * n_info
    n_rels.loc[unmask, MUTAT_REL] = n_mut
    # From the definition of informative bases:
    # n_info := n_ref + n_mut
    # we can estimate the matched bases after adjustment:
    n_ref = n_info - n_mut
    n_rels.loc[unmask, MATCH_REL] = n_ref
    # Compute the factor by which n_mut was scaled for each position.
    with np.errstate(divide="ignore", invalid="ignore"):
        # Division by 0 is possible if no mutations were observed at a
        # given position, resulting in a NaN value at that position.
        scale = n_mut / table_per_pos.loc[unmask, MUTAT_REL].values
    # Replace NaN values with 1 so that missing values do not propagate
    # during multiplication.
    scale = np.nan_to_num(scale, nan=1.)
    # Scale every subtype of mutation by this factor.
    for mut in SUBMUTS:
        n_rels.loc[unmask, mut] = scale * table_per_pos.loc[unmask, mut].values
    logger.routine(f"Ended {action}")
    return n_rels, n_clust
