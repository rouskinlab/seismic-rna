from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .base import (COVER_REL,
                   DELET_REL,
                   INSRT_REL,
                   MATCH_REL,
                   MUTAT_REL,
                   SUBST_REL,
                   SUB_A_REL,
                   SUB_C_REL,
                   SUB_G_REL,
                   SUB_T_REL,
                   INFOR_REL,
                   R_ADJ_TITLE,
                   R_OBS_TITLE,
                   TABLE_RELS)
from ..cluster.data import ClusterMutsDataset
from ..core.batch import END5_COORD, END3_COORD, accum_fits
from ..core.header import ORDER_NAME, Header, make_header
from ..core.mu import (calc_p_ends_given_noclose,
                       calc_p_noclose,
                       calc_p_noclose_given_ends_frame,
                       calc_params_frame)
from ..core.rel import RelPattern, HalfRelPattern
from ..core.seq import Section
from ..mask.data import MaskMutsDataset
from ..pool.data import PoolDataset
from ..relate.data import RelateDataset

logger = getLogger(__name__)

# These relationships are of all subtypes of mutations.
SUBMUTS = [SUBST_REL,
           SUB_A_REL,
           SUB_C_REL,
           SUB_G_REL,
           SUB_T_REL,
           DELET_REL,
           INSRT_REL]


# Tabulator Classes ####################################################

class Tabulator(ABC):
    """ Base class for tabulating data for multiple tables from a report
    loader. """

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    @classmethod
    def _format_table(cls, counts: pd.DataFrame, header: Header):
        """ Format the count data with the proper header. """
        # Initialize an empty table.
        table = pd.DataFrame(cls.get_null_value(), counts.index, header.index)
        # Count reads with each relationship at each position.
        for rel, rel_counts in counts.items():
            table.loc[:, rel] = rel_counts
        # Add a column for informative relationships.
        table.loc[:, INFOR_REL] = (table.loc[:, MATCH_REL].values
                                   +
                                   table.loc[:, MUTAT_REL].values)
        return table

    def __init__(self, dataset: (RelateDataset
                                 | PoolDataset
                                 | MaskMutsDataset
                                 | ClusterMutsDataset)):
        self.dataset = dataset

    @property
    def top(self):
        return self.dataset.top

    @property
    def sample(self):
        return self.dataset.sample

    @property
    def ref(self):
        return self.dataset.ref

    @property
    def refseq(self):
        return self.dataset.refseq

    @property
    def section(self):
        return self.dataset.section

    @property
    @abstractmethod
    def max_order(self) -> int:
        """ Number of clusters, or 0 if not clustered. """

    @cached_property
    def pos_header(self):
        """ Header of the per-position data. """
        return make_header(rels=TABLE_RELS, max_order=self.max_order)

    @cached_property
    def read_header(self):
        """ Header of the per-read data. """
        return make_header(rels=TABLE_RELS, max_order=self.max_order)

    @cached_property
    def _counts(self):
        return accum_fits(self.dataset.iter_batches(),
                          self.refseq,
                          self.section.unmasked_int,
                          all_patterns(self.dataset.pattern),
                          max_order=self.max_order)

    @property
    def _num_reads(self):
        """ Raw number of reads. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return num_reads

    @property
    def _counts_per_pos(self):
        """ Raw counts per position. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return per_pos

    @property
    def _counts_per_read(self):
        """ Raw counts per read. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return per_read

    @property
    def _end_counts(self):
        """ Raw counts for each pair of end coordinates. """
        num_reads, per_pos, per_read, end_counts = self._counts
        return end_counts

    @cached_property
    def p_ends_given_noclose(self):
        """ Probability of each end coordinate. """
        # Ensure end_counts has 2 dimensions.
        if self._end_counts.ndim == 1:
            end_counts = self._end_counts.values[:, np.newaxis]
        else:
            end_counts = self._end_counts.values
        return calc_p_ends_given_noclose(
            self.section.length,
            (self._end_counts.index.get_level_values(END5_COORD).values
             - self.section.end5),
            (self._end_counts.index.get_level_values(END3_COORD).values
             - self.section.end5),
            end_counts,
        )

    @cached_property
    def p_clust_given_noclose(self):
        """ Probability of each cluster. """
        if isinstance(self._num_reads, pd.Series):
            return (self._num_reads
                    / self._num_reads.groupby(level=ORDER_NAME).sum())
        if isinstance(self._num_reads, int):
            return 1.
        raise TypeError("Number of reads must be Series or int, "
                        f"but got {type(self._num_reads).__name__}")

    @cached_property
    def table_per_pos(self):
        return self._format_table(self._counts_per_pos,
                                  self.pos_header).reindex(self.section.range)

    @cached_property
    def table_per_read(self):
        return self._format_table(self._counts_per_read,
                                  self.read_header)


class FullTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return 0


class PartialTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @cached_property
    def _adjusted(self):
        return adjust_counts(super().table_per_pos,
                             self.p_ends_given_noclose,
                             self.p_clust_given_noclose,
                             self.section,
                             self.dataset.min_mut_gap)

    @cached_property
    def table_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        n_rels, n_ends, n_clust = self._adjusted
        return n_rels


class AvgTabulator(Tabulator, ABC):

    @property
    def max_order(self) -> int:
        return 0


class RelateTabulator(FullTabulator, AvgTabulator):
    pass


class MaskTabulator(PartialTabulator, AvgTabulator):
    pass


class ClustTabulator(PartialTabulator, ABC):

    @property
    def max_order(self):
        return self.dataset.max_order

    @cached_property
    def num_reads_obs(self):
        """ Observed number of reads in each cluster. """
        return self._num_reads

    @cached_property
    def num_reads_adj(self):
        """ Adjusted number of reads in each cluster. """
        n_rels, n_ends, n_clust = self._adjusted
        return n_clust

    @cached_property
    def clust_header(self):
        """ Header of the per-cluster data. """
        return make_header(max_order=self.max_order)

    @cached_property
    def table_per_clust(self):
        """ Observed and adjusted numbers of reads in each cluster. """
        return pd.DataFrame.from_dict({R_OBS_TITLE: self.num_reads_obs,
                                       R_ADJ_TITLE: self.num_reads_adj}).T


# Helper functions #####################################################

def _iter_mut_patterns():
    """ Yield a HalfRelPattern for each type of mutation. """
    yield SUBST_REL, HalfRelPattern.from_counts(count_sub=True)
    yield SUB_A_REL, HalfRelPattern("ca", "ga", "ta")
    yield SUB_C_REL, HalfRelPattern("ac", "gc", "tc")
    yield SUB_G_REL, HalfRelPattern("ag", "cg", "tg")
    yield SUB_T_REL, HalfRelPattern("at", "ct", "gt")
    yield DELET_REL, HalfRelPattern.from_counts(count_del=True)
    yield INSRT_REL, HalfRelPattern.from_counts(count_ins=True)


def _iter_patterns(mask: RelPattern | None = None):
    """ Yield a RelPattern for every type of relationship. """
    # Count everything except for no coverage.
    yield COVER_REL, RelPattern.allc()
    # Count matches to the reference sequence.
    yield MATCH_REL, RelPattern.muts().intersect(mask, invert=True)
    # Count all types of mutations, relative to reference matches.
    yield MUTAT_REL, RelPattern.muts().intersect(mask)
    # Count each type of mutation, relative to reference matches.
    for mut, pattern in _iter_mut_patterns():
        yield mut, RelPattern(pattern, HalfRelPattern.refs()).intersect(mask)


def all_patterns(mask: RelPattern | None = None):
    """ Every RelPattern, keyed by its name. """
    return dict(_iter_patterns(mask))


def adjust_counts(table_per_pos: pd.DataFrame,
                  p_ends_given_noclose: np.ndarray,
                  p_clust_given_noclose: pd.Series | float,
                  section: Section,
                  min_mut_gap: int):
    """ Adjust the given table of masked/clustered counts per position
    to correct for observer bias.

    Parameters
    ----------
    p_mut_given_span_noclose: DataFrame
        Counts of the bits for each type of relation (column) at each
        position (index) in the section of interest.
    section: Section
        The section of interest.
    min_mut_gap: int
        Minimum number of non-mutated bases permitted between mutations.
    """
    # Compute the observed fraction of mutations at each position.
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # NaN values will be zeroed out subsequently.
        p_mut_given_span_noclose = (table_per_pos[MUTAT_REL]
                                    / table_per_pos[INFOR_REL])
    # Fill any missing (NaN) values with zero.
    p_mut_given_span_noclose.fillna(0., inplace=True)
    # Initialize an empty DataFrame of the adjusted counts with the same
    # index and columns as the observed counts.
    n_rels = pd.DataFrame(np.nan, table_per_pos.index, table_per_pos.columns)
    # Adjust the fraction of mutations to correct the observer bias.
    p_mut, p_ends, p_clust = calc_params_frame(section,
                                               p_mut_given_span_noclose,
                                               p_ends_given_noclose,
                                               p_clust_given_noclose,
                                               min_mut_gap)
    # Compute the probability that a read from each cluster would have
    # no two mutations too close.
    p_noclose = calc_p_noclose(p_ends,
                               calc_p_noclose_given_ends_frame(section,
                                                               p_mut,
                                                               min_mut_gap))
    if isinstance(p_clust_given_noclose, pd.Series):
        p_noclose = pd.Series(p_noclose, index=p_clust_given_noclose.index)
    elif isinstance(p_clust_given_noclose, float):
        if p_noclose.size > 1 or p_noclose.ndim > 1:
            raise ValueError(f"p_noclose must be a scalar, but got {p_noclose}")
        p_noclose = float(p_noclose[0] if p_noclose.ndim else p_noclose)
    else:
        raise TypeError("p_clust_given_noclose must be a Series or float, "
                        f"but got {type(p_clust_given_noclose).__name__}")
    # Assume that the observance bias affects the counts of covered and
    # informative bases equally, so p_noclose_given_ends is defined as:
    # p_noclose_given_ends := ninfo_obs / n_info
    # from which we can estimate the informative bases after adjustment:
    # n_info = ninfo_obs / p_noclose_given_ends
    n_info = table_per_pos[INFOR_REL] / p_noclose
    n_rels.loc[:, INFOR_REL] = n_info.values
    n_rels.loc[:, COVER_REL] = (table_per_pos[COVER_REL] / p_noclose).values
    # From the definition of the adjusted fraction of mutations:
    # p_mut := n_mut / n_info
    # we can also estimate the mutated bases after adjustment:
    # n_mut = p_mut * n_info
    n_mut = p_mut * n_info
    n_rels.loc[:, MUTAT_REL] = n_mut.values
    # From the definition of informative bases:
    # n_info := n_ref + n_mut
    # we can estimate the matched bases after adjustment:
    # n_ref = n_info - n_mut
    n_ref = n_info - n_mut
    n_rels.loc[:, MATCH_REL] = n_ref.values
    # Compute the factor by which n_mut was scaled for each position.
    with np.errstate(divide="ignore"):
        # Division by 0 is possible if no mutations were observed at a
        # given position, resulting in a NaN value at that position.
        scale = n_mut / table_per_pos[MUTAT_REL]
    # Replace NaN values with 1 so that missing values do not propagate
    # during multiplication.
    scale = scale.fillna(1.)
    # Scale every subtype of mutation by this factor.
    for mut in SUBMUTS:
        n_rels.loc[:, mut] = (table_per_pos[mut] * scale).values
    # Calculate the count of each end coordinate and cluster.
    n_ends = (p_ends[:, :, np.newaxis]
              if isinstance(p_noclose, pd.Series)
              else p_ends) / np.asarray(p_noclose)
    n_clust = p_clust / p_noclose
    return n_rels, n_ends, n_clust


def tabulate_loader(dataset: (RelateDataset
                              | PoolDataset
                              | MaskMutsDataset
                              | ClusterMutsDataset)):
    """ Return a new Dataset, choosing the subclass based on the type
    of the argument `dataset`. """
    if isinstance(dataset, (RelateDataset, PoolDataset)):
        return RelateTabulator(dataset)
    if isinstance(dataset, MaskMutsDataset):
        return MaskTabulator(dataset)
    if isinstance(dataset, ClusterMutsDataset):
        return ClustTabulator(dataset)
    raise TypeError(f"Invalid dataset type: {type(dataset).__name__}")

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
