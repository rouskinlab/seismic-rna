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
                   TABLE_RELS)
from ..cluster.data import ClusterMutsDataset
from ..core.batch import END5_COORD, END3_COORD, accum_fits
from ..core.dims import triangular
from ..core.header import ORDER_NAME, Header, make_header
from ..core.mu import (calc_p_ends_observed,
                       calc_p_noclose,
                       calc_p_noclose_given_ends,
                       calc_params)
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
        return calc_p_ends_observed(
            self.section.length,
            (self._end_counts.index.get_level_values(END5_COORD).values
             - self.section.end5),
            (self._end_counts.index.get_level_values(END3_COORD).values
             - self.section.end5),
            end_counts,
        )

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
                             self._num_reads,
                             self.section,
                             self.dataset.min_mut_gap)

    @cached_property
    def table_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        n_rels, n_clust = self._adjusted
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
    def clust_header(self):
        """ Header of the per-cluster data. """
        return make_header(max_order=self.max_order)

    @cached_property
    def table_per_clust(self):
        """ Number of reads in each cluster. """
        n_rels, n_clust = self._adjusted
        n_clust.name = "Number of Reads"
        return n_clust


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


def _get_clusters(n_reads_clust: pd.Series | int):
    """ Determine the clusters. """
    if isinstance(n_reads_clust, int):
        return None
    if isinstance(n_reads_clust, pd.Series):
        return n_reads_clust.index
    raise TypeError("n_reads_clust must be an int or Series, "
                    f"but got {type(n_reads_clust).__name__}")


def _calc_n_reads(n_reads_clust: pd.Series | int):
    """ Total number of reads among all clusters. """
    if isinstance(n_reads_clust, int):
        return n_reads_clust
    if isinstance(n_reads_clust, pd.Series):
        # Calculate the number of reads for each order.
        n_reads_order = n_reads_clust.groupby(level=ORDER_NAME).sum()
        # The numbers of reads should be identical.
        n_reads_order_1 = int(n_reads_order.loc[1])
        if not np.allclose(n_reads_order, n_reads_order_1):
            raise ValueError("Numbers of reads per order should match, "
                             f"but got {n_reads_order}")
        return n_reads_order_1
    raise TypeError("n_reads_clust must be an int or Series, "
                    f"but got {type(n_reads_clust).__name__}")


def _insert_masked(p_mut: pd.Series | pd.DataFrame,
                   section: Section):
    """ 2D array where masked positions are filled with 0. """
    # Fill masked positions with 0.
    p_mut = p_mut.reindex(index=section.range, fill_value=0.)
    # Convert to a 2D NumPy array.
    return p_mut.values.reshape((section.length, -1))


def _order_indices(order: int):
    """ First and last indices of the order in the array. """
    last = triangular(order)
    first = last - order
    return first, last


def adjust_counts(table_per_pos: pd.DataFrame,
                  p_ends_given_noclose: np.ndarray,
                  n_reads_clust: pd.Series | int,
                  section: Section,
                  min_mut_gap: int):
    """ Adjust the given table of masked/clustered counts per position
    to correct for observer bias.

    Parameters
    ----------
    p_mut_given_noclose: DataFrame
        Counts of the bits for each type of relation (column) at each
        position (index) in the section of interest.
    section: Section
        The section of interest.
    min_mut_gap: int
        Minimum number of non-mutated bases permitted between mutations.
    """
    # Determine which positions are unmasked.
    unmask = section.unmasked_bool
    # Calculate the number of reads with no two mutations too close.
    n_reads_noclose = _calc_n_reads(n_reads_clust)
    # Calculate the fraction of reads with no two mutations too close in
    # each cluster.
    p_clust_given_noclose = np.atleast_1d(n_reads_clust / n_reads_noclose)
    if p_clust_given_noclose.ndim != 1:
        raise ValueError("p_clust_given_noclose must have 1 dimension, "
                         f"but got {p_clust_given_noclose.ndim}")
    # Calculate the fraction of mutations at each position among reads
    # with no two mutations too close.
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # resulting NaN values are zeroed by nan_to_num.
        p_mut_given_noclose = np.nan_to_num(_insert_masked(
            table_per_pos[MUTAT_REL] / table_per_pos[INFOR_REL],
            section
        ))
    # Determine the clusters.
    clusters = _get_clusters(n_reads_clust)
    # Calculate the parameters.
    if clusters is None:
        # Calculate the parameters.
        p_mut, p_ends, p_clust = calc_params(
            p_mut_given_noclose,
            p_ends_given_noclose,
            p_clust_given_noclose,
            min_mut_gap
        )
        # Compute the probability that reads would have no two mutations
        # too close.
        p_noclose_given_clust = calc_p_noclose(
            p_ends,
            calc_p_noclose_given_ends(p_mut, min_mut_gap)
        )
        # Drop the cluster dimension from the parameters.
        if p_mut.shape != (section.length, 1):
            raise ValueError(f"p_mut must have shape {(section.length, 1)}, "
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
        p_noclose_given_clust = float(p_noclose_given_clust[0])
    else:
        # Calculate the parameters for each order separately.
        p_mut = np.empty_like(p_mut_given_noclose)
        p_clust = np.empty_like(p_clust_given_noclose)
        p_noclose_given_clust = np.empty_like(p_clust_given_noclose)
        for order in clusters.get_level_values(ORDER_NAME):
            i, j = _order_indices(order)
            # Calculate the parameters for each cluster.
            p_mut[:, i: j], p_ends, p_clust[i: j] = calc_params(
                p_mut_given_noclose[:, i: j],
                p_ends_given_noclose[:, :, i: j],
                p_clust_given_noclose[i: j],
                min_mut_gap
            )
            # Compute the probability that reads from each cluster would
            # have no two mutations too close.
            p_noclose_given_clust[i: j] = calc_p_noclose(
                p_ends,
                calc_p_noclose_given_ends(p_mut[:, i: j], min_mut_gap)
            )
    # Remove masked positions from the mutation rates.
    p_mut = p_mut[unmask]
    # Create the table of adjusted counts.
    # Initialize an empty DataFrame of the adjusted counts with the same
    # index and columns as the observed counts.
    n_rels = pd.DataFrame(np.nan, table_per_pos.index, table_per_pos.columns)
    # Calculate the probability that a read from any cluster would have
    # no two mutations too close.
    p_noclose = np.vdot(p_noclose_given_clust, p_clust)
    # Calculate the total number of reads.
    n_reads_total = n_reads_noclose / p_noclose
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
    with np.errstate(divide="ignore"):
        # Division by 0 is possible if no mutations were observed at a
        # given position, resulting in a NaN value at that position.
        scale = n_mut / table_per_pos.loc[unmask, MUTAT_REL].values
    # Replace NaN values with 1 so that missing values do not propagate
    # during multiplication.
    scale = np.nan_to_num(scale, nan=1.)
    # Scale every subtype of mutation by this factor.
    for mut in SUBMUTS:
        n_rels.loc[unmask, mut] = scale * table_per_pos.loc[unmask, mut].values
    # Calculate the number of reads in each cluster.
    n_clust = pd.Series(p_clust * n_reads_total, index=clusters)
    return n_rels, n_clust


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
