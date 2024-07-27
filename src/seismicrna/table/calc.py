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
                   UNAMB_REL,
                   TABLE_RELS)
from ..cluster.data import ClusterMutsDataset, load_cluster_dataset
from ..core.batch import END5_COORD, END3_COORD, accum_fits
from ..core.data import MutsDataset, UnbiasDataset
from ..core.header import NUM_CLUSTS_NAME, Header, make_header, validate_ks
from ..core.mu import (calc_p_ends_observed,
                       calc_p_noclose,
                       calc_p_noclose_given_ends,
                       calc_params)
from ..core.rel import RelPattern, HalfRelPattern
from ..core.seq import Section
from ..mask.data import load_mask_dataset
from ..pool.data import load_relate_dataset

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
        table.loc[:, UNAMB_REL] = (table.loc[:, MATCH_REL].values
                                   +
                                   table.loc[:, MUTAT_REL].values)
        return table

    def __init__(self,
                 dataset: MutsDataset | ClusterMutsDataset | UnbiasDataset):
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
    def ks(self) -> list[int] | None:
        """ Numbers of clusters. """

    @cached_property
    def pos_header(self):
        """ Header of the per-position data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def read_header(self):
        """ Header of the per-read data. """
        return make_header(rels=TABLE_RELS, ks=self.ks)

    @cached_property
    def _counts(self):
        return accum_fits(self.dataset.iter_batches(),
                          self.refseq,
                          self.section.unmasked_int,
                          all_patterns(self.dataset.pattern),
                          ks=self.ks)

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
        end5s = (self._end_counts.index.get_level_values(END5_COORD).values
                 - self.section.end5)
        end3s = (self._end_counts.index.get_level_values(END3_COORD).values
                 - self.section.end5)
        return calc_p_ends_observed(self.section.length,
                                    end5s,
                                    end3s,
                                    end_counts)

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
        table_per_pos = super().table_per_pos
        if self.dataset.min_mut_gap > 0:
            if self.section.length > np.sqrt(1.e9):
                logger.warning("Using bias correction on a section with "
                               f"{self.section.length} positions requires "
                               ">1 GB of memory. If this is impractical, you "
                               "can (at the cost of lower accuracy) disable "
                               "bias correction using --min-mut-gap 0.")
            try:
                return adjust_counts(table_per_pos,
                                     self.p_ends_given_noclose,
                                     self._num_reads,
                                     self.section,
                                     self.dataset.min_mut_gap,
                                     self.dataset.quick_unbias,
                                     self.dataset.quick_unbias_thresh)
            except Exception as error:
                logger.warning(f"Bias correction failed: {error}")
        return table_per_pos, self._num_reads

    @cached_property
    def table_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        n_rels, n_clust = self._adjusted
        return n_rels


class AvgTabulator(Tabulator, ABC):

    @property
    def ks(self):
        return None


class RelateTabulator(FullTabulator, AvgTabulator):
    pass


class MaskTabulator(PartialTabulator, AvgTabulator):
    pass


class ClustTabulator(PartialTabulator):

    @property
    def ks(self):
        ks = self.dataset.ks
        if not ks:
            raise ValueError(f"Cannot tabulate {self.dataset} with no clusters")
        return ks

    @cached_property
    def clust_header(self):
        """ Header of the per-cluster data. """
        return make_header(ks=self.ks)

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


def _insert_masked(p_mut: pd.Series | pd.DataFrame,
                   section: Section):
    """ 2D array where masked positions are filled with 0. """
    # Fill masked positions with 0.
    p_mut = p_mut.reindex(index=section.range, fill_value=0.)
    # Convert to a 2D NumPy array.
    return p_mut.values.reshape((section.length, -1))


def adjust_counts(table_per_pos: pd.DataFrame,
                  p_ends_given_noclose: np.ndarray,
                  n_reads_clust: pd.Series | int,
                  section: Section,
                  min_mut_gap: int,
                  quick_unbias: bool,
                  quick_unbias_thresh: float):
    """ Adjust the given table of masked/clustered counts per position
    to correct for observer bias. """
    # Determine which positions are unmasked.
    unmask = section.unmasked_bool
    # Calculate the fraction of mutations at each position among reads
    # with no two mutations too close.
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # resulting NaN values are zeroed by nan_to_num.
        p_mut_given_noclose = np.nan_to_num(_insert_masked(
            table_per_pos[MUTAT_REL] / table_per_pos[UNAMB_REL],
            section
        ))
    if isinstance(n_reads_clust, int):
        # There is only one cluster, so the probability that each read
        # belongs to that cluster is 1.
        p_clust_given_noclose = np.array([1.])
        # Calculate the parameters.
        p_mut, p_ends, p_clust = calc_params(
            p_mut_given_noclose,
            p_ends_given_noclose,
            p_clust_given_noclose,
            min_mut_gap,
            quick_unbias=quick_unbias,
            quick_unbias_thresh=quick_unbias_thresh
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
                p_ends_given_noclose[:, :, ki],
                p_clust_given_noclose,
                min_mut_gap
            )
            # Compute the probability that reads from each cluster would
            # have no two mutations too close.
            p_noclose_given_clust[ki] = calc_p_noclose(
                p_ends,
                calc_p_noclose_given_ends(p_mut[:, ki], min_mut_gap)
            )
            # Compute the probability that reads from any cluster would
            # have no two mutations too close.
            p_noclose = float(np.vdot(p_noclose_given_clust[ki], p_clust[ki]))
            # Compute the number of reads in each cluster.
            n_clust.loc[k] = (n_reads_noclose / p_noclose) * p_clust[ki]
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
    n_info = table_per_pos.loc[unmask, UNAMB_REL].values / p_noclose_given_clust
    n_rels.loc[unmask, COVER_REL] = n_cov
    n_rels.loc[unmask, UNAMB_REL] = n_info
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
    return n_rels, n_clust


def tabulate_loader(dataset: MutsDataset | ClusterMutsDataset | UnbiasDataset):
    """ Return a new Dataset, choosing the subclass based on the type
    of the argument `dataset`. """
    if load_relate_dataset.is_dataset_type(dataset):
        return RelateTabulator(dataset)
    if load_mask_dataset.is_dataset_type(dataset):
        return MaskTabulator(dataset)
    if load_cluster_dataset.is_dataset_type(dataset):
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
