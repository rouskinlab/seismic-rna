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
from ..clust.data import ClustMerger
from ..core.batch import accum_fits
from ..core.header import Header, make_header
from ..core.mu import calc_f_obs_frame, calc_mu_adj_frame
from ..core.rel import RelPattern, HalfRelPattern
from ..core.seq import Section
from ..mask.data import MaskMerger
from ..relate.data import RelateLoader

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

    def __init__(self, dataset: RelateLoader | MaskMerger | ClustMerger):
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
        num_reads, per_pos, per_read = self._counts
        return num_reads

    @property
    def _counts_per_pos(self):
        """ Raw counts per position. """
        num_reads, per_pos, per_read = self._counts
        return per_pos

    @property
    def _counts_per_read(self):
        """ Raw counts per read. """
        num_reads, per_pos, per_read = self._counts
        return per_read

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
                             self.section,
                             self.dataset.min_mut_gap)

    @cached_property
    def table_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        counts_adj, f_obs = self._adjusted
        return counts_adj


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
        counts_adj, f_obs = self._adjusted
        return self.num_reads_obs / f_obs

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


def list_mut_patterns():
    return list(_iter_mut_patterns())


def _iter_patterns(mask: RelPattern | None = None):
    """ Yield a RelPattern for every type of relationship. """
    # Count everything except for no coverage.
    yield COVER_REL, RelPattern.allc()
    # Count matches to the reference sequence.
    yield MATCH_REL, RelPattern.muts().intersect(mask, invert=True)
    # Count all types of mutations, relative to reference matches.
    yield MUTAT_REL, RelPattern.muts().intersect(mask)
    # Count each type of mutation, relative to reference matches.
    for mut, pattern in list_mut_patterns():
        yield mut, RelPattern(pattern, HalfRelPattern.refs()).intersect(mask)


def all_patterns(mask: RelPattern | None = None):
    """ Every RelPattern, keyed by its name. """
    return dict(_iter_patterns(mask))


def adjust_counts(counts_obs: pd.DataFrame,
                  section: Section,
                  min_mut_gap: int):
    """
    Adjust the given table of masked/clustered bit counts per position
    to correct for observer bias. The table is mutated in-place.

    Parameters
    ----------
    counts_obs: DataFrame
        Counts of the bits for each type of relation (column) at each
        position (index) in the section of interest.
    section: Section
        The section of interest.
    min_mut_gap: int
        Minimum number of non-mutated bases permitted between mutations.
    """
    # Initialize an empty DataFrame of the adjusted counts with the same
    # index and columns as the observed counts.
    counts_adj = pd.DataFrame(np.nan,
                              index=counts_obs.index,
                              columns=counts_obs.columns)
    # Compute the observed fraction of mutations at each position.
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # NaN values will be zeroed out subsequently.
        fmuts_obs = counts_obs[MUTAT_REL] / counts_obs[INFOR_REL]
    # Fill any missing (NaN) values with zero.
    fmuts_obs.fillna(0., inplace=True)
    # Adjust the fraction of mutations to correct the observer bias.
    fmuts_adj = calc_mu_adj_frame(fmuts_obs, section, min_mut_gap)
    # Compute the fraction of all reads that would be observed.
    f_obs = calc_f_obs_frame(fmuts_adj, section, min_mut_gap)
    # Assume that the observance bias affects the counts of covered and
    # informative bases equally, so that we can define f_obs as follows:
    # f_obs := ninfo_obs / ninfo_adj
    # from which we can estimate the informative bases after adjustment:
    # ninfo_adj = ninfo_obs / f_obs
    ninfo_adj = counts_obs.loc[:, INFOR_REL] / f_obs
    counts_adj.loc[:, INFOR_REL] = ninfo_adj.values
    counts_adj.loc[:, COVER_REL] = (counts_obs.loc[:, COVER_REL] / f_obs).values
    # From the definition of the adjusted fraction of mutations:
    # fmuts_adj := nmuts_adj / ninfo_adj
    # we can also estimate the mutated bases after adjustment:
    # nmuts_adj = fmuts_adj * ninfo_adj
    nmuts_adj = fmuts_adj * ninfo_adj
    counts_adj.loc[:, MUTAT_REL] = nmuts_adj.values
    # From the definition of informative bases:
    # ninfo_adj := nrefs_adj + nmuts_adj
    # we can estimate the matched bases after adjustment:
    # nrefs_adj = ninfo_adj - nmuts_adj
    nrefs_adj = ninfo_adj - nmuts_adj
    counts_adj.loc[:, MATCH_REL] = nrefs_adj.values
    # Compute the factor by which nmuts was adjusted for each position.
    with np.errstate(divide="ignore"):
        # Division by 0 is possible if no mutations were observed at a
        # given position, resulting in a NaN value at that position.
        adj_factor = nmuts_adj / counts_obs.loc[:, MUTAT_REL]
    # Replace NaN values with 1 so that missing values do not propagate
    # during multiplication.
    adj_factor = adj_factor.fillna(1.)
    # Adjust every subtype of mutation by this factor.
    for mut in SUBMUTS:
        counts_adj.loc[:, mut] = (counts_obs.loc[:, mut] * adj_factor).values
    logger.debug(f"Adjusted counts\n\nfrom\n{counts_obs}\n\nto\n{counts_adj}")
    return counts_adj, f_obs


def tabulate_loader(dataset: RelateLoader | MaskMerger | ClustMerger):
    """ Return a new Dataset, choosing the subclass based on the type
    of the argument `dataset`. """
    if isinstance(dataset, RelateLoader):
        return RelateTabulator(dataset)
    if isinstance(dataset, MaskMerger):
        return MaskTabulator(dataset)
    if isinstance(dataset, ClustMerger):
        return ClustTabulator(dataset)
    raise TypeError(f"Invalid dataset type: {type(dataset).__name__}")

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
