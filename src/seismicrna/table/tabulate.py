from abc import ABC, abstractmethod
from functools import cache, cached_property
from logging import getLogger
from typing import Any

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
                   CLUST_INDEX_NAMES,
                   R_ADJ_TITLE,
                   R_OBS_TITLE,
                   READ_TITLE,
                   REL_NAME,
                   TABLE_RELS)
from ..cluster.load import ClustLoader
from ..core.batch import accum_fits
from ..core.mu import calc_f_obs_series, calc_mu_adj_series
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

    def __init__(self, loader: RelateLoader | MaskMerger | ClustLoader):
        self.loader = loader

    @property
    def top(self):
        return self.loader.top

    @property
    def sample(self):
        return self.loader.sample

    @property
    def ref(self):
        return self.loader.ref

    @property
    def refseq(self):
        return self.loader.refseq

    @property
    def section(self):
        return self.loader.section

    @property
    @abstractmethod
    def columns(self):
        """ Columns of the table. """

    @cached_property
    @abstractmethod
    def _counts(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """ Raw counts per position and per read. """

    @property
    def _counts_per_pos(self):
        """ Raw counts per position. """
        per_pos, per_read = self._counts
        return per_pos

    @property
    def _counts_per_read(self):
        """ Raw counts per read. """
        per_pos, per_read = self._counts
        return per_read

    @cached_property
    @abstractmethod
    def table_per_pos(self):
        """ Count relationships at every position in the section. """

    @cached_property
    @abstractmethod
    def table_per_read(self):
        """ Count relationships in every read. """


class EnsembleTabulator(Tabulator, ABC):

    @classmethod
    def columns(cls):
        return pd.Index(TABLE_RELS, name=REL_NAME)

    @classmethod
    def _tabulate(cls, counts: pd.DataFrame):
        """ Counts with the proper columns. """
        # Initialize an empty table.
        table = pd.DataFrame(cls.get_null_value(), counts.index, cls.columns())
        # Count reads with each relationship at each position.
        for rel, rel_counts in counts.items():
            table[rel] = rel_counts
        # Add a column for informative relationships.
        table[INFOR_REL] = table[MATCH_REL] + table[MUTAT_REL]
        return table

    @cached_property
    def _counts(self):
        return accum_fits(self.section.unmasked_int,
                          self.refseq,
                          name_patterns(self.loader.pattern),
                          self.loader.iter_batches())

    @cached_property
    def table_per_pos(self):
        return self._tabulate(self._counts_per_pos).reindex(self.section.range)

    @cached_property
    def table_per_read(self):
        return self._tabulate(self._counts_per_read)


class NullableTabulator(Tabulator, ABC):

    @classmethod
    def get_null_value(cls):
        return np.nan


class RelTabulator(EnsembleTabulator):

    @classmethod
    def get_null_value(cls):
        return 0


class MaskTabulator(EnsembleTabulator, NullableTabulator):

    @cached_property
    def table_per_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass, then adjust for observer bias.
        return adjust_counts(super().table_per_pos,
                             self.section,
                             self.loader.min_mut_gap)


class ClusterTabulator(NullableTabulator):

    @property
    def ord_clust(self):
        """ Order and number of each cluster. """
        return self.loader.clusters

    @cached_property
    def columns(self):
        return pd.MultiIndex.from_tuples([(order, cluster, rel)
                                          for order, cluster in self.ord_clust
                                          for rel in TABLE_RELS],
                                         names=CLUST_INDEX_NAMES)

    @cached_property
    def bit_counts(self):
        return dict(
            (rel, ClustBitCounter(
                self.section,
                self.ord_clust,
                self._loader.iter_batches_processed(bit_caller=bc, **kwargs)
            )) for rel, bc, kwargs in self.iter_bit_callers()
        )

    @staticmethod
    def _slice_rel(counts: pd.DataFrame, rel: str, values: Any | None = None):
        """ Slice one relationship from all clusters. """
        slicer = slice(None), slice(None), rel
        if values is None:
            return counts.loc[:, slicer]
        counts.loc[:, slicer] = values

    @cache
    def table_per_pos(self):
        """ DataFrame of the bit count for each position and caller. """
        # Initialize an empty DataFrame of observed counts.
        counts = pd.DataFrame(self.get_null_value(),
                              index=self.section.unmasked,
                              columns=self.columns)
        # Fill the observed counts one relationship at a time.
        for rel in TABLE_RELS:
            try:
                counter = self.bit_counts[rel]
            except KeyError:
                # The relationship was not among those tabulated.
                pass
            else:
                self._slice_rel(counts, rel, counter.n_affi_per_pos.values)
        # Count the informative bits.
        self._slice_rel(counts, INFOR_REL,
                        self._slice_rel(counts, MATCH_REL).values
                        + self._slice_rel(counts, MUTAT_REL).values)
        # Adjust the counts for each cluster.
        for ok in self.ord_clust:
            # Adjust the counts to correct for observer bias.
            counts.loc[:, ok] = adjust_counts(counts.loc[:, ok],
                                              self.section,
                                              self.loader.min_mut_gap).values
        return counts

    def table_per_read(self):
        """ DataFrame of the bit count for each read and caller. """

        # Initialize an empty DataFrame.
        counts = pd.DataFrame(self.get_null_value(),
                              index=pd.Index(self.get_read_names(),
                                             name=READ_TITLE),
                              columns=self.columns)
        # Fill in the DataFrame one relationship at a time.
        for rel in TABLE_RELS:
            try:
                counter = self.bit_counts[rel]
            except KeyError:
                # The relationship was not among those tabulated.
                pass
            else:
                self._slice_rel(counts, rel, counter.n_affi_per_read.values)
        # Count the informative bits.
        self._slice_rel(counts, INFOR_REL,
                        self._slice_rel(counts, MATCH_REL).values
                        + self._slice_rel(counts, MUTAT_REL).values)
        return counts

    def tabulate_by_clust(self):
        """ Return the adjusted number of reads in each cluster as a
        Series with dimension (clusters). """
        return pd.DataFrame.from_dict({
            R_OBS_TITLE: self.loader.n_reads_obs,
            R_ADJ_TITLE: self.loader.n_reads_adj,
        })


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


@cache
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


@cache
def name_patterns(mask: RelPattern | None = None):
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
    logger.debug(f"Adjusting mutation counts of table\n{counts_obs}")
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
    fmuts_adj = calc_mu_adj_series(fmuts_obs, section, min_mut_gap)
    # Compute the fraction of all reads that would be observed.
    f_obs = calc_f_obs_series(fmuts_adj, section, min_mut_gap)
    # Assume that the observance bias affects the counts of covered and
    # informative bases equally, so that we can define f_obs as follows:
    # f_obs := ninfo_obs / ninfo_adj
    # from which we can estimate the informative bases after adjustment:
    # ninfo_adj = ninfo_obs / f_obs
    ninfo_adj = counts_obs.loc[:, INFOR_REL] / f_obs
    counts_adj.loc[:, INFOR_REL] = ninfo_adj
    counts_adj.loc[:, COVER_REL] = counts_obs.loc[:, COVER_REL] / f_obs
    # From the definition of the adjusted fraction of mutations:
    # fmuts_adj := nmuts_adj / ninfo_adj
    # we can also estimate the mutated bases after adjustment:
    # nmuts_adj = fmuts_adj * ninfo_adj
    nmuts_adj = fmuts_adj * ninfo_adj
    counts_adj.loc[:, MUTAT_REL] = nmuts_adj
    # From the definition of informative bases:
    # ninfo_adj := nrefs_adj + nmuts_adj
    # we can estimate the matched bases after adjustment:
    # nrefs_adj = ninfo_adj - nmuts_adj
    nrefs_adj = ninfo_adj - nmuts_adj
    counts_adj.loc[:, MATCH_REL] = nrefs_adj
    # Compute the factor by which nmuts was adjusted for each position.
    with np.errstate(divide="ignore"):
        # Division by 0 is possible if no mutations were observed at a
        # given position, resulting in a NaN value at that position.
        adj_factor = nmuts_adj / counts_obs.loc[:, MUTAT_REL]
    # So that the multiplication works properly, replace NaN values with
    # 1 so that they do not propagate, then cast to a DataFrame so that
    # the dimensions are correct.
    adj_factor = adj_factor.fillna(1.).to_frame().values
    # Adjust every subtype of mutation by this factor.
    counts_adj.loc[:, SUBMUTS] = counts_obs.loc[:, SUBMUTS] * adj_factor
    return counts_adj


def tabulate_loader(dataset: RelateLoader | MaskMerger | ClustLoader):
    """ Return a new DataLoader, choosing the subclass based on the type
    of the argument `dataset`. """
    if isinstance(dataset, RelateLoader):
        return RelTabulator(dataset)
    if isinstance(dataset, MaskMerger):
        return MaskTabulator(dataset)
    if isinstance(dataset, ClustLoader):
        return ClusterTabulator(dataset)
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
