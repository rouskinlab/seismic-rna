from abc import ABC, abstractmethod
from functools import cache, cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .base import (TOTAL_REL, DELET_REL, INSRT_REL, MATCH_REL, MUTAT_REL,
                   SUBST_REL, SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL,
                   CLUST_INDEX_NAMES, REL_NAME, R_ADJ_TITLE, R_OBS_TITLE,
                   REL_CODES)
from ..cluster.names import READ_NAME
from ..cluster.load import ClustLoader
from ..core.bitcall import BitCaller, SemiBitCaller
from ..core.bitvect import BitCounter, ClustBitCounter
from ..core.mu import calc_f_obs_df, calc_mu_adj_df
from ..core.sect import Section
from ..mask.load import MaskLoader
from ..relate.load import RelateLoader

logger = getLogger(__name__)


# Tabulator Classes ####################################################

class Tabulator(ABC):
    """ Base class for tabulating data for multiple tables from a report
    loader. """

    def __init__(self, loader: RelateLoader | MaskLoader | ClustLoader,
                 rel_codes: str):
        self._loader = loader
        self._rel_codes = rel_codes

    @property
    def out_dir(self):
        """ Output directory. """
        return self._loader.out_dir

    @property
    def sample(self):
        """ Name of the sample. """
        return self._loader.sample

    @property
    def ref(self):
        """ Name of the reference. """
        return self._loader.ref

    @property
    def sect(self):
        """ Name of the section covered by the table. """
        return self._loader.sect

    @property
    def section(self):
        """ Section covered by the table. """
        return self._loader.section

    def iter_bit_callers(self):
        """ Yield a BitCaller for each relationship in the table. """
        yield from iter_bit_callers(self.section, self._rel_codes)

    @cached_property
    def rels(self):
        """ Relationships counted in the table. """
        return [rel for rel, _, _ in self.iter_bit_callers()]

    @property
    def seq_array(self):
        """ Array of the bases of the section at unmasked positions. """
        return self._loader.seq.to_str_array()

    @property
    @abstractmethod
    def columns(self):
        """ Columns of the table. """

    @cached_property
    @abstractmethod
    def bit_counts(self) -> dict[str, BitCounter]:
        """ Return BitCounters for the batches of relation vectors. """

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    def tabulate_by_pos(self):
        return pd.DataFrame.from_dict({rel: counter.n_affi_per_pos
                                       for rel, counter
                                       in self.bit_counts.items()})

    def tabulate_by_read(self):
        return pd.DataFrame.from_dict({rel: counter.n_affi_per_read
                                       for rel, counter
                                       in self.bit_counts.items()})


class RelTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return 0

    @property
    def columns(self):
        return pd.Index(self.rels, name=REL_NAME)

    @cached_property
    def bit_counts(self):
        bit_counts = dict()
        for name, bit_caller, kwargs in self.iter_bit_callers():
            if kwargs:
                bit_caller = BitCaller.inter(bit_caller, **kwargs)
            bit_counts[name] = BitCounter(
                self.section,
                bit_caller.iter(self._loader.iter_batches_processed())
            )
        return bit_counts


class MaskTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def columns(self):
        return pd.Index(self.rels, name=REL_NAME)

    @cached_property
    def bit_counts(self):
        return {
            rel: BitCounter(self.section,
                            self._loader.iter_batches_processed(bit_caller=bc,
                                                                **kwargs))
            for rel, bc, kwargs in self.iter_bit_callers()
        }

    def tabulate_by_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass.
        counts_obs = super().tabulate_by_pos()
        # Adjust the counts to correct for observer bias.
        return adjust_counts(counts_obs,
                             self._loader.section,
                             self._loader.min_mut_gap)


class ClustTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def clusters(self):
        """ Order and number of each cluster. """
        return self._loader.clusters

    def get_read_names(self):
        return self._loader.import_loader.get_read_names()

    @cached_property
    def columns(self):
        return pd.MultiIndex.from_tuples([(order, cluster, rel)
                                          for order, cluster in self.clusters
                                          for rel in self.rels],
                                         names=CLUST_INDEX_NAMES)

    @cached_property
    def bit_counts(self):
        return dict(
            (rel, ClustBitCounter(
                self.section,
                self.clusters,
                self._loader.iter_batches_processed(bit_caller=bc, **kwargs)
            )) for rel, bc, kwargs in self.iter_bit_callers()
        )

    @cache
    def tabulate_by_pos(self):
        """ DataFrame of the bit count for each position and caller. """
        # Initialize an empty DataFrame.
        counts_obs = pd.DataFrame(self.get_null_value(),
                                  index=self._loader.section.unmasked,
                                  columns=self.columns)
        # Fill in the DataFrame column by column.
        for col in self.columns:
            nk, k, rel = col
            counts_obs[col] = self.bit_counts[rel].n_affi_per_pos.loc[:, (nk, k)]
        # Adjust the counts to correct for observer bias.
        counts_adj = dict()
        for order, k in self.clusters:
            for rel, adj in adjust_counts(counts_obs.loc[:, (order, k)],
                                          self._loader.section,
                                          self._loader.min_mut_gap).items():
                counts_adj[order, k, rel] = adj
        counts_adj = pd.DataFrame.from_dict(counts_adj)
        counts_adj.columns.rename(counts_obs.columns.names, inplace=True)
        return counts_adj

    def tabulate_by_read(self):
        """ DataFrame of the bit count for each read and caller. """
        # Initialize an empty DataFrame.
        counts = pd.DataFrame(self.get_null_value(),
                              index=pd.Index(self.get_read_names(),
                                             name=READ_NAME),
                              columns=self.columns)
        # Fill in the DataFrame column by column.
        for col in self.columns:
            # Split the column label into the clustering order, cluster
            # number, and type of relationship.
            order, k, rel = col
            # Get the bit counter for the relationship and select the
            # column for the clustering order and cluster number.
            counts[col] = self.bit_counts[rel].n_affi_per_read.loc[:, (order, k)]
        return counts

    def tabulate_by_clust(self):
        """ Return the adjusted number of reads in each cluster as a
        Series with dimension (clusters). """
        return pd.DataFrame.from_dict({
            R_OBS_TITLE: self._loader.n_reads_obs,
            R_ADJ_TITLE: self._loader.n_reads_adj,
        })


# Helper functions #####################################################


def iter_mut_semi_callers():
    """ Yield a SemiBitCaller for each type of mutation to tabulate. """
    yield SUBST_REL, SemiBitCaller.from_counts(count_sub=True)
    yield SUB_A_REL, SemiBitCaller("ca", "ga", "ta")
    yield SUB_C_REL, SemiBitCaller("ac", "gc", "tc")
    yield SUB_G_REL, SemiBitCaller("ag", "cg", "tg")
    yield SUB_T_REL, SemiBitCaller("at", "ct", "gt")
    yield DELET_REL, SemiBitCaller.from_counts(count_del=True)
    yield INSRT_REL, SemiBitCaller.from_counts(count_ins=True)


def _iter_bit_callers(section: Section):
    """ Yield a BitCaller for every type of relationship. """
    # Call reference matches.
    refc = SemiBitCaller.from_counts(count_ref=True)
    # Call mutations.
    mutc = SemiBitCaller.from_counts(count_sub=True,
                                     count_del=True,
                                     count_ins=True)
    # Create a standard bit caller for all matches and mutations.
    bitc = BitCaller(section, mutc, refc)
    # Count all base calls (everything but the bytes 0 and 255).
    yield TOTAL_REL, bitc, dict(merge=True)
    # Count matches to the reference sequence.
    yield MATCH_REL, bitc, dict(invert=True)
    # Count all types of mutations, relative to reference matches.
    yield MUTAT_REL, bitc, dict()
    # Count each type of mutation, relative to reference matches.
    for mut, mutc in iter_mut_semi_callers():
        yield mut, BitCaller(section, mutc, refc), dict()


def iter_bit_callers(section: Section, rel_codes: str):
    """ Yield a BitCaller for each type of relationship to tabulate. """
    # Determine which types of relationships to call.
    if rel_codes:
        # Use only the selected types of relationships.
        call_rels = {REL_CODES[code] for code in rel_codes}
    else:
        # Use all types of relationships.
        call_rels = set(REL_CODES.values())
    # Yield a BitCaller for each selected type of relationship.
    for rel, bit_caller, kwargs in _iter_bit_callers(section):
        if rel in call_rels:
            yield rel, bit_caller, kwargs


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
    # Compute the observed fraction of mutations at each position.
    nrefs_obs = counts_obs.loc[:, MATCH_REL]
    nmuts_obs = counts_obs.loc[:, MUTAT_REL]
    ninfo_obs = nrefs_obs + nmuts_obs
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # NaN values will be zeroed out subsequently.
        fmuts_obs = nmuts_obs / ninfo_obs
    # Fill missing NaN values with zero.
    fmuts_obs.fillna(0., inplace=True)
    # Adjust the fraction of mutations to correct the observer bias.
    mu_adj = calc_mu_adj_df(fmuts_obs.to_frame(),
                            section, min_mut_gap).squeeze(axis=1)
    # Compute the fraction of reads that would be observed.
    f_obs = float(calc_f_obs_df(mu_adj.to_frame(), section, min_mut_gap)[0])
    # Initialize new data for the adjusted counts.
    counts_adj = dict()
    # Assume that the total number of base calls including unobserved
    # reads is the number observed divided by the fraction of the total
    # reads that were observed.
    counts_adj[TOTAL_REL] = counts_obs.loc[:, TOTAL_REL] / f_obs
    # Two conditions must be met:
    # : The number of informative bases (matches + mutations) after
    #   adjustment must equal the number before adjustment divided
    #   by the fraction of reads that were observed:
    #   (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : The fraction of mutations at each position after adjustment
    #   must equal the adjusted number of mutations divided by the
    #   adjusted number of informative bases:
    #   nmuts_adj / (nrefs_adj + nmuts_adj) = mu_adj
    # The two unknown variables (nrefs_adj and nmuts_adj) can be
    # solved using the above system of two equations. By solving for
    # (nrefs_adj + nmuts_adj) in both equations, we get:
    # : (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : (nrefs_adj + nmuts_adj) = nmuts_adj / mu_adj, if mu_adj > 0
    # Setting both right hand sides equal, nmuts_adj is solved:
    # : (nrefs_obs + nmuts_obs) / f_obs = nmuts_adj / mu_adj
    # : nmuts_adj = (nrefs_obs + nmuts_obs) * (mu_adj / f_obs)
    nmuts_adj = ninfo_obs * (mu_adj / f_obs)
    # Then, plugging this solution into the second equation:
    # : nrefs_adj = nmuts_adj / mu_adj - nmuts_adj
    #             = nmuts_adj * (1 / mu_adj - 1), if mu_adj > 0
    # If, on the other hand, mu_adj = 0, then
    # : nmuts_adj = mu_adj * (nrefs_adj + nmuts_adj) = 0
    # Thus, we solve for nrefs_adj using the first equation instead:
    # : (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : nrefs_adj = nrefs_obs / f_obs
    with np.errstate(divide="ignore"):
        # Ignore division-by-zero warnings in np.reciprocal because
        # positions where mu_adj = 0 are just skipped by np.where.
        nrefs_adj = np.where(mu_adj > 0.,
                             nmuts_adj * (np.reciprocal(mu_adj) - 1.),
                             nrefs_obs / f_obs)
    counts_adj[MATCH_REL] = nrefs_adj
    counts_adj[MUTAT_REL] = nmuts_adj
    # Compute the factor by which nmuts was adjusted:
    nmuts_fac = nmuts_adj / nmuts_obs
    # Adjust every type of mutation by this factor.
    for mut in (SUBST_REL, SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL,
                DELET_REL, INSRT_REL):
        counts_adj[mut] = counts_obs.loc[:, mut] * nmuts_fac
    # Assemble the adjusted counts into a new DataFrame.
    return pd.DataFrame.from_dict(counts_adj)


def tabulate_loader(report_loader: RelateLoader | MaskLoader | ClustLoader,
                    rel_codes: str):
    """ Return a new DataLoader, choosing the subclass based on the type
    of the argument `report_loader`. """
    if isinstance(report_loader, RelateLoader):
        return RelTabulator(report_loader, rel_codes)
    if isinstance(report_loader, MaskLoader):
        return MaskTabulator(report_loader, rel_codes)
    if isinstance(report_loader, ClustLoader):
        return ClustTabulator(report_loader, rel_codes)
    raise TypeError(f"Invalid loader type: {type(report_loader).__name__}")
