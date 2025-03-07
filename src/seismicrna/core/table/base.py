from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any, Generator, Iterable

import numpy as np
import pandas as pd

from .. import path
from ..batch import RB_INDEX_NAMES
from ..dataset import LoadFunction
from ..header import REL_NAME, Header, parse_header
from ..mu import winsorize
from ..rel import HalfRelPattern, RelPattern
from ..rna import RNAProfile
from ..seq import DNA, SEQ_INDEX_NAMES, Region, index_to_pos, index_to_seq

# General fields
READ_TITLE = "Read Name"

# Count relationships
COVER_REL = "Covered"
INFOR_REL = "Informative"
MATCH_REL = "Matched"
MUTAT_REL = "Mutated"
DELET_REL = "Deleted"
INSRT_REL = "Inserted"
SUBST_REL = "Subbed"
SUB_A_REL = "Subbed-A"
SUB_C_REL = "Subbed-C"
SUB_G_REL = "Subbed-G"
SUB_T_REL = "Subbed-T"

# One-letter codes for each type of relationship
REL_CODES = {"v": COVER_REL,
             "n": INFOR_REL,
             "e": MATCH_REL,
             "m": MUTAT_REL,
             "s": SUBST_REL,
             "a": SUB_A_REL,
             "c": SUB_C_REL,
             "g": SUB_G_REL,
             "t": SUB_T_REL,
             "d": DELET_REL,
             "i": INSRT_REL}
REL_NAMES = {name: code for code, name in REL_CODES.items()}

# Columns of each relation-based table
TABLE_RELS = list(REL_CODES.values())


@cache
def get_pattern(rel: str):
    """ Get a RelPattern from the name of its relationship. """
    if rel == COVER_REL:
        return RelPattern.allc()
    if rel == MATCH_REL:
        return RelPattern.matches()
    if rel == MUTAT_REL:
        return RelPattern.muts()
    if rel == SUBST_REL:
        yes = HalfRelPattern.from_counts(count_sub=True)
    elif rel == SUB_A_REL:
        yes = HalfRelPattern("ca", "ga", "ta")
    elif rel == SUB_C_REL:
        yes = HalfRelPattern("ac", "gc", "tc")
    elif rel == SUB_G_REL:
        yes = HalfRelPattern("ag", "cg", "tg")
    elif rel == SUB_T_REL:
        yes = HalfRelPattern("at", "ct", "gt")
    elif rel == DELET_REL:
        yes = HalfRelPattern.from_counts(count_del=True)
    elif rel == INSRT_REL:
        yes = HalfRelPattern.from_counts(count_ins=True)
    else:
        raise ValueError(rel)
    return RelPattern(yes, HalfRelPattern.refs())


def get_subpattern(rel: str, subpattern: RelPattern | None = None):
    """ Get a RelPattern, optionally with masking. """
    if rel == COVER_REL:
        # Still match everything except the non-covered relationship.
        return get_pattern(rel)
    if rel == MATCH_REL and subpattern is not None:
        # Because select_muts specifies the types of mutations, it must
        # be inverted before intersecting with matches.
        subpattern = subpattern.invert()
    # All other relationships are mutations needing no special handling.
    return get_pattern(rel).intersect(subpattern)


def all_patterns(subpattern: RelPattern | None = None):
    """ Every RelPattern, keyed by its name. """
    return {rel: get_subpattern(rel, subpattern)
            for rel in REL_CODES.values()
            if rel != INFOR_REL}


def get_rel_name(rel_code: str):
    """ Get the name of a relationship from its code. """
    try:
        return REL_CODES[rel_code]
    except KeyError:
        raise ValueError(f"Relationship code must be one of {list(REL_CODES)}, "
                         f"but got {repr(rel_code)}")


def _get_denom_rel(rel: str):
    """ Get the relationship that serves as the denominator. """
    return COVER_REL if rel == COVER_REL or rel == INFOR_REL else INFOR_REL


def _get_denom_cols(numer_cols: pd.Index):
    """ Get the denominator columns based on the numerator columns. """
    if isinstance(numer_cols, pd.MultiIndex):
        return pd.MultiIndex.from_arrays(
            [list(map(_get_denom_rel, numer_cols.get_level_values(name)))
             if name == REL_NAME
             else numer_cols.get_level_values(name)
             for name in numer_cols.names],
            names=numer_cols.names
        )
    if numer_cols.name is not None and numer_cols.name != REL_NAME:
        raise ValueError(f"Expected index to be named {repr(REL_NAME)}, "
                         f"but got {repr(numer_cols.name)}")
    return pd.Index(list(map(_get_denom_rel, numer_cols)),
                    name=numer_cols.name)


class Table(path.HasRefFilePath, ABC):
    """ Table base class. """

    @classmethod
    @abstractmethod
    def get_load_function(cls) -> LoadFunction:
        """ LoadFunction for all Dataset types for this Table. """

    @classmethod
    @abstractmethod
    def get_by_read(cls) -> bool:
        """ Whether the table contains data for each read. """

    @classmethod
    def get_ext(cls):
        return path.CSVZIP_EXT if cls.get_by_read() else path.CSV_EXT

    @classmethod
    @cache
    def get_auto_path_fields(cls):
        return {path.STEP: cls.get_step(),
                path.TABLE: cls.get_step(),
                **super().get_auto_path_fields()}

    @classmethod
    @abstractmethod
    def get_header_type(cls) -> type[Header]:
        """ Type of the header for the table. """

    @classmethod
    def get_header_depth(cls):
        return cls.get_header_type().get_num_levels()

    @classmethod
    @abstractmethod
    def get_index_depth(cls) -> int:
        """ Number of columns in the index. """

    @classmethod
    def get_index_cols(cls) -> list[int]:
        """ Column(s) of the file to use as the index. """
        return list(range(cls.get_index_depth()))

    @property
    @abstractmethod
    def _attrs(self):
        """ Source of the table's attributes. """

    @property
    def top(self) -> Path:
        """ Path of the table's output directory. """
        return self._attrs.top

    @property
    def branches(self) -> dict[str, str]:
        """ Branches of the workflow. """
        return self._attrs.branches

    @property
    def sample(self) -> str:
        """ Name of the table's sample. """
        return self._attrs.sample

    @property
    def ref(self) -> str:
        """ Name of the table's reference. """
        return self._attrs.ref

    @property
    def reg(self) -> str:
        """ Name of the table's region. """
        return self._attrs.region.name

    @cached_property
    def refseq(self) -> DNA:
        """ Reference sequence. """
        return self._attrs.refseq

    @cached_property
    def path(self):
        """ Path of the table's file. """
        field_values = dict()
        for field in path.get_fields_in_seg_types(self.get_path_seg_types(),
                                                  include_top=True):
            try:
                field_values[field] = getattr(self, field)
            except AttributeError:
                # If the field name is not an attribute of the table,
                # then it has a default value that will be supplied by
                # self.build_path().
                pass
        return self.build_path(field_values)

    @abstractmethod
    def _get_header(self) -> Header:
        """ Get the Header for the table's data. """

    @cached_property
    def header(self):
        """ Header for the table's data. """
        header = self._get_header()
        if not isinstance(header, self.get_header_type()):
            raise TypeError(f"Expected {self.get_header_type().__name__}, "
                            f"but got {type(header).__name__}")
        return header

    @cached_property
    @abstractmethod
    def data(self) -> pd.DataFrame | pd.Series:
        """ Table's data. """

    def __str__(self):
        return f"{type(self).__name__}: {self.path}"


class RelTypeTable(Table, ABC):
    """ Table with multiple types of relationships. """

    @classmethod
    def get_header_rows(cls) -> list[int]:
        """ Row(s) of the file to use as the columns. """
        return list(range(cls.get_header_depth()))

    @classmethod
    def _format_data(cls,
                     data: pd.DataFrame, *,
                     precision: int | None,
                     squeeze: bool):
        """ Perform final formatting on fetched data. """
        if precision is not None:
            # Round the data to the desired number of decimal places.
            data = data.round(precision)
        if squeeze:
            # If the DataFrame has exactly one column, then return that
            # column as a Series; otherwise, raise an error.
            if data.columns.size != 1:
                raise ValueError("Only data with exactly 1 column can be "
                                 f"squeezed, but got {data.columns}")
            data = data[data.columns[0]]
        return data

    @cached_property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        """ Table's data. """

    def _get_header(self):
        return parse_header(self.data.columns)

    @abstractmethod
    def _fetch_data(self,
                    columns: pd.Index,
                    exclude_masked: bool = False) -> pd.DataFrame:
        """ Fetch data from the table. """

    def fetch_count(self, *,
                    exclude_masked: bool = False,
                    squeeze: bool = False,
                    **kwargs) -> pd.Series | pd.DataFrame:
        """ Fetch counts of one or more columns. """
        return self._format_data(self._fetch_data(self.header.select(**kwargs),
                                                  exclude_masked),
                                 precision=None,
                                 squeeze=squeeze)

    def fetch_ratio(self, *,
                    exclude_masked: bool = False,
                    squeeze: bool = False,
                    precision: int | None = None,
                    quantile: float = 0.,
                    **kwargs) -> pd.Series | pd.DataFrame:
        """ Fetch ratios of one or more columns. """
        # Fetch the data for the numerator.
        numer = self._fetch_data(self.header.select(**kwargs), exclude_masked)
        # Fetch the data for the denominator.
        denom = self._fetch_data(_get_denom_cols(numer.columns), exclude_masked)
        # Compute the ratio of the numerator and the denominator.
        return self._format_data(winsorize(numer / denom.values, quantile),
                                 precision=precision,
                                 squeeze=squeeze)


class PositionTable(RelTypeTable, ABC):
    """ Table indexed by position. """

    MASK = "pos-mask"

    @classmethod
    def get_file_seg_type(cls):
        return path.PositionTableSeg

    @classmethod
    def get_by_read(cls):
        return False

    @classmethod
    def get_index_depth(cls):
        return len(SEQ_INDEX_NAMES)

    @cached_property
    def range(self):
        return self.data.index

    @cached_property
    def range_int(self):
        return index_to_pos(self.range)

    @cached_property
    def seq(self):
        return index_to_seq(self.range)

    @property
    def end5(self):
        return int(self.range_int[0])

    @property
    def end3(self):
        return int(self.range_int[-1])

    @cached_property
    def region(self):
        """ Region covered by the table. """
        # Infer masked positions from the table.
        masked_bool = self.data.isna().all(axis=1)
        unmasked_bool = ~masked_bool
        unmasked = self.data.index[unmasked_bool]
        unmasked_int = index_to_pos(unmasked)
        region = Region(self.ref,
                        self.seq,
                        seq5=self.end5,
                        end5=self.end5,
                        end3=self.end3,
                        name=self.reg)
        region.add_mask(self.MASK, unmasked_int, complement=True)
        return region

    def _fetch_data(self,
                    columns: pd.Index,
                    exclude_masked: bool = False):
        return (self.data.loc[self.region.unmasked, columns] if exclude_masked
                else self.data.loc[:, columns])

    @abstractmethod
    def _iter_profiles(self, *,
                       regions: Iterable[Region] | None,
                       quantile: float,
                       rel: str,
                       k: int | None,
                       clust: int | None) -> Generator[RNAProfile, Any, Any]:
        """ Yield RNA mutational profiles from the table. """

    def iter_profiles(self, *,
                      regions: Iterable[Region] | None = None,
                      quantile: float = 0.,
                      rel: str = MUTAT_REL,
                      k: int | None = None,
                      clust: int | None = None):
        """ Yield RNA mutational profiles from the table. """
        yield from self._iter_profiles(regions=regions,
                                       quantile=quantile,
                                       rel=rel,
                                       k=k,
                                       clust=clust)

    def _compute_ci(self,
                    confidence: float,
                    use_ratio: bool, *,
                    exclude_masked: bool = False,
                    **kwargs):
        """ Calculate the confidence intervals of counts or ratios.

        Parameters
        ----------
        confidence: float
            Confidence level; must be in [0, 1).
        use_ratio: bool
            Compute confidence intervals of the ratio, not the count.
        exclude_masked: bool = False
            Exclude masked positions from the output.
        **kwargs
            Keyword arguments for fetch methods.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Lower and upper bounds of the confidence interval.
        """
        from scipy.stats import binom
        if not 0. <= confidence < 1.:
            raise ValueError(
                f"confidence level must be in [0, 1), but got {confidence}")
        # Fetch the probability of each relationship.
        p = self.fetch_ratio(exclude_masked=True, **kwargs)
        # Fetch the number of reads for each relationship.
        n = np.asarray(np.round(self._fetch_data(_get_denom_cols(p.columns),
                                                 exclude_masked=True).values),
                       dtype=int)
        # Model the counts using binomial distributions.
        lo, up = binom.interval(confidence, n, p.values)
        # Copy the confidence interval bounds into two DataFrames.
        index = self.region.unmasked if exclude_masked else self.data.index
        lower = pd.DataFrame(np.nan, index=index, columns=p.columns)
        upper = pd.DataFrame(np.nan, index=index, columns=p.columns)
        lower.loc[self.region.unmasked] = lo / n if use_ratio else lo
        upper.loc[self.region.unmasked] = up / n if use_ratio else up
        return lower, upper

    def ci_count(self, confidence: float, **kwargs):
        """ Confidence intervals of counts, under these simplifications:

        - Counts are independent of each other.
        - Counts follow binomial distributions.
        - Coverage counts are constant.

        Parameters
        ----------
        confidence: float
            Confidence level; must be in [0, 1).
        **kwargs
            Keyword arguments for fetch methods.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Lower and upper bounds of the confidence interval.
        """
        return self._compute_ci(confidence, False, **kwargs)

    def ci_ratio(self, confidence: float, **kwargs):
        """ Confidence intervals of ratios, under these simplifications:

        - Ratios are independent of each other.
        - Ratios follow beta distributions.
        - Coverage counts are constant.

        Parameters
        ----------
        confidence: float
            Confidence level; must be in [0, 1).
        **kwargs
            Keyword arguments for fetch methods.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Lower and upper bounds of the confidence interval.
        """
        return self._compute_ci(confidence, True, **kwargs)

    def _resample_clust(self,
                        k: int,
                        clust: int,
                        n_cov: np.ndarray, *,
                        seed: int | None = None,
                        tol: float = 1.e-3):
        """ Resample mutations for one cluster.

        Parameters
        ----------
        k: int
            Number of clusters.
        clust: int
            Number of the cluster.
        n_cov: numpy.ndarray
            Number of reads covering each position (1-D array).
        seed: int | None = None
            Seed for the random number generator.
        tol: float = 1.e-3
            Tolerance for floating-point rounding error.
        """
        fetch_kwargs = dict(k=k, clust=clust, exclude_masked=True)
        rng = np.random.default_rng(seed)
        if n_cov.ndim != 1:
            raise ValueError("n_cov must have exactly 1 dimension, "
                             f"but got {n_cov.ndim}")
        resampled = {COVER_REL: n_cov}
        # Resample unambiguous reads.
        p_inf = self.fetch_ratio(rel=INFOR_REL,
                                 squeeze=True,
                                 **fetch_kwargs).values
        n_inf = rng.binomial(n_cov, p_inf)
        resampled[INFOR_REL] = n_inf
        # Resample mutations and matches.
        p_mut = self.fetch_ratio(rel=MUTAT_REL,
                                 squeeze=True,
                                 **fetch_kwargs).values
        n_mut = rng.binomial(n_inf, p_mut)
        resampled[MUTAT_REL] = n_mut
        resampled[MATCH_REL] = n_inf - n_mut

        def resample_multi(n_any_rel: np.ndarray,
                           p_any_rel: np.ndarray,
                           rel_kinds: list[str],
                           tolerance: float):
            """ Resample multiple kinds of relationships at once using a
            multinomial distribution.

            Parameters
            ----------
            n_any_rel: numpy.ndarray
                Number of mutations of any kind at each position.
            p_any_rel: numpy.ndarray
                Probability that each position has any kind of mutation.
            rel_kinds: list[str]
                Kinds of relationships to resample.
            tolerance: float
                Tolerance for floating-point rounding error.
            """
            if n_any_rel.ndim != 1:
                raise ValueError("n_any_rel must have exactly 1 dimension, "
                                 f"but got {n_any_rel.ndim}")
            if p_any_rel.ndim != 1:
                raise ValueError("p_any_rel must have exactly 1 dimension, "
                                 f"but got {p_any_rel.ndim}")
            if tolerance < 0.:
                raise ValueError(f"tolerance must be ≥ 0, but got {tolerance}")
            # Fetch the probability of each given kind of relationship.
            p_given_rel = self.fetch_ratio(rel=rel_kinds, **fetch_kwargs).values
            # Make each probability conditional on a relationship of any
            # kind by dividing by the probability of any relationship.
            with np.errstate(divide="ignore"):
                p_given_rel = np.nan_to_num(p_given_rel
                                            / p_any_rel.reshape(-1, 1))
            # Compute the probability of another type of relationship.
            p_other_rel = 1. - p_given_rel.sum(axis=1)
            # Due to rounding errors in floating-point arithmetic, these
            # probabilities may have small negative values: zero them.
            p_other_rel[np.logical_and(p_other_rel < 0.,
                                       p_other_rel >= -tolerance)] = 0.
            # Join all probabilities into one array.
            p_every_rel = np.hstack([p_given_rel, p_other_rel.reshape(-1, 1)])
            # To compensate, re-normalize the sum of each row to 1.
            p_every_rel /= p_every_rel.sum(axis=1).reshape(-1, 1)
            # Resample each kind of relationship.
            n_every_rel = rng.multinomial(n_any_rel, p_every_rel)
            # Extract the number of each kind of relationship.
            n_each_rel = pd.DataFrame(n_every_rel[:, :len(rel_kinds)],
                                      index=self.region.unmasked,
                                      columns=rel_kinds)
            return n_each_rel

        # Resample each kind of mutation.
        mut_kinds = [SUBST_REL, DELET_REL, INSRT_REL]
        n_mut_kinds = resample_multi(n_mut, p_mut, mut_kinds, tol)
        for kind, n_kind in n_mut_kinds.items():
            resampled[kind] = n_kind.values
        # Resample each kind of substitution.
        sub_kinds = [SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL]
        n_sub = n_mut_kinds[SUBST_REL].values
        p_sub = self.fetch_ratio(rel=SUBST_REL,
                                 squeeze=True,
                                 **fetch_kwargs).values
        n_sub_kinds = resample_multi(n_sub, p_sub, sub_kinds, tol)
        for kind, n_kind in n_sub_kinds.items():
            resampled[kind] = n_kind.values
        # Return the resampled counts.
        return resampled

    def resample(self,
                 fraction: float = 1., *,
                 exclude_masked: bool = False,
                 seed: int | None = None,
                 max_seed: int = 2 ** 32):
        """ Resample the reads and return a new DataFrame.

        Parameters
        ----------
        fraction: float = 1.
            Number of reads to resample, expressed as a fraction of the
            original number of reads. Must be ≥ 0; may be > 1.
        exclude_masked: bool = False
            Exclude positions that have been masked.
        seed: int | None = None
            Seed for the random number generator.
        max_seed: int = 2 ** 32
            Maximum seed to pass to the next random number generator.
        """
        rng = np.random.default_rng(seed)
        # Initialize the resampled data.
        resampled = pd.DataFrame(np.nan,
                                 index=(self.region.unmasked if exclude_masked
                                        else self.data.index),
                                 columns=self.data.columns)
        for k, clust in self.header.clusts:
            # Get the read coverage for the cluster.
            n_cov = np.asarray(np.round(
                fraction * self.fetch_count(rel=COVER_REL,
                                            k=k,
                                            clust=clust,
                                            squeeze=True,
                                            exclude_masked=True).values
            ), dtype=int)
            # Resample the other relationships for the cluster.
            resampled_clust = self._resample_clust(k,
                                                   clust,
                                                   n_cov,
                                                   seed=rng.integers(max_seed))
            # Record the resampled data for the cluster.
            for rel, n_rel in resampled_clust.items():
                key = (rel, k, clust) if k > 0 else rel
                resampled.loc[self.region.unmasked, key] = n_rel
        return resampled


class ReadTable(RelTypeTable, ABC):
    """ Table indexed by read. """

    @classmethod
    def get_file_seg_type(cls):
        return path.ReadTableSeg

    @classmethod
    def get_by_read(cls):
        return True

    @classmethod
    def get_index_depth(cls):
        return len(RB_INDEX_NAMES)

    @property
    def reads(self):
        return self.data.index.values

    def _fetch_data(self,
                    columns: pd.Index,
                    exclude_masked: bool = False):
        return self.data.loc[:, columns]


class AbundanceTable(Table, ABC):
    """ Table of abundances. """

    @classmethod
    def get_file_seg_type(cls):
        return path.AbundanceTableSeg

    @classmethod
    def get_by_read(cls):
        return False

    @classmethod
    def get_index_depth(cls):
        return cls.get_header_depth()

    @cached_property
    @abstractmethod
    def data(self) -> pd.Series:
        """ Table's data. """

    @cached_property
    @abstractmethod
    def proportions(self) -> pd.Series:
        """ Proportion of each item. """
