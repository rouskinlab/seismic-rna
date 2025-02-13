from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from itertools import chain
from typing import Iterable

import numpy as np
import pandas as pd

from .validate import (require_equal,
                       require_array_equal,
                       require_index_equals,
                       require_atleast,
                       require_atmost,
                       require_isinstance)

# Index level names
REL_NAME = "Relationship"
NUM_CLUSTS_NAME = "K"
CLUST_NAME = "Cluster"

# Profile name prefixes
AVERAGE_PREFIX = "average"
CLUSTER_PREFIX = "cluster"

# Header selection keys
K_CLUST_KEY = "k_clust_list"

# No-cluster lists
NO_K = 0
NO_KS = [NO_K]
NO_CLUST = 0, 0
NO_CLUSTS = [NO_CLUST]


def validate_k_clust(k: int, clust: int):
    """ Validate a pair of k and cluster numbers.

    Parameters
    ----------
    k: int
        Number of clusters
    clust: int
        Cluster number

    Returns
    -------
    None
        If the k and cluster numbers form a valid pair.

    Raises
    ------
    TypeError
        If k or clust is not an integer.
    ValueError
        If k and clust do not form a valid pair.
    """
    require_atleast("k", k, 0, classes=int)
    min_clust = 1 if k > 0 else 0
    require_atleast("clust", clust, min_clust, classes=int)
    require_atmost("clust", clust, k, "k", classes=int)


def validate_ks(ks: Iterable):
    """ Validate and sort numbers of clusters.

    Parameters
    ----------
    ks: Iterable
        Numbers of clusters

    Returns
    -------
    list[int]
        Sorted numbers of clusters

    Raises
    ------
    ValueError
        If any k is not positive or is repeated.
    """
    validated = set()
    for k in map(int, ks):
        require_atleast("k", k, 1, classes=int)
        if k in validated:
            raise ValueError(f"Duplicate k: {k}")
        validated.add(k)
    return sorted(validated)


def deduplicate_rels(rels: Iterable):
    """ Remove duplicate relationships while preserving their order.

    Parameters
    ----------
    rels: Iterable
        Relationships

    Returns
    -------
    list[str]
        Relationships with duplicates removed, in the original order.
    """
    ordered = list()
    unordered = set()
    for rel in map(str, rels):
        if rel not in unordered:
            ordered.append(rel)
            unordered.add(rel)
    return ordered


def format_clust_name(k: int, clust: int):
    """ Format a pair of k and cluster numbers into a name.

    Parameters
    ----------
    k: int
        Number of clusters
    clust: int
        Cluster number

    Returns
    -------
    str
        Name specifying k and clust, or "average" if k is 0.
    """
    validate_k_clust(k, clust)
    return f"{CLUSTER_PREFIX} {k}-{clust}" if k > 0 else AVERAGE_PREFIX


def format_clust_names(clusts: Iterable[tuple[int, int]],
                       allow_duplicates: bool = False):
    """ Format pairs of k and clust into a list of names.

    Parameters
    ----------
    clusts: Iterable[tuple[int, int]]
        Zero or more pairs of k and cluster numbers.
    allow_duplicates: bool = False
        Allow k and clust pairs to be duplicated.

    Returns
    -------
    list[str]
        List of names of the pairs of k and clust.

    Raises
    ------
    ValueError
        If `allow_duplicates` is False and clusts has duplicates.
    """
    if not allow_duplicates:
        counts = Counter(clusts := list(clusts))
        if dups := [clust for clust, count in counts.items() if count > 1]:
            raise ValueError(f"Duplicate clusters: {dups}")
    return [format_clust_name(k, clust) for k, clust in clusts]


def list_clusts(k: int):
    """ List all cluster numbers for one k.

    Parameters
    ----------
    k: int
        Number of clusters (≥ 1)

    Returns
    -------
    list[int]
        List of cluster numbers.
    """
    require_atleast("k", k, 1, classes=int)
    return list(range(1, k + 1))


def list_k_clusts(k: int):
    """ List k and cluster numbers as 2-tuples for one k.

    Parameters
    ----------
    k: int
        Number of clusters (≥ 1)

    Returns
    -------
    list[tuple[int, int]]
        List wherein each item is a tuple of the number of clusters
        and the cluster number.
    """
    return [(k, clust) for clust in list_clusts(k)]


def list_ks_clusts(ks: Iterable[int]):
    """ List k and cluster numbers as 2-tuples.

    Parameters
    ----------
    ks: Iterable[int]

    Returns
    -------
    list[tuple[int, int]]
        List wherein each item is a tuple of the number of clusters
        and the cluster number.
    """
    return list(chain(*map(list_k_clusts, validate_ks(ks))))


class Header(ABC):
    """ Header for a table. """

    @classmethod
    @abstractmethod
    def clustered(cls) -> bool:
        """ Whether the header has clusters. """

    @classmethod
    @abstractmethod
    def levels(cls):
        """ Levels of the index. """
        return dict()

    @classmethod
    def num_levels(cls):
        """ Number of levels. """
        return len(cls.levels())

    @classmethod
    def level_keys(cls):
        """ Level keys of the index. """
        return list(cls.levels().keys())

    @classmethod
    def level_names(cls):
        """ Level names of the index. """
        return list(cls.levels().values())

    @property
    @abstractmethod
    def ks(self) -> list[int]:
        """ Numbers of clusters. """

    @cached_property
    @abstractmethod
    def clusts(self) -> list[tuple[int, int]]:
        """ Tracks of data: clusters for clustered data, otherwise one
        track of the average. """

    @cached_property
    def names(self):
        """ Formatted name of each track. """
        return format_clust_names(self.clusts, not self.clustered())

    @cached_property
    def signature(self):
        """ Signature of the header, which will generate an identical
        header if passed as keyword arguments to `make_header`. """
        return dict()

    @cached_property
    @abstractmethod
    def index(self) -> pd.Index:
        """ Index of the header. """

    @abstractmethod
    def iter_clust_indexes(self):
        """ For each cluster, yield an Index/MultiIndex of every column
        that is part of the cluster. """

    @property
    def size(self):
        """ Number of items in the Header. """
        return self.index.size

    def select(self, **kwargs) -> pd.Index:
        """ Select and return items from the header as an Index. """
        index = self.index
        selected = np.ones(index.size, dtype=bool)
        # Handle combinations of k and clust in a list of tuples.
        if value := kwargs.pop(K_CLUST_KEY, None):
            require_isinstance(K_CLUST_KEY, value, list)
            k_name = self.levels().get('k')
            clust_name = self.levels().get('clust')
            combo_selected = np.zeros(index.size, dtype=bool)
            for k, clust in value:
                # Find rows that match each specified combination.
                k_match = index.get_level_values(k_name) == k
                clust_match = index.get_level_values(clust_name) == clust
                combo_selected |= (k_match & clust_match)
            if np.all(selected):
                selected &= combo_selected
            else:
                selected |= combo_selected
        # Handle relationship selection.
        for key, name in self.levels().items():
            if value := kwargs.pop(key, None):
                level_values = index.get_level_values(name)
                equal_values = np.isin(level_values, np.atleast_1d(value))
                if not np.any(equal_values):
                    expect = np.unique(level_values).tolist()
                    raise ValueError(f"{key} must be in {repr(expect)}, "
                                     f"but got {repr(value)}")
                selected &= equal_values
        # Check if any extra keyword arguments were given; allow extra
        # arguments only if their values are falsy (e.g. None, 0).
        if extras := {k: v for k, v in kwargs.items() if v}:
            raise TypeError(
                f"Extra keyword arguments for {type(self)}: {extras}"
            )
        return index[selected]

    def modified(self, **kwargs):
        """ Return a new header with a possibly modified signature.

        Parameters
        ----------
        **kwargs
            Keyword arguments for modifying the signature of the header.
            Each argument given here will be passed to `make_header` and
            override the attribute (if any) with the same name in this
            header's signature. Attributes of this header's signature
            that are not overriden will also be passed to `make_header`.

        Returns
        -------
        Header
            New header with a possibly modified signature.
        """
        return make_header(**(self.signature | kwargs))

    def get_clust_header(self):
        """ Corresponding ClustHeader. """
        return self.modified(rels=None)

    def get_rel_header(self):
        """ Corresponding RelHeader. """
        return self.modified(ks=None)

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, Header):
            return NotImplemented
        if type(other) is not type(self):
            return False
        for key, value in self.signature.items():
            other_value = other.signature[key]
            require_isinstance(key, other_value, type(value))
            if value != other_value:
                return False
        return True

    def __str__(self):
        return f"{type(self).__name__}({self.signature})"

    def __repr__(self):
        return str(self)


class RelHeader(Header):
    """ Header of relationships. """

    @classmethod
    def clustered(cls):
        return False

    @classmethod
    def levels(cls):
        return super().levels() | dict(rel=REL_NAME)

    def __init__(self, *, rels: Iterable[str], **kwargs):
        """
        Parameters
        ----------
        rels: Iterable[str]
            One or more relationships in the header.
        """
        super().__init__(**kwargs)
        self._rels = deduplicate_rels(rels)

    @property
    def rels(self):
        """ Relationships. """
        return self._rels

    @property
    def ks(self):
        return NO_KS

    @property
    def clusts(self):
        return NO_CLUSTS

    @cached_property
    def signature(self):
        return super().signature | dict(rels=self.rels)

    @cached_property
    def index(self):
        return pd.Index(self.rels, name=REL_NAME)

    def iter_clust_indexes(self):
        yield self.index


class ClustHeader(Header):
    """ Header of clusters. """

    @classmethod
    def clustered(cls):
        return True

    @classmethod
    def levels(cls):
        return super().levels() | dict(k=NUM_CLUSTS_NAME, clust=CLUST_NAME)

    def __init__(self, *, ks: Iterable[int], **kwargs):
        super().__init__(**kwargs)
        self._ks = validate_ks(ks)

    @property
    def ks(self):
        return self._ks

    @cached_property
    def signature(self):
        return super().signature | dict(ks=self.ks)

    @cached_property
    def clusts(self):
        return list_ks_clusts(self.ks)

    @cached_property
    def index(self):
        return pd.MultiIndex.from_tuples(self.clusts, names=self.level_names())

    def iter_clust_indexes(self):
        for k, clust in self.clusts:
            yield self.select(k=k, clust=clust)


class RelClustHeader(ClustHeader, RelHeader):
    """ Header of relationships and clusters. """

    @cached_property
    def index(self):
        return pd.MultiIndex.from_tuples([(rel, k, clust)
                                          for rel in self.rels
                                          for k, clust in self.clusts],
                                         names=self.level_names())


def make_header(*,
                rels: Iterable[str] | None = None,
                ks: Iterable[int] | None = None):
    """ Make a new Header of an appropriate type.

    Parameters
    ----------
    rels: Iterable[str] | None = None
        Relationships in the header
    ks: Iterable[int] | None = None
        Numbers of clusters

    Returns
    -------
    Header
        Header of the appropriate type.
    """
    if ks is not None and ks != NO_KS:
        if rels is not None:
            return RelClustHeader(rels=rels, ks=ks)
        return ClustHeader(ks=ks)
    if rels is not None:
        return RelHeader(rels=rels)
    raise TypeError("Must give rels, ks, or both, but got neither")


def parse_header(index: pd.Index | pd.MultiIndex):
    """ Parse an Index into a Header of an appropriate type.

    Parameters
    ----------
    index: pd.Index | pd.MultiIndex
        Index to parse.

    Returns
    -------
    Header
        New Header whose index is `index`.
    """
    if isinstance(index, pd.MultiIndex):
        names = index.names
        if REL_NAME in names:
            rels = list(map(str,
                            deduplicate_rels(index.get_level_values(REL_NAME))))
        else:
            rels = None
        if NUM_CLUSTS_NAME in names:
            ks = list(map(int,
                          np.unique(index.get_level_values(NUM_CLUSTS_NAME))))
        else:
            ks = None
    elif index.name is None or index.name == REL_NAME:
        rels = deduplicate_rels(index.values)
        ks = None
    else:
        raise ValueError(f"index must be named {repr(REL_NAME)}, "
                         f"but got {repr(index.name)}")
    header = make_header(rels=rels, ks=ks)
    # Verify that the index of the new header matches the given index.
    if isinstance(index, pd.MultiIndex):
        require_equal("index names",
                      sorted(index.names),
                      sorted(header.level_names()))
        for name in header.level_names():
            require_array_equal(
                f"values of index level {repr(name)}",
                np.asarray(index.get_level_values(name).values,
                           dtype=str if name == REL_NAME else int),
                header.index.get_level_values(name).values,
                "values of header index"
            )
    else:
        require_index_equals("index", index, header.index, "header index")
    return header
