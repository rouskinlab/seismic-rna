from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property, partial
from itertools import chain, product
from typing import Iterable

import numpy as np
import pandas as pd

# Index level names
REL_NAME = "Relationship"
ORDER_NAME = "Order"
CLUST_NAME = "Cluster"

# Profile name prefixes
AVERAGE_PREFIX = "average"
CLUSTER_PREFIX = "cluster"


def validate_order_clust(order: int, clust: int, allow_zero: bool = False):
    """ Validate a pair of order and cluster numbers.

    Parameters
    ----------
    order: int
        Order number
    clust: int
        Cluster number
    allow_zero: bool = False
        Allow order and cluster to be 0.

    Returns
    -------
    None
        If the order and cluster numbers form a valid pair.

    Raises
    ------
    TypeError
        If order or cluster is not an integer.
    ValueError
        If the order and cluster numbers do not form a valid pair.
    """
    if not isinstance(order, int):
        raise TypeError(f"order must be an int, but got {type(order).__name__}")
    if not isinstance(clust, int):
        raise TypeError(f"clust must be an int, but got {type(clust).__name__}")
    min_order = int(not allow_zero)
    if order < min_order:
        raise ValueError(f"order must be ≥ {min_order}, but got {order}")
    min_clust = int(order > 0)
    if clust < min_clust:
        raise ValueError(f"clust must be ≥ {min_clust}, but got {clust}")
    if clust > order:
        raise ValueError(f"clust must be ≤ order, but got {clust} > {order}")


def format_clust_name(order: int, clust: int, allow_zero: bool = False):
    """ Format a pair of order and cluster numbers into a name.

    Parameters
    ----------
    order: int
        Order number
    clust: int
        Cluster number
    allow_zero: bool = False
        Allow order and cluster to be 0.

    Returns
    -------
    str
        Name specifying the order and cluster numbers, or "average" if
        the order number is 0.
    """
    validate_order_clust(order, clust, allow_zero)
    return f"{CLUSTER_PREFIX} {order}-{clust}" if order > 0 else AVERAGE_PREFIX


def format_clust_names(clusts: Iterable[tuple[int, int]],
                       allow_zero: bool = False,
                       allow_duplicates: bool = False):
    """ Format pairs of order and cluster numbers into a list of names.

    Parameters
    ----------
    clusts: Iterable[tuple[int, int]]
        Zero or more pairs of order and cluster numbers.
    allow_zero: bool = False
        Allow order and cluster to be 0.
    allow_duplicates: bool = False
        Allow order and cluster pairs to be duplicated.

    Returns
    -------
    list[str]
        List of names of the pairs of order and cluster numbers.

    Raises
    ------
    ValueError
        If `allow_duplicates` is False and an order-cluster pair occurs
        more than once.
    """
    if not allow_duplicates:
        counts = Counter(clusts := list(clusts))
        if dups := [clust for clust, count in counts.items() if count > 1]:
            raise ValueError(f"Duplicate clusters: {dups}")
    return [format_clust_name(order, clust, allow_zero=allow_zero)
            for order, clust in clusts]


def list_clusts(order: int):
    """ List all cluster numbers for one order.

    Parameters
    ----------
    order: int
        Number of clusters (≥ 0)

    Returns
    -------
    list[int]
        List of cluster numbers.
    """
    if order < 0:
        raise ValueError(f"order must be ≥ 0, but got {order}")
    return list(range(1, order + 1)) if order > 0 else [0]


def list_orders(max_order: int, min_order: int = 1):
    """ List order numbers from `min_order` to `max_order`.

    Parameters
    ----------
    max_order: int
        Maximum number of clusters (≥ 0)
    min_order: int = 1
        Minimum number of clusters (≥ 1)

    Returns
    -------
    list[int]
        List of numbers of clusters
    """
    if min_order < 1:
        raise ValueError(f"min_order must be ≥ 1, but got {min_order}")
    if max_order < 0:
        raise ValueError(f"max_order must be ≥ 0, but got {max_order}")
    if max_order == 0:
        if min_order != 1:
            raise ValueError("If max_order = 0, then min_order must be 1, "
                             f"but got {min_order}")
        return [0]
    if max_order < min_order:
        raise ValueError("If not 0, max_order must be ≥ min_order, "
                         f"but got {max_order} < {min_order}")
    return list(range(min_order, max_order + 1))


def index_clusts(order: int):
    """ Index of cluster numbers for one order.

    Parameters
    ----------
    order: int
        Number of clusters (≥ 0)

    Returns
    -------
    pd.Index
        Index of cluster numbers
    """
    return pd.Index(list_clusts(order), name=CLUST_NAME)


def index_orders(max_order: int, min_order: int = 1):
    """ Index of order numbers from `min_order` to `max_order`.

    Parameters
    ----------
    max_order: int
        Maximum number of clusters (≥ 0)
    min_order: int = 1
        Minimum number of clusters (≥ 1)

    Returns
    -------
    pd.Index
        Index of order numbers
    """
    return pd.Index(list_orders(max_order, min_order), name=ORDER_NAME)


def list_order_clusts(order: int):
    """ List order and cluster numbers as 2-tuples for one order.

    Parameters
    ----------
    order: int
        Number of clusters (≥ 0)

    Returns
    -------
    list[tuple[int, int]]
        List wherein each item is a tuple of the order of clustering
        (i.e. number of clusters) and the cluster number.
    """
    return list(product([order], list_clusts(order)))


def list_orders_clusts(max_order: int, min_order: int = 1):
    """ List order and cluster numbers as 2-tuples for every order from
    `min_order` to `max_order`.

    Parameters
    ----------
    max_order: int
        Maximum number of clusters (≥ 0)
    min_order: int = 1
        Minimum number of clusters (≥ 1)

    Returns
    -------
    list[tuple[int, int]]
        List wherein each item is a tuple of the order of clustering
        (i.e. number of clusters) and the cluster number.
    """
    return list(chain(*map(list_order_clusts,
                           list_orders(max_order, min_order))))


def index_order_clusts(order: int):
    """ List order and cluster numbers as a MultiIndex for one order.

    Parameters
    ----------
    order: int
        Number of clusters (≥ 0)

    Returns
    -------
    pd.MultiIndex
        Index wherein each item is a tuple of the order of clustering
        (i.e. number of clusters) and the cluster number.
    """
    return pd.MultiIndex.from_tuples(list_order_clusts(order),
                                     names=ClustHeader.level_names())


def index_orders_clusts(max_order: int, min_order: int = 1):
    """ List order and cluster numbers as a MultiIndex for every order
    from `min_order` to `max_order`.

    Parameters
    ----------
    max_order: int
        Maximum number of clusters (≥ 0)
    min_order: int = 1
        Minimum number of clusters (≥ 1)

    Returns
    -------
    pd.MultiIndex
        Index wherein each item is a tuple of the order of clustering
        (i.e. number of clusters) and the cluster number.
    """
    return pd.MultiIndex.from_tuples(list_orders_clusts(max_order, min_order),
                                     names=ClustHeader.level_names())


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
    def max_order(self) -> int:
        """ Maximum number of clusters (≥ 1) if clustered, else 0. """

    @property
    @abstractmethod
    def min_order(self) -> int:
        """ Minimum number of clusters (≥ 1) if clustered, else 1. """

    @cached_property
    def orders(self):
        """ Index of order numbers. """
        return index_orders(self.max_order, self.min_order)

    @cached_property
    def clusts(self):
        """ Order and cluster numbers of the header. """
        return index_orders_clusts(self.max_order, self.min_order)

    @cached_property
    def names(self):
        """ Formatted name of each cluster. """
        return format_clust_names(self.clusts, not self.clustered())

    @cached_property
    def signature(self):
        """ Signature of the header, which will generate an identical
        header if passed as keyword arguments to `make_header`. """
        return dict(max_order=self.max_order, min_order=self.min_order)

    @cached_property
    @abstractmethod
    def index(self) -> pd.Index:
        """ Index of the header. """

    def iter_clust_indexes(self):
        """ For each cluster in the header, yield an Index or MultiIndex
        of every column in the header that is part of the cluster. """
        if self.max_order == 0:
            # There are no clusters, so just yield all relationships.
            yield self.index
        else:
            # For each cluster in the header, select the corresponding
            # order and cluster in the index.
            for order, clust in self.clusts:
                yield self.select(order=order, clust=clust)

    @property
    def size(self):
        """ Number of items in the Header. """
        return self.index.size

    def select(self, **kwargs) -> pd.Index:
        """ Select and return items from the header as an Index. """
        index = self.index
        selected = np.ones(index.size, dtype=bool)
        for key, name in self.levels().items():
            if value := kwargs.pop(key, None):
                level_values = index.get_level_values(name)
                equal_values = np.isin(level_values, np.atleast_1d(value))
                if not np.any(equal_values):
                    expect = np.unique(level_values).tolist()
                    raise ValueError(f"Expected {key} to be one of {expect}, "
                                     f"but got {repr(value)}")
                selected &= equal_values
        # Check if any extra keyword arguments were given; allow extra
        # arguments only if their values are falsy (e.g. None, 0).
        if extras := {k: v for k, v in kwargs.items() if v}:
            raise TypeError("Unexpected keyword arguments for "
                            f"{type(self).__name__}: {extras}")
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

    @property
    @abstractmethod
    def _descripts(self):
        """ Description keys. """
        return dict()

    @cached_property
    def _descripts_str(self):
        """ String of description keys. """
        return ", ".join(f"{k}={v}" for k, v in self._descripts.items())

    def __eq__(self, other):
        if not isinstance(other, Header):
            return NotImplemented
        if not isinstance(other, type(self)):
            return False
        if self is other:
            return True
        for key, value in self.signature.items():
            other_value = other.signature[key]
            if not isinstance(other_value, type(value)):
                raise TypeError(f"For key {repr(key)}, types of {repr(value)} "
                                f"and {repr(other_value)} differ")
            if isinstance(value, int):
                if value != other_value:
                    return False
            elif isinstance(value, np.ndarray):
                if not np.array_equal(value, other_value):
                    return False
            elif isinstance(value, pd.Index):
                if not value.equals(other_value):
                    return False
            else:
                raise TypeError(f"Invalid type for {repr(key)}: {repr(value)}")
        return True

    def __str__(self):
        return f"{type(self).__name__}({self._descripts_str})"

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
        # Convert rels to an array of strings.
        rels = np.asarray(rels, dtype=str)
        if rels.size == 0:
            raise ValueError("Got no relationships for header")
        # Find the unique relationships.
        _, uniq_indexes = np.unique(rels, return_index=True)
        # Get the unique relationships in their original order.
        self._rels = rels[np.sort(uniq_indexes)]

    @property
    def rels(self) -> np.ndarray:
        """ Relationships. """
        return self._rels

    @property
    def max_order(self):
        return 0

    @property
    def min_order(self):
        return 1

    @cached_property
    def signature(self):
        return super().signature | dict(rels=self.rels)

    @cached_property
    def index(self):
        return pd.Index(self.rels, name=REL_NAME)

    @property
    def _descripts(self):
        return super()._descripts | dict(rels=self.rels.tolist())


class ClustHeader(Header):
    """ Header of order and cluster numbers. """

    @classmethod
    def clustered(cls):
        return True

    @classmethod
    def levels(cls):
        return super().levels() | dict(order=ORDER_NAME, clust=CLUST_NAME)

    def __init__(self, *, max_order: int, min_order: int = 1, **kwargs):
        super().__init__(**kwargs)
        if max_order < 1:
            raise ValueError(f"max_order must be ≥ 1, but got {max_order}")
        self._max_order = max_order
        self._min_order = min_order

    @property
    def max_order(self):
        return self._max_order

    @property
    def min_order(self):
        return self._min_order

    @cached_property
    def index(self):
        return self.clusts

    @property
    def _descripts(self):
        orders = dict(max_order=self.max_order)
        if self.min_order is not None:
            orders.update(min_order=self.min_order)
        return super()._descripts | orders


class RelClustHeader(ClustHeader, RelHeader):
    """ Header of relationships with order and cluster numbers. """

    @cached_property
    def index(self):
        return pd.MultiIndex.from_arrays(
            [np.repeat(self.rels, self.clusts.size),
             *map(partial(np.tile, reps=self.rels.size),
                  map(self.clusts.get_level_values,
                      ClustHeader.level_names()))],
            names=self.level_names()
        )


def make_header(*,
                rels: Iterable[str] = (),
                max_order: int = 0,
                min_order: int = 1):
    """ Make a new Header of an appropriate type.

    Parameters
    ----------
    rels: Iterable[str]
        Relationships in the header.
    max_order: int = 0
        Maximum number of clusters (≥ 1), or 0 if not clustered.
    min_order: int = 1
        Minimum number of clusters (≥ 1), or 1 if not clustered.

    Returns
    -------
    Header
        Header of the appropriate type.
    """
    if rels := list(rels):
        if max_order > 0:
            return RelClustHeader(rels=rels,
                                  max_order=max_order,
                                  min_order=min_order)
        return RelHeader(rels=rels)
    if max_order > 0:
        return ClustHeader(max_order=max_order, min_order=min_order)
    raise TypeError(f"No header for rels={rels} and max_order={max_order}")


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
    kwargs = dict()
    try:
        if isinstance(index, pd.MultiIndex):
            rels = index.get_level_values(REL_NAME).values
        elif index.name is None or index.name == REL_NAME:
            rels = index.values
        else:
            raise ValueError(f"Expected index named {repr(REL_NAME)}, "
                             f"but got {repr(index.name)}")
        kwargs.update(rels=rels)
    except KeyError:
        pass
    try:
        orders = np.asarray(index.get_level_values(ORDER_NAME).values,
                            dtype=int)
        kwargs.update(max_order=np.max(orders), min_order=np.min(orders))
    except (KeyError, ValueError):
        pass
    header = make_header(**kwargs)
    # Verify that the index of the new header matches the given index.
    if isinstance(index, pd.MultiIndex):
        if set(index.names) != set(header.level_names()):
            raise ValueError(f"Expected index names {header.level_names()}, "
                             f"but got {index.names}")
        for name in header.level_names():
            if not np.array_equal(
                    np.asarray(index.get_level_values(name).values,
                               dtype=str if name == REL_NAME else int),
                    header.index.get_level_values(name).values
            ):
                raise ValueError(f"Invalid index level {repr(name)}: expected "
                                 f"{header.index.get_level_values(name)}, "
                                 f"but got {index.get_level_values(name)}")
    else:
        if not index.equals(header.index):
            raise ValueError(f"Indexes do not match: {header.index} ≠ {index}")
    return header

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
