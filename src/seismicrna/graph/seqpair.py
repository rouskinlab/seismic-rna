from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain, product
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .base import CartesianGraph, TwoTableSeqGraph, PRECISION
from .color import SeqColorMap
from .seq import get_table_params
from .write import TwoTableGraphWriter
from ..cluster.names import (CLS_NAME, ORD_NAME, ORD_CLS_NAME,
                             validate_order_cluster)
from ..core import docdef
from ..core.cli import (arg_input_file, opt_rels,
                        opt_y_ratio, opt_quantile, opt_arrange,
                        opt_csv, opt_html, opt_pdf,
                        opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.sect import POS_NAME
from ..table.base import REL_CODES
from ..table.load import (PosTableLoader, RelPosTableLoader, MaskPosTableLoader,
                          ClustPosTableLoader, find_tables, get_clusters)

logger = getLogger(__name__)

# Sample index level name
SAMPLE_NAME = "Sample"
ROW_NAME = "Row"
COL_NAME = "Column"


# Base Sequence Pair Graph #############################################

class SeqPairGraph(CartesianGraph, TwoTableSeqGraph, ABC):
    """ Graph that compares a pair of samples of the same sequence. """

    def __init__(self, *args,
                 rel: str, y_ratio: bool, quantile: float,
                 order1: int = 0, cluster1: int = 0,
                 order2: int = 0, cluster2: int = 0,
                 **kwargs):
        super().__init__(*args, **kwargs)
        if rel not in REL_CODES:
            raise ValueError(f"rel must be one of {list(REL_CODES)}, "
                             f"but got '{rel}'")
        self.rel_code = rel
        self.y_ratio = y_ratio
        self.quantile = quantile
        validate_order_cluster(order1, cluster1, allow_zero=True)
        self._order1 = [order1] if order1 > 0 else None
        self._cluster1 = [cluster1] if cluster1 > 0 else None
        validate_order_cluster(order2, cluster2, allow_zero=True)
        self._order2 = [order2] if order2 > 0 else None
        self._cluster2 = [cluster2] if cluster2 > 0 else None

    @property
    def rel(self):
        """ Return the relationship to graph. """
        return REL_CODES[self.rel_code]

    @property
    def rels(self):
        return [self.rel]

    @classmethod
    def sources(cls) -> dict[type[PosTableLoader], str]:
        """ Names of the sources of data. """
        return {RelPosTableLoader: "Related",
                MaskPosTableLoader: "Masked",
                ClustPosTableLoader: "Clustered"}

    def _get_source(self, table: PosTableLoader):
        table_type = type(table)
        try:
            return self.sources()[table_type]
        except KeyError:
            raise TypeError(f"Invalid table for {self}: {table_type.__name__}")

    @property
    def source1(self):
        """ Source of the first data set. """
        return self._get_source(self.table1)

    @property
    def source2(self):
        """ Source of the second data set. """
        return self._get_source(self.table2)

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @property
    def title(self):
        return (f"{self.quantity} of {self.rel} bases "
                f"in {self.source1} reads from {self.sample1} "
                f"vs. {self.source2} reads from {self.sample2} "
                f"per position in {self.ref}:{self.sect}")

    @property
    def subject1(self):
        if isinstance(self.table1, ClustPosTableLoader):
            osym = 'x' if self._order1 is None else self._order1[0]
            ksym = 'x' if self._cluster1 is None else self._cluster1[0]
            return f"{self.source1}-{osym}-{ksym}"
        return self.source1

    @property
    def subject2(self):
        if isinstance(self.table2, ClustPosTableLoader):
            osym = 'x' if self._order2 is None else self._order2[0]
            ksym = 'x' if self._cluster2 is None else self._cluster2[0]
            return f"{self.source2}-{osym}-{ksym}"
        return self.source2

    @property
    def subject(self):
        return f"{self.subject1}__and__{self.subject2}"

    @property
    def quantity(self):
        return "Ratio" if self.y_ratio else "Count"

    @property
    def predicate(self):
        return f"{self.rel_code}-{self.quantity}"

    @classmethod
    @abstractmethod
    def graph_type(cls):
        """ Type of the graph. """

    @property
    def graph_filename(self):
        return f"{self.subject}_{self.predicate}_{self.graph_type()}".lower()

    def _select_table(self, table: PosTableLoader,
                      order: int | None = None,
                      cluster: int | None = None):
        # Determine which relationships, orders, and clusters to use.
        select = dict(rels=self.rels)
        if order is not None:
            if isinstance(table, ClustPosTableLoader):
                select.update(dict(order=order))
                if cluster is not None:
                    select.update(dict(cluster=cluster))
            else:
                logger.warning(f"{self} cannot use order {order} with {table}")
        elif cluster is not None:
            logger.warning(f"{self} cannot use cluster {cluster} with no order")
        return table.select(self.y_ratio, self.quantile, PRECISION, **select)

    @cached_property
    def data1(self):
        return self._select_table(self.table1, self._order1, self._cluster1)

    @cached_property
    def data2(self):
        return self._select_table(self.table2, self._order2, self._cluster2)

    @property
    def nrows(self):
        return self.data2.columns.size

    @property
    def ncols(self):
        return self.data1.columns.size

    @cached_property
    def clusters1(self):
        return get_clusters(self.data1.columns, allow_zero=True)

    @cached_property
    def clusters2(self):
        return get_clusters(self.data2.columns, allow_zero=True)


class SeqPairOneAxisGraph(SeqPairGraph, ABC):
    """ Graph of a pair of datasets over the same sequence in which the
    data series are merged in some fashion into another series, and the
    original data are not graphed directly. """

    @property
    def x_title(self) -> str:
        return POS_NAME

    @cached_property
    def col_titles(self):
        return get_titles(sample=self.sample1, source=self.source1,
                          clusters=self.clusters1)

    @cached_property
    def row_titles(self):
        return get_titles(sample=self.sample2, source=self.source2,
                          clusters=self.clusters2)

    @classmethod
    @abstractmethod
    def _trace_function(cls) -> Callable:
        """ Function to generate the graph's traces. """

    @property
    @abstractmethod
    def _merge_data(self) -> Callable:
        """ Function to merge the two datasets into one. """

    def _iter_prod(self):
        """ Yield the product of the columns of data tables 1 and 2
        with titles for rows and columns. """
        for (col, (_, s1)), (row, (_, s2)) in product(zip(self.col_titles,
                                                          self.data1.items(),
                                                          strict=True),
                                                      zip(self.row_titles,
                                                          self.data2.items(),
                                                          strict=True)):
            yield (row, col), (s1, s2)

    @cached_property
    def data(self):
        # Take the Cartesian product of the columns in datasets 1 and 2,
        # and merge each pair using self._merge_data().
        data = pd.DataFrame.from_dict({titles: self._merge_data(vals1, vals2)
                                       for titles, (vals1, vals2)
                                       in self._iter_prod()})
        # Indicate that the column index levels now correspond to the
        # rows and columns of the graph.
        data.columns.rename([ROW_NAME, COL_NAME], inplace=True)
        return data

    def get_traces(self):
        for (i, row), (j, col) in product(enumerate(self.row_titles, start=1),
                                          enumerate(self.col_titles, start=1)):
            for trace in self._trace_function()(self.data.loc[:, (row, col)],
                                                self.cmap):
                yield (i, j), trace


class SeqPairTwoAxisGraph(SeqPairGraph, ABC):
    """ Graph of a pair of datasets over the same sequence in which one
    dataset is graphed on each axis. """

    @cached_property
    def data(self):
        # Join data tables 1 and 2 horizontally.
        data = pd.concat([self.data1, self.data2], axis=1, join="inner")

        def make_columns(columns: pd.Index | pd.MultiIndex, sample: str):
            # Sample names.
            samples = [sample] * columns.size
            # Clusters (if any).
            try:
                orders = columns.get_level_values(ORD_NAME).to_list()
                clusters = columns.get_level_values(CLS_NAME).to_list()
            except KeyError:
                orders = [0] * columns.size
                clusters = [0] * columns.size
            cols = {SAMPLE_NAME: samples, ORD_NAME: orders, CLS_NAME: clusters}
            return pd.MultiIndex.from_arrays(list(cols.values()),
                                             names=list(cols.keys()))

        # Add the sample names and order/cluster numbers to the columns
        # for tables 1 and 2.
        cols1 = make_columns(self.data1.columns, self.sample1)
        cols2 = make_columns(self.data2.columns, self.sample2)
        # Drop the order/cluster numbers if they are all zero.
        names = [SAMPLE_NAME] + [name for name in ORD_CLS_NAME
                                 if ((cols1.get_level_values(name) != 0).any()
                                     and
                                     (cols2.get_level_values(name) != 0).any())]
        # Make the columns for the merged table.
        data.columns = pd.MultiIndex.from_arrays(
            [(cols1.get_level_values(name).to_list()
              + cols2.get_level_values(name).to_list())
             for name in names],
            names=names
        )
        return data


# Sequence Pair Graph Writer ###########################################

class SeqPairGraphWriter(TwoTableGraphWriter, ABC):

    @property
    @abstractmethod
    def graph_type(self) -> type[SeqPairGraph]:
        """ Type of the graph to iterate. """

    def iter(self, rels_sets: tuple[str, ...],
             arrange: str, y_ratio: bool, quantile: float):
        _, _, csparams1 = get_table_params(self.table1, arrange)
        _, _, csparams2 = get_table_params(self.table2, arrange)
        for cparams1, cparams2 in product(csparams1, csparams2):
            cparams = dict()
            if (param := cparams1.get("order")) is not None:
                cparams["order1"] = param
            if (param := cparams1.get("cluster")) is not None:
                cparams["cluster1"] = param
            if (param := cparams2.get("order")) is not None:
                cparams["order2"] = param
            if (param := cparams2.get("cluster")) is not None:
                cparams["cluster2"] = param
            for rels in rels_sets:
                yield self.graph_type(table1=self.table1, table2=self.table2,
                                      rel=rels, y_ratio=y_ratio,
                                      quantile=quantile, **cparams)


# Helper functions #####################################################

def get_titles(*, sample: str = "", source: str,
               clusters: Iterable[tuple[int, int]]):
    """ Return a list of column titles given the name of a sample and
    its source, as well as one or more cluster numbers. """
    sample_source = f"{sample}, {source}" if sample else source
    return [f"{sample_source} {order}-{cluster}" if order > 0 else sample_source
            for order, cluster in clusters]


class SeqPairGraphRunner(object):
    params = [
        arg_input_file,
        opt_rels,
        opt_y_ratio,
        opt_quantile,
        opt_arrange,
        opt_csv,
        opt_html,
        opt_pdf,
        opt_max_procs,
        opt_parallel,
    ]

    @classmethod
    @abstractmethod
    def writer_type(cls) -> type[SeqPairGraphWriter]:
        """ Type of SeqPairGraphWriter. """

    @classmethod
    @docdef.auto()
    def run(cls,
            input_file: tuple[str, ...],
            rels: tuple[str, ...], *,
            y_ratio: bool,
            quantile: float,
            arrange: str,
            csv: bool,
            html: bool,
            pdf: bool,
            max_procs: int,
            parallel: bool) -> list[Path]:
        """ Run the graph seqdiff module. """
        tables = list(find_tables(input_file))
        if len(tables) % 2 != 0:
            raise ValueError(
                f"Number of files must be even, but got {len(tables)}")
        writers = [cls.writer_type()(table1_file=t1, table2_file=t2)
                   for t1, t2 in zip(tables[0::2], tables[1::2], strict=True)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs, parallel, pass_n_procs=False,
                                    kwargs=dict(rels_sets=rels,
                                                y_ratio=y_ratio,
                                                quantile=quantile,
                                                arrange=arrange,
                                                csv=csv, html=html, pdf=pdf))))
