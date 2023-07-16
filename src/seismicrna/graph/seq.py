import os
from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import (PRECISION, find_tables, GraphWriter, CartesianGraph,
                   OneTableSeqGraph, OneSampGraph)
from .color import ColorMap, RelColorMap, SeqColorMap
from ..cluster.names import (ENSEMBLE_NAME, CLS_NAME, ORD_NAME, ORD_CLS_NAME,
                             fmt_clust_name, validate_order_cluster)
from ..core import docdef
from ..core.cli import (arg_input_file, opt_rels,
                        opt_stack, opt_arrange, opt_y_ratio,
                        opt_csv, opt_html, opt_pdf, opt_max_procs, opt_parallel,
                        CLUST_INDIV, CLUST_UNITE, CLUST_ORDER)
from ..core.parallel import dispatch
from ..core.sect import BASE_NAME, POS_NAME
from ..core.seq import BASES
from ..table.base import REL_CODES, REL_NAME
from ..table.load import (RelTypeTableLoader, RelPosTableLoader,
                          MaskPosTableLoader, ClustPosTableLoader)

logger = getLogger(__name__)

# Number of digits to which to round decimals.

params = [
    arg_input_file,
    opt_rels,
    opt_stack,
    opt_arrange,
    opt_y_ratio,
    opt_csv,
    opt_html,
    opt_pdf,
    opt_max_procs,
    opt_parallel,
]


@command(__name__.split(os.path.extsep)[-1], params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positional attributes. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_file: tuple[str, ...],
        rels: str, *,
        stack: bool,
        arrange: str,
        y_ratio: bool,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    writers = list(map(SeqGraphWriter, find_tables(input_file)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels=rels, stack=stack,
                                            arrange=arrange, y_ratio=y_ratio,
                                            csv=csv, html=html, pdf=pdf))))


# Helper Functions #####################################################

def get_base_trace(data: pd.Series, cmap: ColorMap, base_int: int):
    # Validate the base.
    base = chr(base_int)
    if base not in BASES.decode():
        raise ValueError(f"Invalid DNA base: '{base}'")
    # Find the position of every base of that type.
    seq_mask = data.index.get_level_values(BASE_NAME) == base
    # Get the values at those positions, excluding NaN values.
    vals = data.loc[seq_mask].dropna()
    # Set the index of the values to the numerical positions.
    vals.index = vals.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x}: {y}" for x, y in vals.items()]
    # Create a trace comprising all bars for this base type.
    return go.Bar(name=base, x=vals.index, y=vals,
                  marker_color=cmap[base_int],
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_base_traces(data: pd.Series, cmap: ColorMap):
    for base in BASES:
        yield get_base_trace(data, cmap, base)


def get_stack_trace(data: pd.Series, cmap: ColorMap):
    # Get the relationship from the name of the data series.
    rel = data.name
    # Get the sequence and positions.
    bases = data.index.get_level_values(BASE_NAME)
    pos = data.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x} {rel}: {y}"
                 for base, x, y in zip(bases, pos, data, strict=True)]
    # Create a trace comprising all bars for this field.
    return go.Bar(name=rel, x=pos, y=data,
                  marker_color=cmap[rel],
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_stack_traces(data: pd.DataFrame, cmap: ColorMap):
    for rel, series in data.items():
        yield get_stack_trace(series, cmap)


# Sequence Graph Writer ################################################

class SeqGraphWriter(GraphWriter):

    def iter(self, rels: str, arrange: str, stack: bool, y_ratio: bool):
        if type(self.table) in EnsembleSeqGraph.sources():
            if stack:
                yield EnsembleMultiRelSeqGraph(table=self.table,
                                               rels=rels,
                                               y_ratio=y_ratio)
            else:
                for rel in rels:
                    yield EnsembleSingleRelSeqGraph(table=self.table,
                                                    rels=rel,
                                                    y_ratio=y_ratio)
        elif type(self.table) in ClusterSeqGraph.sources():
            if arrange == CLUST_INDIV:
                # One file per cluster, with no subplots.
                clusters_params = [dict(order=order, cluster=cluster)
                                   for order, cluster in self.table.ord_clust]
            elif arrange == CLUST_ORDER:
                # One file per order, with a subplot for each cluster.
                orders = sorted(self.table.orders)
                clusters_params = [dict(order=order) for order in orders]
            elif arrange == CLUST_UNITE:
                # One file, with subplots of all clusters of all orders.
                clusters_params = [dict()]
            else:
                raise ValueError(f"Invalid value for arrange: '{arrange}'")
            for cluster_params in clusters_params:
                if stack:
                    yield ClusterMultiRelSeqGraph(table=self.table,
                                                  rels=rels,
                                                  y_ratio=y_ratio,
                                                  **cluster_params)
                else:
                    for rel in rels:
                        yield ClusterSingleRelSeqGraph(table=self.table,
                                                       rels=rel,
                                                       y_ratio=y_ratio,
                                                       **cluster_params)
        else:
            logger.error(f"{self} cannot graph {self.table}")


# Base Sequence Graph ##################################################

class SeqGraph(CartesianGraph, OneTableSeqGraph, OneSampGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *args,
                 table: (RelPosTableLoader
                         | MaskPosTableLoader
                         | ClustPosTableLoader),
                 rels: str,
                 y_ratio: bool,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table
        self.rel_codes = rels
        self.y_ratio = y_ratio
        # Verify that the table type is valid.
        _ = self.source

    @cached_property
    def rels(self):
        """ List the relationships to graph. """
        return [REL_CODES[rel] for rel in self.rel_codes]

    @property
    def ncols(self):
        # There is only one column of plots in this graph.
        return 1

    @classmethod
    @abstractmethod
    def sources(cls) -> dict[type[RelTypeTableLoader], str]:
        """ Names of the sources of data. """

    @property
    def source(self):
        """ Source of the data. """
        table_type = type(self.table)
        try:
            return self.sources()[table_type]
        except KeyError:
            raise TypeError(f"Invalid table for {self}: {table_type.__name__}")

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @property
    def x_title(self):
        return POS_NAME

    @property
    def y_title(self):
        return "Ratio" if self.y_ratio else "Count"

    @property
    def title(self):
        rels = '/'.join(self.rels)
        return (f"{self.y_title} of {rels} bases in {self.source} reads "
                f"from {self.sample} per position in {self.ref}:{self.sect}")

    @property
    @abstractmethod
    def subject(self):
        """ Subject of the data. """

    @property
    @abstractmethod
    def predicate(self):
        """ Predicate of the data. """

    @property
    def graph_filename(self):
        return f"{self.subject}_{self.predicate}".lower()


# Sequence Graphs by Source ############################################

class EnsembleSeqGraph(SeqGraph, ABC):

    @cached_property
    def data(self):
        return self.table.select(self.y_ratio, PRECISION, rels=self.rels)

    @classmethod
    def sources(cls):
        return {RelPosTableLoader: "Related", MaskPosTableLoader: "Masked"}

    @property
    def nrows(self):
        return 1

    @property
    def subplot_titles(self):
        return [ENSEMBLE_NAME]

    @property
    def subject(self):
        return self.source


class ClusterSeqGraph(SeqGraph, ABC):

    def __init__(self, *args, order: int = 0, cluster: int = 0, **kwargs):
        super().__init__(*args, **kwargs)
        validate_order_cluster(order, cluster, allow_zero=True)
        self._order = [order] if order > 0 else None
        self._cluster = [cluster] if cluster > 0 else None

    @cached_property
    def data(self):
        # Determine which relationships, orders, and clusters to use.
        select = dict(rels=self.rels)
        if self._order is not None:
            select.update(dict(order=self._order))
            if self._cluster is not None:
                select.update(dict(cluster=self._cluster))
        return self.table.select(self.y_ratio, PRECISION, **select)

    @cached_property
    def clusters(self):
        return pd.MultiIndex.from_arrays(
            [self.data.columns.get_level_values(level)
             for level in ORD_CLS_NAME],
            names=ORD_CLS_NAME
        ).drop_duplicates()

    @property
    def nrows(self):
        return self.clusters.size

    @property
    def subplot_titles(self):
        level_values = self.clusters.get_level_values
        return [fmt_clust_name(order, cluster)
                for order, cluster in zip(level_values(ORD_NAME),
                                          level_values(CLS_NAME),
                                          strict=True)]

    @classmethod
    def sources(cls):
        return {ClustPosTableLoader: "Clustered"}

    @property
    def subject(self):
        osymbol = 'x' if self._order is None else self._order[0]
        csymbol = 'x' if self._cluster is None else self._cluster[0]
        return f"{self.source}-{osymbol}-{csymbol}"


# Sequence Graphs by Series Type #######################################

class SingleRelSeqGraph(SeqGraph, ABC):
    """ Bar graph where each bar shows one relationship of the base. """

    @property
    def rel(self):
        """ Relationship of the identity sequence graph. """
        if len(self.rels) != 1:
            raise ValueError(
                f"Expected exactly 1 relationship, but got {len(self.rels)}")
        return self.rels[0]

    @property
    def predicate(self):
        return f"{self.rel_codes}-{self.y_title}"


class MultiRelSeqGraph(SeqGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    relationships for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    @cache
    def get_figure(self):
        fig = super().get_figure()
        # Stack the bars at each position.
        fig.update_layout(barmode="stack")
        return fig

    @property
    def predicate(self):
        return f"{self.rel_codes}-stack-{self.y_title}"


# Instantiable Sequence Graphs #########################################

class EnsembleSingleRelSeqGraph(EnsembleSeqGraph, SingleRelSeqGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        for trace in iter_base_traces(self.data.squeeze(axis=1), self.cmap):
            yield (1, 1), trace


class ClusterSingleRelSeqGraph(ClusterSeqGraph, SingleRelSeqGraph):

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_base_traces(values, self.cmap):
                yield (row, 1), trace


class EnsembleMultiRelSeqGraph(EnsembleSeqGraph, MultiRelSeqGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        data = self.data.copy()
        if data.columns.nlevels != 1:
            raise ValueError(
                f"Expected 1 level of columns, but got {data.columns.names}")
        # Replace the columns with a single index.
        data.columns = data.columns.get_level_values(REL_NAME)
        for trace in iter_stack_traces(self.data, self.cmap):
            yield (1, 1), trace


class ClusterMultiRelSeqGraph(ClusterSeqGraph, MultiRelSeqGraph):

    def get_traces(self):
        for row, ok in enumerate(self.clusters, start=1):
            for trace in iter_stack_traces(self.data.loc[:, ok], self.cmap):
                yield (row, 1), trace
