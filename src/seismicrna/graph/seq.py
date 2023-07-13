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
from ..core import docdef
from ..core.cli import (opt_input_file, opt_rels, opt_stack, opt_yfrac,
                        opt_csv, opt_html, opt_pdf, opt_max_procs, opt_parallel,
                        SUBPLOT_CLUST, SUBPLOT_ORDER, SUBPLOT_NONE)
from ..core.parallel import dispatch
from ..core.sect import BASE_NAME, POS_NAME
from ..core.seq import BASES
from ..table.base import REL_CODES
from ..table.load import (RelTypeTableLoader, RelPosTableLoader,
                          MaskPosTableLoader, ClustPosTableLoader)

logger = getLogger(__name__)

# Number of digits to which to round decimals.

params = [
    opt_input_file,
    opt_rels,
    opt_stack,
    opt_yfrac,
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
        subplot: str,
        stack: bool,
        yf: bool,
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
                                            subplot=subplot, yfrac=yf,
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
    rel = str(data.name)
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
    for field, series in data.items():
        yield get_stack_trace(series, cmap)


# Sequence Graph Writer ################################################

class SeqGraphWriter(GraphWriter):

    def iter(self, rels: str, subplot: str, stack: bool, yfrac: bool):
        if type(self.table) in PopAvgSeqGraph.sources():
            if stack:
                yield PopAvgMultiRelSeqGraph(table=self.table,
                                             rels=rels,
                                             yfrac=yfrac)
            else:
                for rel in rels:
                    yield PopAvgSingleRelSeqGraph(table=self.table,
                                                  rels=rel,
                                                  yfrac=yfrac)
        elif type(self.table) in ClustSeqGraph.sources():
            if subplot == SUBPLOT_NONE:
                # Create one file per cluster.
                clusters_params = [dict(order=order, cluster=cluster)
                                   for order, cluster in self.table.ord_clust]
            elif subplot == SUBPLOT_ORDER:
                # Create one file per order: each cluster is a subplot.
                orders = sorted(self.table.orders)
                clusters_params = [dict(order=order) for order in orders]
            elif subplot == SUBPLOT_CLUST:
                # Create one file: each cluster is a subplot.
                clusters_params = [dict()]
            else:
                raise ValueError(f"Invalid value for subplot: '{subplot}'")
            for cluster_params in clusters_params:
                if stack:
                    yield ClustMultiRelSeqGraph(table=self.table,
                                                rels=rels,
                                                yfrac=yfrac,
                                                **cluster_params)
                else:
                    for rel in rels:
                        yield ClustSingleRelSeqGraph(table=self.table,
                                                     rels=rel,
                                                     yfrac=yfrac,
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
                 yfrac: bool,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table
        self.rels = rels
        self.yratio = yfrac
        # Verify that the table type is valid.
        _ = self.source

    @property
    def nrows(self):
        # Each column of the data is plotted on a separate row.
        return self.data.shape[1]

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
    def xattr(self):
        return POS_NAME

    @property
    def yattr(self):
        return "Ratio" if self.yratio else "Count"

    @property
    def title(self):
        fields = '/'.join(sorted(REL_CODES[c] for c in self.rels))
        return (f"{self.yattr} of {fields} bases in {self.source} reads "
                f"from {self.sample} per position in {self.ref}:{self.sect}")

    @property
    def graph_filename(self):
        return f"{self.source}_{self.rels}_{self.yattr}".lower()


# Sequence Graphs by Source ############################################

class PopAvgSeqGraph(SeqGraph, ABC):

    @cached_property
    def data(self):
        return self.table.select(self.yratio, PRECISION, rels=self.rels)

    @classmethod
    def sources(cls):
        return {RelPosTableLoader: "Related", MaskPosTableLoader: "Masked"}

    @property
    def nrows(self):
        return 1


class ClustSeqGraph(SeqGraph, ABC):

    def __init__(self, *args, order: int = 0, cluster: int = 0, **kwargs):
        super().__init__(*args, **kwargs)
        self._orders = [order] if order > 0 else None
        self._clusters = [cluster] if cluster > 0 else None

    @cached_property
    def data(self):
        return self.table.select(self.yratio, PRECISION,
                                 rels=self.rels,
                                 orders=self._orders,
                                 clusters=self._clusters)

    @classmethod
    def sources(cls):
        return {ClustPosTableLoader: "Clustered"}


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


class MultiRelSeqGraph(SeqGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    relationships for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def get_traces(self):
        yield from iter_stack_traces(self.data, self.cmap)

    @cache
    def get_figure(self):
        fig = super().get_figure()
        # Stack the bars at each position.
        fig.update_layout(barmode="stack")
        return fig


# Instantiable Sequence Graphs #########################################

class PopAvgSingleRelSeqGraph(PopAvgSeqGraph, SingleRelSeqGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        for trace in iter_base_traces(self.data.squeeze(axis=1), self.cmap):
            yield 0, 0, trace


class ClustSingleRelSeqGraph(ClustSeqGraph, SingleRelSeqGraph):

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items()):
            for trace in iter_base_traces(values, self.cmap):
                yield row, 0, trace


class PopAvgMultiRelSeqGraph(PopAvgSeqGraph, MultiRelSeqGraph):
    pass


class ClustMultiRelSeqGraph(ClustSeqGraph, MultiRelSeqGraph):
    pass
