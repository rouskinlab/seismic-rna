import os
from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command
from plotly import graph_objects as go

from .color import RelColorMap, SeqColorMap
from .base import CartesianGraph, OneTableSeqGraph, PRECISION
from .seq import get_table_params
from .traces import iter_seq_base_bar_traces, iter_seq_stack_bar_traces
from .write import OneTableGraphWriter
from ..cluster.names import (ENSEMBLE_NAME, CLS_NAME, ORD_NAME, fmt_clust_name,
                             validate_order_cluster)
from ..core import docdef
from ..core.cli import (arg_input_file, opt_rels, opt_y_ratio, opt_quantile,
                        opt_arrange, opt_csv, opt_html, opt_pdf,
                        opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.sect import POS_NAME
from ..table.base import REL_CODES, REL_NAME
from ..table.load import (PosTableLoader, RelPosTableLoader, MaskPosTableLoader,
                          ClustPosTableLoader, find_tables, get_clusters)

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

# Number of digits to which to round decimals.

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


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positions in a sequence. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_file: tuple[str, ...],
        rels: tuple[str, ...], *,
        y_ratio: bool,
        quantile: float,
        arrange: str,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph seqbar module. """
    writers = list(map(SeqBarGraphWriter, find_tables(input_file)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels_sets=rels, y_ratio=y_ratio,
                                            quantile=quantile, arrange=arrange,
                                            csv=csv, html=html, pdf=pdf))))


# Sequence Bar Graph Writer ############################################

class SeqBarGraphWriter(OneTableGraphWriter):

    def iter(self, rels_sets: tuple[str, ...],
             y_ratio: bool, quantile: float, arrange: str):
        stype, mtype, csparams = get_table_params(self.table, arrange,
                                                  EnsembleSingleRelSeqBarGraph,
                                                  ClusterSingleRelSeqBarGraph,
                                                  EnsembleMultiRelSeqBarGraph,
                                                  ClusterMultiRelSeqBarGraph)
        for cparams in csparams:
            for rels in rels_sets:
                graph_type = stype if len(rels) == 1 else mtype
                yield graph_type(table=self.table, rels=rels,
                                 y_ratio=y_ratio, quantile=quantile,
                                 **cparams)


# Base Sequence Bar Graph ##############################################

class SeqBarGraph(CartesianGraph, OneTableSeqGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *args, rels: str, y_ratio: bool, quantile: float,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.rel_codes = rels
        self.y_ratio = y_ratio
        self.quantile = quantile
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
    def sources(cls) -> dict[type[PosTableLoader], str]:
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
    def predicate(self):
        """ Predicate of the data. """
        return f"{self.rel_codes}-{self.y_title}"

    @property
    def graph_filename(self):
        return f"{self.subject}_{self.predicate}_{COMMAND}".lower()

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


# Sequence Graphs by Source ############################################

class EnsembleSeqBarGraph(SeqBarGraph, ABC):

    @cached_property
    def data(self):
        return self.table.select(self.y_ratio, self.quantile, PRECISION,
                                 rels=self.rels)

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


class ClusterSeqBarGraph(SeqBarGraph, ABC):

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
        return self.table.select(self.y_ratio, self.quantile, PRECISION,
                                 **select)

    @cached_property
    def clusters(self):
        return get_clusters(self.data.columns)

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

class SingleRelSeqBarGraph(SeqBarGraph, ABC):
    """ Bar graph where each bar shows one relationship of the base. """

    @property
    def rel(self):
        """ Relationship of the identity sequence graph. """
        if len(self.rels) != 1:
            raise ValueError(
                f"Expected exactly 1 relationship, but got {len(self.rels)}")
        return self.rels[0]


class MultiRelSeqBarGraph(SeqBarGraph, ABC):
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


# Instantiable Sequence Graphs #########################################

class EnsembleSingleRelSeqBarGraph(EnsembleSeqBarGraph, SingleRelSeqBarGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        for trace in iter_seq_base_bar_traces(self.data.squeeze(axis=1), self.cmap):
            yield (1, 1), trace


class ClusterSingleRelSeqBarGraph(ClusterSeqBarGraph, SingleRelSeqBarGraph):

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_seq_base_bar_traces(values, self.cmap):
                yield (row, 1), trace


class EnsembleMultiRelSeqBarGraph(EnsembleSeqBarGraph, MultiRelSeqBarGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        data = self.data.copy()
        if data.columns.nlevels != 1:
            raise ValueError(
                f"Expected 1 level of columns, but got {data.columns.names}")
        # Replace the columns with a single index.
        data.columns = data.columns.get_level_values(REL_NAME)
        for trace in iter_seq_stack_bar_traces(self.data, self.cmap):
            yield (1, 1), trace


class ClusterMultiRelSeqBarGraph(ClusterSeqBarGraph, MultiRelSeqBarGraph):

    def get_traces(self):
        for row, ok in enumerate(self.clusters, start=1):
            for trace in iter_seq_stack_bar_traces(self.data.loc[:, ok], self.cmap):
                yield (row, 1), trace
