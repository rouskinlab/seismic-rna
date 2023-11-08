import os
from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command
from plotly import graph_objects as go

from .base import PRECISION, CartesianGraph, OneTableSeqGraph, make_subject
from .color import RelColorMap, SeqColorMap
from .seq import get_table_params
from .traces import iter_seq_base_bar_traces, iter_seqbar_stack_traces
from .write import OneTableGraphWriter
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_rels,
                        opt_y_ratio,
                        opt_quantile,
                        opt_arrange,
                        opt_csv,
                        opt_html,
                        opt_pdf,
                        opt_force,
                        opt_max_procs)
from ..core.header import make_header
from ..core.parallel import dispatch
from ..core.seq import POS_NAME
from ..table.base import get_rel_name
from ..table.load import (PosTableLoader,
                          RelPosTableLoader,
                          MaskPosTableLoader,
                          ClustPosTableLoader,
                          find_tables)

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


# Sequence Bar Graph Writer ############################################

class SeqBarGraphWriter(OneTableGraphWriter):

    def iter(self,
             rels_sets: tuple[str, ...],
             y_ratio: bool,
             quantile: float,
             arrange: str):
        stype, mtype, csparams = get_table_params(self.table,
                                                  arrange,
                                                  AverageSingleRelSeqBarGraph,
                                                  ClusterSingleRelSeqBarGraph,
                                                  AverageMultiRelSeqBarGraph,
                                                  ClusterMultiRelSeqBarGraph)
        for cparams in csparams:
            for rels in rels_sets:
                graph_type = stype if len(rels) == 1 else mtype
                yield graph_type(table=self.table,
                                 rels=rels,
                                 y_ratio=y_ratio,
                                 quantile=quantile,
                                 **cparams)


# Base Sequence Bar Graph ##############################################

class SeqBarGraph(CartesianGraph, OneTableSeqGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    def __init__(self, *, rels: str, y_ratio: bool, quantile: float, **kwargs):
        super().__init__(**kwargs)
        self.rel_codes = rels
        self.y_ratio = y_ratio
        self.quantile = quantile

    @property
    def rels(self):
        """ List of relationships to graph. """
        return list(map(get_rel_name, self.rel_codes))

    @property
    def _fetch_kwargs(self):
        """ Keyword arguments for `self.table.fetch()`. """
        return dict(rel=self.rels)

    def _fetch_data(self, **kwargs):
        all_kwargs = self._fetch_kwargs | kwargs
        return (self.table.fetch_ratio(quantile=self.quantile,
                                       precision=PRECISION,
                                       **all_kwargs)
                if self.y_ratio
                else self.table.fetch_count(**all_kwargs))

    @cached_property
    def data(self):
        return self._fetch_data()

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

    @property
    def x_title(self):
        return POS_NAME

    @property
    def y_title(self):
        return "Ratio" if self.y_ratio else "Count"

    @property
    def title(self):
        rels = ", ".join(self.rels)
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

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        figure.update_yaxes(gridcolor="#d0d0d0")


# Sequence Graphs by Source ############################################

class AverageSeqBarGraph(SeqBarGraph, ABC):

    @property
    def row_index(self):
        return None

    @classmethod
    def sources(cls):
        return {RelPosTableLoader: "Related", MaskPosTableLoader: "Masked"}

    @property
    def subject(self):
        return self.source


class ClusterSeqBarGraph(SeqBarGraph, ABC):

    def __init__(self, *,
                 order: int | None = None,
                 clust: int | None = None,
                 **kwargs):
        super().__init__(**kwargs)
        self._order = order
        self._clust = clust

    @cached_property
    def row_index(self):
        return make_header(
            max_order=self.table.header.max_order,
            min_order=self.table.header.min_order
        ).select(
            order=self._order,
            clust=self._clust
        )

    @property
    def _fetch_kwargs(self):
        return super()._fetch_kwargs | dict(order=self._order,
                                            clust=self._clust)

    @classmethod
    def sources(cls):
        return {ClustPosTableLoader: "Clustered"}

    @property
    def subject(self):
        return make_subject(self.source, self._order, self._clust)


# Sequence Graphs by Series Type #######################################

class SingleRelSeqBarGraph(SeqBarGraph, ABC):
    """ Bar graph where each bar shows one relationship of the base. """


class MultiRelSeqBarGraph(SeqBarGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    relationships for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        # Stack the bars at each position.
        figure.update_layout(barmode="stack")


# Instantiable Sequence Graphs #########################################

class AverageSingleRelSeqBarGraph(AverageSeqBarGraph, SingleRelSeqBarGraph):

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


class AverageMultiRelSeqBarGraph(AverageSeqBarGraph, MultiRelSeqBarGraph):

    def get_traces(self):
        if self.nrows != 1:
            raise ValueError(f"Expected 1 series of data, but got {self.nrows}")
        data = self.data.copy()
        if data.columns.nlevels != 1:
            raise ValueError(
                f"Expected 1 level of columns, but got {data.columns.names}")
        for trace in iter_seqbar_stack_traces(self.data, self.cmap):
            yield (1, 1), trace


class ClusterMultiRelSeqBarGraph(ClusterSeqBarGraph, MultiRelSeqBarGraph):

    def get_traces(self):
        for row, (order, clust) in enumerate(self.row_index, start=1):
            for trace in iter_seqbar_stack_traces(self._fetch_data(order=order,
                                                                   clust=clust),
                                                  self.cmap):
                yield (row, 1), trace


@docdef.auto()
def run(input_path: tuple[str, ...],
        rels: tuple[str, ...], *,
        y_ratio: bool,
        quantile: float,
        arrange: str,
        csv: bool,
        html: bool,
        pdf: bool,
        force: bool,
        max_procs: int) -> list[Path]:
    """ Run the graph seqbar module. """
    writers = list(map(SeqBarGraphWriter, find_tables(input_path)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs,
                                parallel=True,
                                pass_n_procs=False,
                                kwargs=dict(rels_sets=rels,
                                            y_ratio=y_ratio,
                                            quantile=quantile,
                                            arrange=arrange,
                                            csv=csv,
                                            html=html,
                                            pdf=pdf,
                                            force=force))))


params = [
    arg_input_path,
    opt_rels,
    opt_y_ratio,
    opt_quantile,
    opt_arrange,
    opt_csv,
    opt_html,
    opt_pdf,
    opt_force,
    opt_max_procs,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positions in a sequence. """
    return run(*args, **kwargs)

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
