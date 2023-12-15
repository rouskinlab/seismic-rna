import os
from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command
from plotly import graph_objects as go

from .base import CartesianGraph, OneRefGraph, OneTableGraph, make_subject
from .color import RelColorMap
from .seq import get_table_params
from .traces import iter_seq_base_bar_traces, iter_seqbar_stack_traces
from .write import OneTableGraphWriter
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_rels,
                        opt_use_ratio,
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
from ..table.base import MUTAT_REL, RelPosTable, MaskPosTable, ClustPosTable
from ..table.load import find_table_files

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


# Receiver Operating Characteristic Graph Writer #######################

class ROCGraphWriter(OneTableGraphWriter):

    def iter(self, arrange: str):
        _, _, csparams = get_table_params(self.table, arrange)
        for cparams in csparams:
            yield ROCGraph(table=self.table, **cparams)


# Receiver Operating Characteristic Graph ##############################

class ROCGraph(CartesianGraph, OneRefGraph, OneTableGraph):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *,
                 order: int | None = None,
                 clust: int | None = None,
                 **kwargs):
        super().__init__(**kwargs)
        self._order = order
        self._clust = clust

    def iter_profiles(self):
        self.table.iter_profiles()

    @cached_property
    def data(self):
        return self.table.fetch_ratio(rel=MUTAT_REL)

    @property
    def x_title(self):
        return "False Positive Rate"

    @property
    def y_title(self):
        return "True Positive Rate"

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
        use_ratio: bool,
        quantile: float,
        arrange: str,
        csv: bool,
        html: bool,
        pdf: bool,
        force: bool,
        max_procs: int) -> list[Path]:
    """ Run the graph seqbar module. """
    table_paths = find_table_files(input_path)
    ct_paths = None
    writers = list(map(SeqBarGraphWriter))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs,
                                parallel=True,
                                pass_n_procs=False,
                                kwargs=dict(rels_sets=rels,
                                            use_ratio=use_ratio,
                                            quantile=quantile,
                                            arrange=arrange,
                                            csv=csv,
                                            html=html,
                                            pdf=pdf,
                                            force=force))))


params = [
    arg_input_path,
    opt_rels,
    opt_use_ratio,
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
