import os
from functools import cached_property
from itertools import chain, product
from logging import getLogger
from pathlib import Path

from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairTwoAxisGraph, SeqPairGraphWriter, get_titles
from .traces import iter_seq_base_scatter_traces
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_rels,
                        opt_y_ratio,
                        opt_quantile,
                        opt_arrange,
                        opt_csv,
                        opt_html,
                        opt_pdf,
                        opt_max_procs,
                        opt_parallel)
from ..core.parallel import dispatch
from ..table.load import find_tables

logger = getLogger(__name__)

# Number of digits to which to round decimals.

params = [
    arg_input_path,
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

COMMAND = __name__.split(os.path.extsep)[-1]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Create scatter plots between pairs of samples at each position
    in a sequence. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_path: tuple[str, ...],
        rels: tuple[str, ...], *,
        y_ratio: bool,
        quantile: float,
        arrange: str,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    tables = list(find_tables(input_path))
    if len(tables) % 2 != 0:
        raise ValueError(f"Number of files must be even, but got {len(tables)}")
    writers = [SeqScatterGraphWriter(table1_file=t1, table2_file=t2)
               for t1, t2 in zip(tables[0::2], tables[1::2], strict=True)]
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels_sets=rels, y_ratio=y_ratio,
                                            quantile=quantile, arrange=arrange,
                                            csv=csv, html=html, pdf=pdf))))


class SeqScatterGraph(SeqPairTwoAxisGraph):

    @classmethod
    def graph_type(cls):
        return COMMAND

    @property
    def x_title(self) -> str:
        return self.sample1

    @property
    def y_title(self):
        return self.sample2

    @cached_property
    def col_titles(self):
        return get_titles(source=self.source1, clusters=self.clusters1)

    @cached_property
    def row_titles(self):
        return get_titles(source=self.source2, clusters=self.clusters2)

    def get_traces(self):
        for (col, (_, vals1)), (row, (_, vals2)) in product(
                enumerate(self.data1.items(), start=1),
                enumerate(self.data2.items(), start=1)
        ):
            for trace in iter_seq_base_scatter_traces(vals1, vals2, self.cmap):
                yield (row, col), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_xaxes(gridcolor="#d0d0d0")
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqScatterGraphWriter(SeqPairGraphWriter):

    @property
    def graph_type(self):
        return SeqScatterGraph

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
