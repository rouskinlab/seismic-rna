import os
from functools import cached_property
from itertools import chain, product
from logging import getLogger
from pathlib import Path

from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraph, SeqPairGraphWriter, get_titles
from .traces import iter_base_bar_traces
from ..core import docdef
from ..core.cli import (arg_input_file, opt_rels, opt_y_ratio, opt_quantile,
                        opt_arrange, opt_csv, opt_html, opt_pdf,
                        opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.sect import POS_NAME
from ..table.load import (find_tables)

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
    """ Create bar graphs of differences between pairs of samples at
    each position in a sequence. """
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
    """ Run the graph pos module. """
    tables = list(find_tables(input_file))
    if len(tables) % 2 != 0:
        raise ValueError(f"Number of files must be even, but got {len(tables)}")
    writers = [SeqDeltaGraphWriter(table1_file=t1, table2_file=t2)
               for t1, t2 in zip(tables[0::2], tables[1::2], strict=True)]
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels_sets=rels,  y_ratio=y_ratio,
                                            quantile=quantile, arrange=arrange,
                                            csv=csv, html=html, pdf=pdf))))


class SeqDeltaGraph(SeqPairGraph):

    @classmethod
    def graph_type(cls):
        return COMMAND

    @property
    def x_title(self) -> str:
        return POS_NAME

    @property
    def y_title(self):
        return f"{self.quantity}-2 minus {self.quantity}-1"

    @cached_property
    def col_titles(self):
        return get_titles(sample=self.sample1, source=self.source1,
                          clusters=self.clusters1)

    @cached_property
    def row_titles(self):
        return get_titles(sample=self.sample2, source=self.source2,
                          clusters=self.clusters2)

    def get_traces(self):
        for (col, (_, vals1)), (row, (_, vals2)) in product(
                enumerate(self.data1.items(), start=1),
                enumerate(self.data2.items(), start=1)
        ):
            for trace in iter_base_bar_traces(vals2 - vals1, self.cmap):
                yield (row, col), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqDeltaGraphWriter(SeqPairGraphWriter):

    @property
    def graph_type(self):
        return SeqDeltaGraph
