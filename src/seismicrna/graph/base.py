from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Any, Callable, Generator, Iterable

import pandas as pd
from click import Argument, Option
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from ..core import path
from ..core.arg import (CLUST_INDIV,
                        CLUST_ORDER,
                        CLUST_UNITE,
                        arg_input_path,
                        opt_rels,
                        opt_use_ratio,
                        opt_quantile,
                        opt_arrange,
                        opt_out_dir,
                        opt_csv,
                        opt_html,
                        opt_pdf,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.header import Header, format_clust_names
from ..core.seq import DNA
from ..core.write import need_write
from ..table.base import (Table,
                          PosTable,
                          RelTable,
                          MaskTable,
                          ClustTable)

logger = getLogger(__name__)

# Define sources.
RELATED = "related"
MASKED = "masked"
CLUSTERED = "clustered"

# String to join sample names.
LINKER = "__and__"


def _write_graph(writer: Callable[[Path], Any],
                 file: Path,
                 force: bool = False):
    """ Write an image or raw data for a graph to a file. """
    if need_write(file, force):
        writer(file)
    return file


def make_index(header: Header, order: int | None, clust: int | None):
    """ Make an index for the rows or columns of a graph. """
    if header.max_order == 0:
        # If there are no clusters, then no clusters must be selected.
        if order or clust:
            raise ValueError(f"Cannot select orders or clusters from {header}")
        return header.clusts
    # If there are any relationship names in the index, then drop them
    # and then select the order(s) and cluster(s) for the index.
    return header.modified(rels=()).select(order=order, clust=clust)


def _index_size(index: pd.Index | None):
    return index.size if index is not None else 1


def _index_titles(index: pd.Index | None):
    return (format_clust_names(index, allow_zero=True, allow_duplicates=False)
            if index is not None
            else None)


def get_source_name(table: Table):
    if isinstance(table, RelTable):
        return RELATED
    if isinstance(table, MaskTable):
        return MASKED
    if isinstance(table, ClustTable):
        return CLUSTERED
    raise TypeError(f"Invalid table type: {type(table).__name__}")


def make_source_sample(source: str, sample: str):
    return f"{source} reads from sample {repr(sample)}"


def make_subject(source: str, order: int | None, clust: int | None):
    if source == RELATED or source == MASKED:
        if order or clust:
            raise ValueError(f"For {source} data, order and clust must both "
                             f"be 0 or None, but got {order} and {clust}")
        return source
    if source == CLUSTERED:
        return "-".join(map(str, [source,
                                  order if order is not None else "x",
                                  clust if clust is not None else "x"]))
    raise ValueError(f"Invalid data source: {repr(source)}")


def arrange_table(table: Table, arrange: str):
    if arrange == CLUST_INDIV:
        # One file per cluster, with no subplots.
        return [dict(order=order, clust=cluster)
                for order, cluster in table.header.clusts]
    elif arrange == CLUST_ORDER:
        # One file per order, with one subplot per cluster.
        return [dict(order=order) for order in sorted(table.header.orders)]
    elif arrange == CLUST_UNITE:
        # One file, with one subplot per cluster for all orders.
        return [dict()]
    raise ValueError(f"Invalid value for arrange: {repr(arrange)}")


class GraphBase(ABC):

    @classmethod
    @abstractmethod
    def graph_kind(cls) -> str:
        """ Kind of graph. """

    @classmethod
    @abstractmethod
    def what(cls) -> str:
        """ What is being graphed. """

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.SectSeg,
                path.GraphSeg)

    def __init__(self, *,
                 out_dir: str | Path,
                 rels: str,
                 use_ratio: bool,
                 quantile: float):
        """
        Parameters
        ----------
        out_dir: str | Path
            Path of the top-level output directory for all graph files.

        rels: str
            Relationship(s) whose data will be pulled from the table(s).
            Each is given as a one-letter code, the definitions of which
            are given in seismicrna.table.base.REL_CODES, as follows:

            - Covered: v
            - Informed: n
            - Matched: r
            - Mutated: m
            - Subbed: s
            - Subbed to A: a
            - Subbed to C: c
            - Subbed to G: g
            - Subbed to T: t
            - Deleted: d
            - Inserted: i

            More than one relationship can be specified by passing a
            string longer than one character; for example, `acgt` would
            mean to use all four types of substitutions.

        use_ratio: bool
            Use the ratio of the number of times the relationship occurs
            to the number of occurrances of another kind of relationship
            (which is Covered for Covered and Informed, and Informed for
            all other relationships), rather than the raw count.

        quantile: float
            If `use_ratio` is True, then normalize the ratios to this
            quantile and then winsorize them to the interval [0, 1].
            Passing 0.0 disables normalization and winsorization.
        """
        self.top = Path(out_dir)
        self.rel_codes = rels
        self.use_ratio = use_ratio
        self.quantile = quantile

    @property
    def data_kind(self):
        """ Kind of data being used: either "ratio" or "count". """
        return "ratio" if self.use_ratio else "count"

    @property
    @abstractmethod
    def source_sample(self) -> str:
        """ Source and sample of the data. """

    @property
    @abstractmethod
    def sample(self) -> str:
        """ Name(s) of the sample(s) from which the data come. """

    @property
    @abstractmethod
    def ref(self) -> str:
        """ Name of the reference sequence from which the data come. """

    @property
    @abstractmethod
    def sect(self) -> str:
        """ Name of the reference section from which the data come. """

    @property
    @abstractmethod
    def seq(self) -> DNA:
        """ Sequence of the section from which the data come. """

    @cached_property
    def details(self) -> list[str]:
        """ Additional details about the graph. """
        return ([f"quantile = {round(self.quantile, 3)}"] if self.use_ratio
                else list())

    @property
    @abstractmethod
    def subject(self):
        """ Subject of the graph. """

    @property
    def predicate(self):
        """ Predicate of the graph. """
        fields = [self.rel_codes, self.data_kind]
        if self.use_ratio:
            fields.append(f"q{round(self.quantile * 100.)}")
        return "-".join(fields)

    @cached_property
    def graph_filename(self):
        """ Name of the graph's output file, without its extension. """
        return "_".join([self.graph_kind(), self.subject, self.predicate])

    def get_path_fields(self):
        """ Path fields. """
        return {path.TOP: self.top,
                path.SAMP: self.sample,
                path.CMD: path.CMD_GRA_DIR,
                path.REF: self.ref,
                path.SECT: self.sect,
                path.GRAPH: self.graph_filename}

    def get_path(self, ext: str):
        """ Path to the output file of the graph. """
        return path.buildpar(*self.get_path_segs(),
                             **self.get_path_fields(),
                             ext=ext)

    @property
    @abstractmethod
    def rel_names(self):
        """ Names of the relationships to graph. """

    @cached_property
    def relationships(self) -> str:
        """ Relationships being graphed as a slash-separated string. """
        return "/".join(self.rel_names)

    @cached_property
    def _fetch_kwargs(self) -> dict[str, Any]:
        """ Keyword arguments for self._fetch_data. """
        return dict(rel=self.rel_names)

    def _fetch_data(self, table: PosTable, **kwargs):
        """ Fetch data from the table. """
        kwargs = self._fetch_kwargs | kwargs
        return (table.fetch_ratio(quantile=self.quantile, **kwargs)
                if self.use_ratio
                else table.fetch_count(**kwargs))

    @cached_property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        """ Data of the graph. """

    @abstractmethod
    def get_traces(self) -> Iterable[tuple[tuple[int, int], go.Trace]]:
        """ Data traces of the graph. """

    @property
    @abstractmethod
    def row_index(self) -> pd.Index | None:
        """ Index of rows of subplots. """

    @property
    @abstractmethod
    def col_index(self) -> pd.Index | None:
        """ Index of columns of subplots. """

    @property
    def nrows(self):
        """ Number of rows of subplots. """
        return _index_size(self.row_index)

    @property
    def ncols(self):
        """ Number of columns of subplots. """
        return _index_size(self.col_index)

    @cached_property
    def row_titles(self):
        """ Titles of the rows. """
        return _index_titles(self.row_index)

    @cached_property
    def col_titles(self):
        """ Titles of the columns. """
        return _index_titles(self.col_index)

    @property
    @abstractmethod
    def x_title(self) -> str:
        """ Title of the x-axis. """

    @property
    @abstractmethod
    def y_title(self) -> str:
        """ Title of the y-axis. """

    @property
    def _subplots_params(self):
        return dict(rows=self.nrows,
                    cols=self.ncols,
                    row_titles=self.row_titles,
                    column_titles=self.col_titles,
                    x_title=self.x_title,
                    y_title=self.y_title,
                    shared_xaxes="all",
                    shared_yaxes="all")

    def _figure_init(self):
        """ Initialize the figure. """
        return make_subplots(**self._subplots_params)

    def _figure_data(self, figure: go.Figure):
        """ Add data to the figure. """
        for (row, col), trace in self.get_traces():
            figure.add_trace(trace, row=row, col=col)

    def _figure_layout(self, figure: go.Figure):
        """ Update the figure's layout. """
        figure.update_layout(title=self.title,
                             plot_bgcolor="#ffffff",
                             paper_bgcolor="#ffffff")
        figure.update_xaxes(linewidth=1,
                            linecolor="#000000",
                            autorange=True)
        figure.update_yaxes(linewidth=1,
                            linecolor="#000000",
                            autorange=True)

    @cached_property
    def figure(self):
        """ Figure object. """
        figure = self._figure_init()
        self._figure_data(figure)
        self._figure_layout(figure)
        return figure

    def write_csv(self, force: bool):
        """ Write the graph's source data to a CSV file. """
        return _write_graph(self.data.to_csv,
                            self.get_path(ext=path.CSV_EXT),
                            force)

    def write_html(self, force: bool):
        """ Write the graph to an HTML file. """
        return _write_graph(self.figure.write_html,
                            self.get_path(ext=path.HTML_EXT),
                            force)

    def write_pdf(self, force: bool):
        """ Write the graph to a PDF file. """
        return _write_graph(self.figure.write_image,
                            self.get_path(ext=path.PDF_EXT),
                            force)

    def write(self, csv: bool, html: bool, pdf: bool, force: bool = False):
        """ Write the selected files. """
        files = list()
        if csv:
            files.append(self.write_csv(force))
        if html:
            files.append(self.write_html(force))
        if pdf:
            files.append(self.write_pdf(force))
        return files

    @cached_property
    def title(self):
        """ Title of the graph. """
        return " ".join(
            [f"{self.what()} of {self.data_kind}s "
             f"of {self.relationships} bases in {self.source_sample} "
             f"over reference {repr(self.ref)} section {repr(self.sect)}"]
            + ([f"({'; '.join(self.details)})"] if self.details else [])
        )


class GraphWriter(ABC):
    """ Write the proper graph(s) for the table(s). """

    @classmethod
    @abstractmethod
    def get_graph_type(cls, *args, **kwargs) -> type[GraphBase]:
        """ Type of the graph to write. """

    def __init__(self, *table_files: Path):
        self.table_files = list(table_files)

    @abstractmethod
    def iter(self, *args, **kwargs) -> Generator[GraphBase, None, None]:
        """ Yield every graph for the table. """

    def write(self,
              *args,
              csv: bool,
              html: bool,
              pdf: bool,
              force: bool,
              **kwargs):
        """ Generate and write every graph for the table. """
        return list(chain(graph.write(csv=csv, html=html, pdf=pdf, force=force)
                          for graph in self.iter(*args, **kwargs)))


class GraphRunner(ABC):

    @classmethod
    @abstractmethod
    def writer_type(cls) -> type[GraphWriter]:
        """ Type of GraphWriter. """

    @classmethod
    def universal_input_params(cls):
        """ Universal parameters controlling the input data. """
        return [arg_input_path,
                opt_rels,
                opt_quantile,
                opt_use_ratio]

    @classmethod
    def universal_output_params(cls):
        """ Universal parameters controlling the output graph. """
        return [opt_arrange,
                opt_out_dir,
                opt_csv,
                opt_html,
                opt_pdf,
                opt_force,
                opt_max_procs,
                opt_parallel]

    @classmethod
    def var_params(cls) -> list[Argument | Option]:
        """ Parameters that can vary among different classes. """
        return list()

    @classmethod
    def params(cls) -> list[Argument | Option]:
        """ Parameters for the command line. """
        return list(chain(cls.universal_input_params(),
                          cls.var_params(),
                          cls.universal_output_params()))

    @classmethod
    @abstractmethod
    def run(cls,
            input_path: tuple[str, ...],
            rels: tuple[str, ...],
            quantile: float,
            use_ratio: bool, *,
            arrange: str,
            out_dir: str,
            csv: bool,
            html: bool,
            pdf: bool,
            force: bool,
            max_procs: int,
            parallel: bool,
            **kwargs) -> list[Path]:
        """ Run graphing. """

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
