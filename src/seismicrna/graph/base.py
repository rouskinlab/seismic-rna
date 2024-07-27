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
from ..core.arg import (NO_GROUP,
                        GROUP_BY_K,
                        GROUP_ALL,
                        arg_input_path,
                        opt_rels,
                        opt_use_ratio,
                        opt_quantile,
                        opt_cgroup,
                        opt_csv,
                        opt_html,
                        opt_svg,
                        opt_pdf,
                        opt_png,
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
from ..table.load import (find_pos_tables,
                          find_read_tables,
                          load_pos_table,
                          load_read_table)

logger = getLogger(__name__)

# Define actions.
ACTION_REL = "related"
ACTION_MASK = "masked"
ACTION_CLUST = "clustered"

# String to join sample names.
LINKER = "__and__"


def make_tracks(header: Header, k: int | None, clust: int | None):
    """ Make an index for the rows or columns of a graph. """
    if not header.clustered():
        # If there are no clusters, then no clusters must be selected.
        if k or clust:
            raise ValueError(f"Cannot select ks or clusters from {header}")
        return header.clusts
    # If there are any relationship names in the index, then drop them
    # and then select the k(s) and cluster(s) for the index.
    return header.get_clust_header().select(k=k, clust=clust).to_list()


def _track_count(tracks: list[tuple[int, int]] | None):
    return len(tracks) if tracks is not None else 1


def _track_titles(tracks: list[tuple[int, int]] | None):
    return (format_clust_names(tracks, allow_duplicates=False)
            if tracks is not None
            else None)


def get_action_name(table: Table):
    if isinstance(table, RelTable):
        return ACTION_REL
    if isinstance(table, MaskTable):
        return ACTION_MASK
    if isinstance(table, ClustTable):
        return ACTION_CLUST
    raise TypeError(f"Invalid table type: {type(table).__name__}")


def make_title_action_sample(action: str, sample: str):
    return f"{action} reads from sample {repr(sample)}"


def make_path_subject(action: str, k: int | None, clust: int | None):
    if action == ACTION_REL or action == ACTION_MASK:
        if k or clust:
            raise ValueError(f"For {action} data, k and clust must both "
                             f"be 0 or None, but got {k} and {clust}")
        return action
    if action == ACTION_CLUST:
        return "-".join(map(str, [action,
                                  k if k is not None else "x",
                                  clust if clust is not None else "x"]))
    raise ValueError(f"Invalid action: {repr(action)}")


def cgroup_table(table: Table, cgroup: str):
    if cgroup == NO_GROUP:
        # One file per cluster, with no subplots.
        return [dict(k=k, clust=clust)
                for k, clust in table.header.clusts]
    elif cgroup == GROUP_BY_K:
        # One file per k, with one subplot per cluster.
        return [dict(k=k, clust=None)
                for k in sorted(table.header.ks)]
    elif cgroup == GROUP_ALL:
        # One file, with one subplot per cluster.
        return [dict(k=None, clust=None)]
    raise ValueError(f"Invalid value for cgroup: {repr(cgroup)}")


class Annotation(object):

    def __init__(self,
                 row: int,
                 col: int,
                 x: float,
                 y: float,
                 text: str,
                 **kwargs):
        self.row = row
        self.col = col
        self.x = x
        self.y = y
        self.text = text
        self.kwargs = kwargs


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
                 use_ratio: bool,
                 quantile: float):
        """
        Parameters
        ----------
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
        self.use_ratio = use_ratio
        self.quantile = quantile

    @property
    @abstractmethod
    def codestring(self):
        """ String of the relationship code(s). """

    @property
    def data_kind(self):
        """ Kind of data being used: either "ratio" or "count". """
        return "ratio" if self.use_ratio else "count"

    @property
    @abstractmethod
    def title_action_sample(self) -> str:
        """ Action and sample for the title. """

    @property
    @abstractmethod
    def top(self) -> Path:
        """ Path of the top-level output directory for all files. """

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
    def path_subject(self):
        """ Subject of the graph. """

    @property
    def predicate(self):
        """ Predicate of the graph. """
        fields = [self.codestring, self.data_kind]
        if self.use_ratio:
            fields.append(f"q{round(self.quantile * 100.)}")
        return "-".join(fields)

    @cached_property
    def graph_filename(self):
        """ Name of the graph's output file, without its extension. """
        return "_".join([self.graph_kind(), self.path_subject, self.predicate])

    def get_path_fields(self):
        """ Path fields. """
        return {path.TOP: self.top,
                path.SAMP: self.sample,
                path.CMD: path.CMD_GRAPH_DIR,
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
    def row_tracks(self) -> list[tuple[int, int]] | None:
        """ Track for each row of subplots. """

    @property
    @abstractmethod
    def col_tracks(self) -> list[tuple[int, int]] | None:
        """ Track for each column of subplots. """

    @property
    def nrows(self):
        """ Number of rows of subplots. """
        return _track_count(self.row_tracks)

    @property
    def ncols(self):
        """ Number of columns of subplots. """
        return _track_count(self.col_tracks)

    @cached_property
    def row_titles(self):
        """ Titles of the rows. """
        return _track_titles(self.row_tracks)

    @cached_property
    def col_titles(self):
        """ Titles of the columns. """
        return _track_titles(self.col_tracks)

    @property
    @abstractmethod
    def x_title(self) -> str:
        """ Title of the x-axis. """

    @property
    @abstractmethod
    def y_title(self) -> str:
        """ Title of the y-axis. """

    @property
    def annotations(self) -> list[Annotation]:
        """ Text annotations for the figure. """
        return list()

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

    def _figure_annot(self, figure: go.Figure):
        """ Annotate the figure. """
        for annotation in self.annotations:
            figure.add_annotation(row=annotation.row,
                                  col=annotation.col,
                                  x=annotation.x,
                                  y=annotation.y,
                                  text=annotation.text,
                                  **annotation.kwargs)

    @cached_property
    def figure(self):
        """ Figure object. """
        figure = self._figure_init()
        self._figure_data(figure)
        self._figure_layout(figure)
        self._figure_annot(figure)
        return figure

    def write_csv(self, force: bool):
        """ Write the graph's source data to a CSV file. """
        file = self.get_path(path.CSV_EXT)
        if need_write(file, force):
            self.data.to_csv(file)
        return file

    def write_html(self, force: bool):
        """ Write the graph to an HTML file. """
        file = self.get_path(path.HTML_EXT)
        if need_write(file, force):
            self.figure.write_html(file)
        return file

    def _write_image(self, ext: str, force: bool):
        """ Write the graph to an image file. """
        file = self.get_path(ext)
        if need_write(file, force):
            self.figure.write_image(file)
        return file

    def write_svg(self, force: bool):
        """ Write the graph to an SVG file. """
        return self._write_image(path.SVG_EXT, force)

    def write_pdf(self, force: bool):
        """ Write the graph to a PDF file. """
        return self._write_image(path.PDF_EXT, force)

    def write_png(self, force: bool):
        """ Write the graph to a PNG file. """
        return self._write_image(path.PNG_EXT, force)

    def write(self,
              csv: bool,
              html: bool,
              svg: bool,
              pdf: bool,
              png: bool,
              force: bool = False):
        """ Write the selected files. """
        files = list()
        if csv:
            files.append(self.write_csv(force))
        if html:
            files.append(self.write_html(force))
        if svg:
            files.append(self.write_svg(force))
        if pdf:
            files.append(self.write_pdf(force))
        if png:
            files.append(self.write_png(force))
        return files

    @cached_property
    def _title_main(self):
        """ Main part of the title, as a list. """
        return [f"{self.what()} of {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over reference {repr(self.ref)} "
                f"section {repr(self.sect)}"]

    @cached_property
    def _title_details(self):
        """ Details of the title, as a list. """
        return [f"({'; '.join(self.details)})"] if self.details else []

    @cached_property
    def title(self):
        """ Title of the graph. """
        return " ".join(self._title_main + self._title_details)


class GraphWriter(ABC):
    """ Write the proper graph(s) for the table(s). """

    @classmethod
    @abstractmethod
    def get_table_loader(cls) -> Callable[[Path], Table]:
        """ Function to load table files. """

    @classmethod
    def load_table_file(cls, table_file: Path):
        """ Load one table file. """
        loader = cls.get_table_loader()
        return loader(table_file)

    def __init__(self, *table_files: Path):
        self.table_files = list(table_files)

    @abstractmethod
    def iter_graphs(self, *args, **kwargs) -> Generator[GraphBase, None, None]:
        """ Yield every graph for the table. """

    def write(self,
              *args,
              csv: bool,
              html: bool,
              svg: bool,
              pdf: bool,
              png: bool,
              force: bool,
              **kwargs):
        """ Generate and write every graph for the table. """
        return list(chain(graph.write(csv=csv,
                                      html=html,
                                      svg=svg,
                                      pdf=pdf,
                                      png=png,
                                      force=force)
                          for graph in self.iter_graphs(*args, **kwargs)))


class PosGraphWriter(GraphWriter, ABC):

    @classmethod
    def get_table_loader(cls):
        return load_pos_table


class ReadGraphWriter(GraphWriter, ABC):

    @classmethod
    def get_table_loader(cls):
        return load_read_table


class GraphRunner(ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[GraphWriter]:
        """ Type of GraphWriter. """

    @classmethod
    def universal_input_params(cls):
        """ Universal parameters controlling the input data. """
        return [arg_input_path,
                opt_rels,
                opt_use_ratio,
                opt_quantile]

    @classmethod
    def universal_output_params(cls):
        """ Universal parameters controlling the output graph. """
        return [opt_cgroup,
                opt_csv,
                opt_html,
                opt_svg,
                opt_pdf,
                opt_png,
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
    def get_table_finder(cls) -> Callable[[tuple[str, ...]], Generator]:
        """ Function to find and filter table files. """

    @classmethod
    def list_table_files(cls, input_path: tuple[str, ...]):
        """ Find, filter, and list all table files from input files. """
        finder = cls.get_table_finder()
        return list(finder(input_path))

    @classmethod
    @abstractmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            rels: tuple[str, ...],
            use_ratio: bool,
            quantile: float,
            cgroup: str,
            csv: bool,
            html: bool,
            svg: bool,
            pdf: bool,
            png: bool,
            force: bool,
            max_procs: int,
            parallel: bool,
            **kwargs) -> list[Path]:
        """ Run graphing. """


class PosGraphRunner(GraphRunner, ABC):

    @classmethod
    def get_table_finder(cls):
        return find_pos_tables


class ReadGraphRunner(GraphRunner, ABC):

    @classmethod
    def get_table_finder(cls):
        return find_read_tables

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
