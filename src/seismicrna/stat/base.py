from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from .color import ColorMap, get_cmap
from ..core import path
from ..core.header import format_clust_names
from ..core.seq import DNA
from ..core.write import need_write
from ..table.base import Table, PosTable

logger = getLogger(__name__)

# Number of digits behind the decimal point to be kept.
PRECISION = 6


def _write_graph(writer: Callable[[Path], Any],
                 file: Path,
                 force: bool = False):
    """ Write an image or raw data for a graph to a file. """
    if need_write(file, force):
        writer(file)
    return file


def _index_size(index: pd.Index | None):
    return index.size if index is not None else 1


def _index_titles(index: pd.Index | None):
    return format_clust_names(index) if index is not None else None


def make_subject(source: str, order: int | None, clust: int | None):
    return "-".join(map(str, [source,
                              order if order is not None else "x",
                              clust if clust is not None else "x"]))


class GraphBase(ABC):

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return path.CmdSeg, path.GraphSeg

    @classmethod
    def col_row_sep(cls):
        """ Separator between column and row title. """
        return " "

    def __init__(self, *, cmap: str | None = None):
        self._cmap_name = cmap

    @cached_property
    @abstractmethod
    def data(self) -> Any:
        """ Data of the graph. """

    @property
    @abstractmethod
    def title(self) -> str:
        """ Title of the graph. """

    @classmethod
    @abstractmethod
    def get_cmap_type(cls) -> type[ColorMap]:
        """ Type of the color map. """

    @property
    def cmap(self) -> ColorMap:
        """ Color map of the graph. """
        return get_cmap(self.get_cmap_type(), self._cmap_name)

    @abstractmethod
    def get_traces(self) -> Iterable[tuple[tuple[int, int], go.Trace]]:
        """ Data traces of the graph. """

    @property
    @abstractmethod
    def top(self) -> Path:
        """ Output directory. """

    @property
    @abstractmethod
    def graph_filename(self):
        """ Name of the graph file. """

    def get_path_fields(self):
        """ Path fields. """
        return {path.TOP: self.top,
                path.CMD: path.CMD_GRA_DIR,
                path.GRAPH: self.graph_filename}

    def get_path(self, ext: str):
        """ Path to the output file of the graph. """
        return path.buildpar(*self.get_path_segs(),
                             **self.get_path_fields(),
                             ext=ext)

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
    def _subplots_params(self):
        return dict(rows=self.nrows,
                    cols=self.ncols,
                    row_titles=self.row_titles,
                    column_titles=self.col_titles)

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


class SampleGraph(GraphBase, ABC):
    """ Graph of one or more samples. """

    @property
    @abstractmethod
    def sample(self):
        """ Name of the sample. """
        return ""

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return super().get_path_segs() + (path.SampSeg,)

    def get_path_fields(self):
        return super().get_path_fields() | {path.SAMP: self.sample}


class OneSampleGraph(SampleGraph, ABC):
    """ Graph of one sample. """


class TwoSampleGraph(SampleGraph, ABC):
    """ Graph of two samples. """

    @property
    @abstractmethod
    def sample1(self) -> str:
        """ Name of the first sample. """

    @property
    @abstractmethod
    def sample2(self) -> str:
        """ Name of the second sample. """

    @property
    def sample(self):
        return f"{self.sample1}_vs_{self.sample2}"


class OneRefGraph(GraphBase, ABC):
    """ Graph of one reference.  """

    @property
    @abstractmethod
    def ref(self) -> str:
        """ Name of the reference sequence. """

    @property
    @abstractmethod
    def sect(self) -> str:
        """ Name of the section of the reference sequence. """

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return super().get_path_segs() + (path.RefSeg, path.SectSeg)

    def get_path_fields(self):
        return super().get_path_fields() | {path.REF: self.ref,
                                            path.SECT: self.sect}


class OneSeqGraph(OneRefGraph, ABC):
    """ Graph of one reference with an explicit sequence. """

    @property
    @abstractmethod
    def seq(self) -> DNA:
        """ Reference sequence as a DNA object. """


class OneTableGraph(OneSampleGraph, OneRefGraph, ABC):
    """ Graph of data from one TableLoader. """

    def __init__(self, *, table: Table | PosTable, **kwargs):
        super().__init__(**kwargs)
        if not isinstance(table, self.get_table_type()):
            raise TypeError(f"{type(self).__name__} expected table "
                            f"of type '{self.get_table_type().__name__}', "
                            f"but got type '{type(table).__name__}'")
        self._table = table

    @classmethod
    @abstractmethod
    def get_table_type(cls) -> type[Table | PosTable]:
        """ Type of TableLoader for this graph. """

    @property
    def table(self):
        """ Table of data. """
        return self._table

    @property
    def top(self):
        return self.table.top

    @property
    def sample(self):
        return self.table.sample

    @property
    def ref(self):
        return self.table.ref

    @property
    def sect(self):
        return self.table.sect

    @property
    def col_index(self):
        return None


class TwoTableGraph(TwoSampleGraph, OneRefGraph, ABC):

    def __init__(self, *,
                 table1: Table | PosTable,
                 table2: Table | PosTable,
                 **kwargs):
        super().__init__(**kwargs)
        for table, table_type in zip((table1, table2),
                                     (self.get_table1_type(),
                                      self.get_table2_type()),
                                     strict=True):
            if not isinstance(table, table_type):
                raise TypeError(f"{type(self).__name__} expected table "
                                f"of type '{table_type.__name__}', "
                                f"but got type '{type(table).__name__}'")
        self._table1 = table1
        self._table2 = table2

    @classmethod
    @abstractmethod
    def get_table1_type(cls) -> type[Table | PosTable]:
        """ Type of the first TableLoader for this graph. """

    @classmethod
    @abstractmethod
    def get_table2_type(cls) -> type[Table | PosTable]:
        """ Type of the second TableLoader for this graph. """

    @property
    def table1(self):
        """ First table of data. """
        return self._table1

    @property
    def table2(self):
        """ Second table of data. """
        return self._table2

    def _get_common_attribute(self, name: str):
        """ Get the common attribute for tables 1 and 2. """
        attr1 = getattr(self.table1, name)
        attr2 = getattr(self.table2, name)
        if attr1 != attr2:
            raise ValueError(f"Attribute {repr(name)} differs between tables "
                             f"1 ({repr(attr1)}) and 2 ({repr(attr2)})")
        return attr1

    @property
    def top(self):
        return self._get_common_attribute("top")

    @property
    def ref(self):
        return self._get_common_attribute("ref")

    @property
    def sect(self):
        return self._get_common_attribute("sect")

    @property
    def sample1(self):
        return self.table1.sample

    @property
    def sample2(self):
        return self.table2.sample

    @property
    def sample(self):
        return (self.sample1 if self.sample1 == self.sample2
                else f"{self.sample1}__and__{self.sample2}")


class OneTableSeqGraph(OneTableGraph, OneSeqGraph, ABC):

    @classmethod
    def get_table_type(cls):
        return PosTable

    @property
    def seq(self):
        return self.table.seq


class TwoTableSeqGraph(TwoTableGraph, OneSeqGraph, ABC):

    @classmethod
    def get_table1_type(cls):
        return PosTable

    @classmethod
    def get_table2_type(cls):
        return PosTable

    @property
    def seq(self):
        return self._get_common_attribute("seq")


class CartesianGraph(GraphBase, ABC):
    """ Graph with one pair of x and y axes per subplot. """

    @property
    @abstractmethod
    def x_title(self) -> str:
        """ Title of the x-axis. """

    @property
    @abstractmethod
    def y_title(self) -> str:
        """ Title of the y-axis. """

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        figure.update_xaxes(linewidth=1,
                            linecolor="#000000",
                            autorange=True)
        figure.update_yaxes(linewidth=1,
                            linecolor="#000000",
                            autorange=True)

    @property
    def _subplots_params(self):
        return super()._subplots_params | dict(x_title=self.x_title,
                                               y_title=self.y_title,
                                               shared_xaxes="all",
                                               shared_yaxes="all")

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
