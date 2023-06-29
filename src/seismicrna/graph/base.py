from abc import ABC, abstractmethod
from logging import getLogger
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd
from plotly import graph_objects as go

from .color import ColorMap, get_cmap
from ..core import path
from ..core.seq import DNA
from ..table.base import Table
from ..table.load import load, TableLoader, PosTableLoader

logger = getLogger(__name__)


# Number of digits behind the decimal point to be kept.
PRECISION = 6


def find_tables(tables: tuple[str, ...]):
    """ Return a file for each given file/directory of a table. """
    return path.find_files_chain(map(Path, tables), [path.MutTabSeg])


class GraphBase(ABC):

    def __init__(self, cmap: str | None = None):
        self._cmap_name = cmap

    @classmethod
    @abstractmethod
    def get_data_type(cls) -> type | tuple[type, ...]:
        """ Type of the data. """
        return tuple()

    @abstractmethod
    def _get_data(self) -> pd.DataFrame | pd.Series:
        """ Get the data of the graph. """
        return pd.DataFrame()

    @cached_property
    def data(self) -> Any:
        """ Data of the graph. """
        data = self._get_data()
        if not isinstance(data, self.get_data_type()):
            raise TypeError(f"{self.__class__.__name__} expected data "
                            f"of type '{self.get_data_type().__name__}', "
                            f"but got type '{type(data).__name__}'")
        return data

    @property
    @abstractmethod
    def title(self):
        """ Title of the graph. """
        return ""

    @classmethod
    @abstractmethod
    def get_cmap_type(cls):
        """ Type of the color map. """
        return ColorMap

    @property
    def cmap(self) -> ColorMap:
        """ Color map of the graph. """
        return get_cmap(self.get_cmap_type(), self._cmap_name)

    @abstractmethod
    def get_traces(self):
        """ Data traces of the graph. """
        return list()

    @property
    @abstractmethod
    def out_dir(self):
        """ Output directory. """
        return Path()

    @property
    @abstractmethod
    def graph_filename(self):
        return ""

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return path.ModSeg, path.GraphSeg

    def get_path_fields(self):
        """ Path fields. """
        return {path.TOP: self.out_dir,
                path.MOD: path.MOD_GRAPH,
                path.GRAPH: self.graph_filename}

    def get_path(self, ext: str):
        """ Path to the output file of the graph. """
        return path.buildpar(*self.get_path_segs(),
                             **self.get_path_fields(),
                             ext=ext)

    def _init_figure(self):
        """ Initialize the figure. """
        return go.Figure(data=self.get_traces())

    def _layout_figure(self, fig: go.Figure):
        """ Update the figure's layout. """
        fig.update_layout(title=self.title,
                          plot_bgcolor="#ffffff",
                          paper_bgcolor="#ffffff")
        fig.update_xaxes(linewidth=1,
                         linecolor="#000000",
                         autorange=True)
        fig.update_yaxes(gridcolor="#d0d0d0",
                         linewidth=1,
                         linecolor="#000000",
                         autorange=True)
        return fig

    @cache
    def get_figure(self):
        """ Figure object. """
        return self._layout_figure(self._init_figure())

    def write_csv(self):
        """ Write the graph's source data to a CSV file. """
        file = self.get_path(ext=path.CSV_EXT)
        if file.is_file():
            logger.warning(f"File exists: {file}")
        else:
            self.data.to_csv(file)
        return file

    def write_html(self):
        """ Write the graph to an HTML file. """
        file = self.get_path(ext=path.HTML_EXT)
        if file.is_file():
            logger.warning(f"File exists: {file}")
        else:
            self.get_figure().write_html(file)
        return file

    def write_pdf(self):
        """ Write the graph to a PDF file. """
        file = self.get_path(ext=path.PDF_EXT)
        if file.is_file():
            logger.warning(f"File exists: {file}")
        else:
            self.get_figure().write_image(file)
        return file

    def write(self, csv: bool, html: bool, pdf: bool):
        """ Write the selected files. """
        files = list()
        if csv:
            files.append(self.write_csv())
        if html:
            files.append(self.write_html())
        if pdf:
            files.append(self.write_pdf())
        return files


class OneSampGraph(GraphBase, ABC):
    """ Graph of one sample. """

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
        return {**super().get_path_fields(), path.SAMP: self.sample}


class OneRefGraph(GraphBase, ABC):
    """ Graph of one reference.  """

    @property
    @abstractmethod
    def ref(self):
        """ Name of the reference sequence. """
        return ""

    @property
    @abstractmethod
    def sect(self):
        """ Name of the section of the reference sequence. """
        return ""

    @classmethod
    def get_path_segs(cls):
        """ Path segments. """
        return super().get_path_segs() + (path.RefSeg, path.SectSeg)

    def get_path_fields(self):
        return {**super().get_path_fields(),
                path.REF: self.ref,
                path.SECT: self.sect}


class OneSeqGraph(OneRefGraph, ABC):
    """ Graph of one reference with an explicit sequence. """

    @property
    @abstractmethod
    def seq(self):
        """ Reference sequence as a DNA object. """
        return DNA(b"")


class OneTableGraph(OneSampGraph, OneRefGraph, ABC):
    """ Graph of data from one TableLoader. """

    def __init__(self):
        super().__init__()
        self._table = None

    @classmethod
    @abstractmethod
    def get_table_type(cls):
        """ Type of TableLoader for this graph. """
        return TableLoader

    @property
    def table(self) -> (TableLoader | PosTableLoader | Table):
        """ Table of data. """
        if self._table is None:
            raise TypeError("table not set")
        return self._table

    @table.setter
    def table(self, table: TableLoader):
        if not isinstance(table, self.get_table_type()):
            raise TypeError(f"{self.__class__.__name__} expected table "
                            f"of type '{self.get_table_type().__name__}', "
                            f"but got type '{type(table).__name__}'")
        self._table = table

    @property
    def out_dir(self):
        return self.table.out_dir

    @property
    def sample(self):
        return self.table.sample

    @property
    def ref(self):
        return self.table.ref

    @property
    def sect(self):
        return self.table.sect


class OneTableSeqGraph(OneTableGraph, OneSeqGraph, ABC):

    @classmethod
    def get_table_type(cls):
        return PosTableLoader

    @property
    def seq(self):
        return self.table.seq


class CartesianGraph(GraphBase, ABC):
    """ Graph with one pair of x and y axes. """

    @abstractmethod
    def get_xattr(self):
        """ Name of the x-axis attribute. """
        return ""

    @abstractmethod
    def get_yattr(self):
        """ Name of the y-axis attribute. """
        return ""

    def _layout_figure(self, fig: go.Figure):
        fig = super()._layout_figure(fig)
        fig.update_layout(xaxis=dict(title=self.get_xattr()),
                          yaxis=dict(title=self.get_yattr()))
        return fig


class GraphWriter(ABC):
    """ Write the proper types of graphs for a given table. """

    def __init__(self, table_file: Path):
        self.table_file = table_file

    @cached_property
    def table(self):
        """ The table providing the data for the graph(s). """
        return load(self.table_file)

    @abstractmethod
    def iter(self, *args, **kwargs):
        """ Yield every graph for the table. """
        yield GraphBase()

    def write(self, *args, csv: bool, html: bool, pdf: bool, **kwargs):
        """ Generate and write every graph for the table. """
        # Get the paths for every graph.
        paths = list()
        for graph in self.iter(*args, **kwargs):
            paths.extend(graph.write(csv=csv, html=html, pdf=pdf))
        return paths
