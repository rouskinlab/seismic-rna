from abc import ABC, abstractmethod
from functools import cache, cached_property
from logging import getLogger
from pathlib import Path
from typing import Any, Callable, Iterable

from plotly import graph_objects as go
from plotly.subplots import make_subplots

from .color import ColorMap, get_cmap
from ..core import path
from ..core.seq import DNA
from ..table.base import Table
from ..table.load import TableLoader, PosTableLoader

logger = getLogger(__name__)

# Number of digits behind the decimal point to be kept.
PRECISION = 6


def _write_graph(writer: Callable[[Path], Any], file: Path):
    """ Write an image or raw data for a graph to a file. """
    if file.is_file():
        logger.warning(f"File exists: {file}")
    else:
        writer(file)
    return file


class GraphBase(ABC):

    def __init__(self, cmap: str | None = None):
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
    def out_dir(self) -> Path:
        """ Output directory. """

    @property
    @abstractmethod
    def graph_filename(self):
        """ Name of the graph file. """

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

    @property
    @abstractmethod
    def nrows(self) -> int:
        """ Number of rows of subplots. """

    @property
    @abstractmethod
    def ncols(self) -> int:
        """ Number of columns of subplots. """

    @property
    def row_titles(self) -> list[str] | None:
        """ Title of each row. """
        return None

    @property
    def col_titles(self) -> list[str] | None:
        """ Title of each column. """
        return None

    @property
    def subplot_titles(self) -> list[str] | None:
        """ Titles of the subplots. """
        return None

    def _get_subplots_params(self):
        return dict(rows=self.nrows,
                    cols=self.ncols,
                    row_titles=self.row_titles,
                    column_titles=self.col_titles)

    def _figure_init(self):
        """ Initialize the figure. """
        return make_subplots(**self._get_subplots_params())

    def _figure_data(self, fig: go.Figure):
        """ Add data to the figure. """
        for (row, col), trace in self.get_traces():
            fig.add_trace(trace, row=row, col=col)

    def _figure_layout(self, fig: go.Figure):
        """ Update the figure's layout. """
        fig.update_layout(title=self.title,
                          plot_bgcolor="#ffffff",
                          paper_bgcolor="#ffffff")

    @cache
    def get_figure(self):
        """ Figure object. """
        fig = self._figure_init()
        self._figure_data(fig)
        self._figure_layout(fig)
        return fig

    def write_csv(self):
        """ Write the graph's source data to a CSV file. """
        return _write_graph(self.data.to_csv,
                            self.get_path(ext=path.CSV_EXT))

    def write_html(self):
        """ Write the graph to an HTML file. """
        return _write_graph(self.get_figure().write_html,
                            self.get_path(ext=path.HTML_EXT))

    def write_pdf(self):
        """ Write the graph to a PDF file. """
        return _write_graph(self.get_figure().write_image,
                            self.get_path(ext=path.PDF_EXT))

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
        return {**super().get_path_fields(), path.SAMP: self.sample}


class OneSampleGraph(SampleGraph, ABC):
    """ Graph of one sample. """


class TwoSampleGraph(SampleGraph, ABC):
    """ Graph of two samples. """

    @property
    @abstractmethod
    def sample1(self):
        """ Name of the first sample. """
        return ""

    @property
    @abstractmethod
    def sample2(self):
        """ Name of the second sample. """
        return ""

    @property
    def sample(self):
        return f"{self.sample1}_vs_{self.sample2}"


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


class OneTableGraph(OneSampleGraph, OneRefGraph, ABC):
    """ Graph of data from one TableLoader. """

    def __init__(self, *args, table: TableLoader, **kwargs):
        super().__init__(*args, **kwargs)
        if not isinstance(table, self.get_table_type()):
            raise TypeError(f"{self.__class__.__name__} expected table "
                            f"of type '{self.get_table_type().__name__}', "
                            f"but got type '{type(table).__name__}'")
        self._table = table

    @classmethod
    @abstractmethod
    def get_table_type(cls):
        """ Type of TableLoader for this graph. """
        return TableLoader

    @property
    def table(self) -> (TableLoader | PosTableLoader | Table):
        """ Table of data. """
        return self._table

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


class TwoTableGraph(TwoSampleGraph, OneRefGraph, ABC):

    def __init__(self, *args, table1: TableLoader, table2: TableLoader,
                 **kwargs):
        super().__init__(*args, **kwargs)
        for table, table_type in zip([table1,
                                      table2],
                                     [self.get_table1_type(),
                                      self.get_table2_type()]):
            if not isinstance(table, table_type):
                raise TypeError(f"{self.__class__.__name__} expected table "
                                f"of type '{table_type.__name__}', "
                                f"but got type '{type(table).__name__}'")
        self._table1 = table1
        self._table2 = table2

    @classmethod
    @abstractmethod
    def get_table1_type(cls):
        """ Type of the first TableLoader for this graph. """
        return TableLoader

    @classmethod
    @abstractmethod
    def get_table2_type(cls):
        """ Type of the second TableLoader for this graph. """
        return TableLoader

    @property
    def table1(self) -> (TableLoader | PosTableLoader | Table):
        """ First table of data. """
        return self._table1

    @property
    def table2(self) -> (TableLoader | PosTableLoader | Table):
        """ Second table of data. """
        return self._table2

    def _get_common_attribute(self, name: str):
        """ Get the common attribute for tables 1 and 2. """
        attr1 = self.table1.__getattribute__(name)
        attr2 = self.table2.__getattribute__(name)
        if attr1 != attr2:
            raise ValueError(f"Attribute '{name}' differs between tables 1 "
                             f"({repr(attr1)}) and 2 ({repr(attr2)})")
        return attr1

    @property
    def out_dir(self):
        return self._get_common_attribute("out_dir")

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
        if self.sample1 == self.sample2:
            return self.sample1
        return f"{self.sample1}__and__{self.sample2}"


class OneTableSeqGraph(OneTableGraph, OneSeqGraph, ABC):

    @classmethod
    def get_table_type(cls):
        return PosTableLoader

    @property
    def seq(self):
        return self.table.seq


class TwoTableSeqGraph(TwoTableGraph, OneSeqGraph, ABC):

    @classmethod
    def get_table1_type(cls):
        return PosTableLoader

    @classmethod
    def get_table2_type(cls):
        return PosTableLoader

    @property
    def seq(self):
        return self._get_common_attribute("seq")


class CartesianGraph(GraphBase, ABC):
    """ Graph with one pair of x and y axes. """

    @property
    @abstractmethod
    def x_title(self) -> str:
        """ Title of the x-axis. """

    @property
    @abstractmethod
    def y_title(self) -> str:
        """ Title of the y-axis. """

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_xaxes(linewidth=1,
                         linecolor="#000000",
                         autorange=True)
        fig.update_yaxes(linewidth=1,
                         linecolor="#000000",
                         autorange=True)

    def _get_subplots_params(self):
        return {**super()._get_subplots_params(),
                **dict(shared_xaxes="all", shared_yaxes="all",
                       x_title=self.x_title, y_title=self.y_title)}
