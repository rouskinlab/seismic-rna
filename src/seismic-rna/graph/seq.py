from abc import ABC, abstractmethod
from functools import cache
from itertools import chain
from logging import getLogger
import os
from pathlib import Path

from click import command
import pandas as pd
from plotly import graph_objects as go

from .base import (PRECISION, find_tables, GraphWriter, CartesianGraph,
                   OneTableSeqGraph, OneSampGraph)
from .color import RelColorMap, SeqColorMap
from ..core import docdef
from ..core.cli import (opt_table, opt_table_cols, opt_stack, opt_yfrac,
                        opt_csv, opt_html, opt_pdf, opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.sect import BASE_NAME, POS_NAME
from ..core.seq import BASES
from ..core.types import get_subclasses_module
from ..table.base import RelTypeTable, REL_CODES
from ..table.load import (TableLoader, RelPosTableLoader,
                          MaskPosTableLoader, ClustPosTableLoader)

logger = getLogger(__name__)

# Number of digits to which to round decimals.

params = [
    opt_table,
    opt_table_cols,
    opt_stack,
    opt_yfrac,
    opt_csv,
    opt_html,
    opt_pdf,
    opt_max_procs,
    opt_parallel,
]


@command(__name__.split(os.path.extsep)[-1], params=params)
def cli(*args, **kwargs):
    """ Create bar graphs of positional attributes. """
    return run(*args, **kwargs)


@docdef.auto()
def run(table: tuple[str, ...],
        table_cols: str,
        stack: bool,
        yfrac: bool, *,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph pos module. """
    writers = list(map(SeqGraphWriter, find_tables(table)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(rels=table_cols, stack=stack,
                                            yfrac=yfrac, csv=csv,
                                            html=html, pdf=pdf))))


class SeqGraphWriter(GraphWriter):

    def iter(self, rels: str, stack: str, yfrac: bool):
        for rel in rels:
            for graph_type in get_subclasses_module(SerialSeqGraph, __name__):
                if type(self.table) in graph_type.sources():
                    yield graph_type(table=self.table, codes=rel, yfrac=yfrac)
        if stack:
            for graph_type in get_subclasses_module(StackedSeqGraph, __name__):
                if type(self.table) in graph_type.sources():
                    yield graph_type(table=self.table, codes=stack, yfrac=yfrac)


# Base Sequence Graphs #################################################

class SeqGraph(CartesianGraph, OneTableSeqGraph, OneSampGraph, ABC):
    """ Bar graph wherein each bar represents one sequence position. """

    def __init__(self, *args,
                 table: (RelPosTableLoader
                         | MaskPosTableLoader
                         | ClustPosTableLoader),
                 codes: str,
                 yfrac: bool,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table
        self.codes = codes
        self.yfrac = yfrac

    @classmethod
    @abstractmethod
    def sources(cls) -> dict[type[TableLoader], str]:
        """ Names of the sources of data. """
        return dict()

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

    def get_xattr(self):
        return POS_NAME

    def get_yattr(self):
        return "Fraction" if self.yfrac else "Count"

    @property
    def title(self):
        fields = '/'.join(sorted(REL_CODES[c] for c in self.codes))
        return (f"{self.get_yattr()}s of {fields} bases in {self.source} reads "
                f"from {self.sample} per position in {self.ref}:{self.sect}")

    @property
    @abstractmethod
    def sort_codes(self):
        return "".join(sorted(self.codes))

    @property
    def graph_filename(self):
        return f"{self.source}_{self.sort_codes}_{self.get_yattr()}".lower()

    def get_table_field(self, table: RelTypeTable | TableLoader, code: str):
        return (table.fract_rel(code).round(PRECISION) if self.yfrac
                else table.count_rel(code))


# Sequence Graphs by Source ############################################

class PopAvgSeqGraph(SeqGraph, ABC):

    @classmethod
    def sources(cls) -> dict[type[TableLoader], str]:
        return {RelPosTableLoader: "Related",
                MaskPosTableLoader: "Masked"}


class ClustSeqGraph(SeqGraph, ABC):

    @classmethod
    def sources(cls) -> dict[type[TableLoader], str]:
        return {ClustPosTableLoader: "Clustered"}

    @classmethod
    def get_data_type(cls):
        return pd.DataFrame


# Sequence Graphs by Series Type #######################################

class SerialSeqGraph(SeqGraph, ABC):
    """ Bar graph with a single series of data. """

    def get_traces(self):
        traces = list()
        # Construct a trace for each type of base.
        for base in BASES:
            letter = chr(base)
            # Find the position of every base of that type.
            seq_mask = self.data.index.get_level_values(BASE_NAME) == letter
            # Get the values at those positions, excluding NaN values.
            vals = self.data.loc[seq_mask].dropna()
            # Set the index of the values to the numerical positions.
            vals.index = vals.index.get_level_values(POS_NAME)
            # Check if there are any values to graph.
            if vals.size > 0:
                # Define the text shown on hovering over a bar.
                hovertext = [f"{letter}{x}: {y}" for x, y in vals.items()]
                # Create a trace comprising all bars for this base type.
                traces.append(go.Bar(name=letter, x=vals.index, y=vals,
                                     marker_color=self.cmap[base],
                                     hovertext=hovertext,
                                     hoverinfo="text"))
        return traces

    @property
    def code(self):
        """ The code of the field to graph. """
        if len(self.codes) != 1:
            raise ValueError(f"Expected 1 code, but got {len(self.codes)}")
        return self.codes[0]

    @property
    def sort_codes(self):
        return "-".join(["serial", super().sort_codes])

    def _get_data(self) -> pd.Series:
        return self.get_table_field(self.table, self.code)


class StackedSeqGraph(SeqGraph, ABC):
    """ Stacked bar graph wherein each stacked bar represents multiple
    outcomes for a base in a sequence. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    @property
    def sort_codes(self):
        return "-".join(["stacked", super().sort_codes])

    def _get_data(self):
        data = dict()
        for code in self.codes:
            series = self.get_table_field(self.table, code)
            data[series.name] = series
        return pd.DataFrame.from_dict(data)

    def get_traces(self):
        traces = list()
        # Construct a trace for each field.
        for field, vals in self.data.items():
            # Get the sequence and positions.
            bases = vals.index.get_level_values(BASE_NAME)
            pos = vals.index.get_level_values(POS_NAME)
            # Define the text shown on hovering over a bar.
            hovertext = [f"{base}{x} {field}: {y}"
                         for base, x, y in zip(bases, pos, vals, strict=True)]
            # Create a trace comprising all bars for this field.
            traces.append(go.Bar(name=field, x=pos, y=vals,
                                 marker_color=self.cmap[field],
                                 hovertext=hovertext,
                                 hoverinfo="text"))
        return traces

    @cache
    def get_figure(self):
        fig = super().get_figure()
        # Stack the bars at each position.
        fig.update_layout(barmode="stack")
        return fig


# Instantiable Sequence Graphs #########################################

class PopAvgSerialSeqGraph(PopAvgSeqGraph, SerialSeqGraph):

    @classmethod
    def get_data_type(cls):
        return pd.Series


class PopAvgStackedSeqGraph(PopAvgSeqGraph, StackedSeqGraph):

    @classmethod
    def get_data_type(cls):
        return pd.DataFrame
