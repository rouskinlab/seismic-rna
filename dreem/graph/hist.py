import os
from abc import ABC, abstractmethod
from functools import cache
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import (find_tables, GraphWriter, CartesianGraph, OneTableGraph,
                   OneSampGraph, OneTableSectGraph)
from .color import RelColorMap
from ..core import docdef
from ..core.cli import (opt_table, opt_rels,
                        opt_yfrac, opt_hist_bins,
                        opt_csv, opt_html, opt_pdf,
                        opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..table.base import Table
from ..table.load import (RelReadTableLoader,
                          MaskReadTableLoader, ClusterReadTableLoader)

# Number of digits to which to round decimals.
PRECISION = 6

params = [
    opt_table,
    opt_rels,
    opt_hist_bins,
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
        fields: str,
        hist_bins: int,
        yfrac: bool, *,
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the graph read module. """
    writers = list(map(ReadHistogramWriter, find_tables(table)))
    return list(chain(*dispatch([writer.write for writer in writers],
                                max_procs, parallel, pass_n_procs=False,
                                kwargs=dict(fields=fields, count=yfrac,
                                            group=group, bins=hist_bins,
                                            csv=csv, html=html, pdf=pdf))))


class ReadHistogramWriter(GraphWriter):

    def iter(self, fields: str, bins: int, count: bool, group: bool):
        if isinstance(self.table, RelReadTableLoader):
            if group:
                yield RelReadHist(table=self.table, codes=fields,
                                  xfrac=bins, yfrac=count)
            else:
                for field in fields:
                    yield RelReadHist(table=self.table, codes=field,
                                      xfrac=bins, yfrac=count)
        elif isinstance(self.table, MaskReadTableLoader):
            if group:
                yield MaskReadHist(table=self.table, codes=fields,
                                   xfrac=bins, yfrac=count)
            else:
                for field in fields:
                    yield MaskReadHist(table=self.table, codes=field,
                                       xfrac=bins, yfrac=count)


class ReadHistogram(CartesianGraph, OneTableGraph, OneSampGraph, ABC):
    """ Histogram of an attribute of reads. """

    def __init__(self, *args,
                 table: (Table | RelReadTableLoader
                         | MaskReadTableLoader | ClusterReadTableLoader),
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.table = table

    @classmethod
    @abstractmethod
    def get_source(cls):
        """ Step from which the data came. """
        return ""

    @classmethod
    def get_table_type(cls):
        return RelReadTableLoader, MaskReadTableLoader, ClusterReadTableLoader

    @classmethod
    def get_data_type(cls):
        return pd.DataFrame

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap


class FieldReadHist(ReadHistogram, ABC):
    """ Read histogram of fields from a table. """

    def __init__(self, *args, codes: str, xfrac: int, yfrac: bool, **kwargs):
        super().__init__(*args, **kwargs)
        self.codes = codes
        self.xfrac = xfrac
        self.yfrac = yfrac

    def get_xattr(self):
        return "Fraction" if self.xfrac else "Count"

    def get_yattr(self):
        return "Read" + ("Fraction" if self.yfrac else "Count")

    @property
    def title(self):
        fields = '/'.join(sorted(Table.REL_CODES[c] for c in self.codes))
        return (f"{self.get_yattr()}s of {self.get_source()} {fields} bases "
                f"per read from sample {self.sample} over reference {self.ref}")

    @property
    def graph_filename(self):
        sort_codes = "".join(sorted(self.codes))
        fname = f"{self.get_source()}_{self.get_yattr()}_{sort_codes}_per-read"
        return fname.lower()

    def get_table_field(self, field_code: str):
        """ Load the data for one field from the table. """
        return (self.table.get_rel_frac(field_code).round(PRECISION)
                if self.xfrac else self.table.get_rel_count(field_code))

    @cache
    def _find_data_max(self):
        """ Find the maximum value in the data table. """
        return max(np.nanmax(self.get_table_field(code)) for code in self.codes)

    @cache
    def _get_bins(self):
        """ Get the bin edges for the histogram. """
        if self.xfrac < 0:
            raise ValueError(f"xfrac must be ≥ 0, but got {self.xfrac}")
        if self.xfrac == 0:
            # Each bin has width 1, with the maximum being the smallest
            # integer larger than every value in the data set.
            max_bin = int(self._find_data_max()) + 1
            n_bins = max_bin
        else:
            # The maximum value is 1, and the number of bins is xfrac.
            max_bin = 1
            n_bins = self.xfrac
        return np.linspace(0, max_bin, n_bins + 1)

    def _get_data(self):
        data = dict()
        bins = self._get_bins()
        for code in self.codes:
            # Compute the histogram for the field.
            series = self.get_table_field(code)
            data[series.name], _ = np.histogram(series, bins=bins)
        index = pd.MultiIndex.from_arrays([bins[:-1], bins[1:]],
                                          names=["[Lower Bound]",
                                                 "(Upper Bound)"])
        data = pd.DataFrame.from_dict(data)
        data.index = index
        return data

    def get_traces(self):
        traces = list()
        # Construct a trace for each field.
        for field, vals in self.data.items():
            # Define the text shown on hovering over a bar.
            hovertext = [(f"{round(lo, PRECISION)} ≤ {field} "
                          f"{self.get_xattr()} < {round(up, PRECISION)}: {y}")
                         for (lo, up), y in vals.items()]
            # Create a trace comprising all bars for this field.
            traces.append(go.Bar(name=field, x=vals.index.get_level_values(0),
                                 y=vals, marker_color=self.cmap[field],
                                 hovertext=hovertext, hoverinfo="text"))
        return traces


class SectReadHist(ReadHistogram, OneTableSectGraph, ABC):
    """ Bar graph of the positions in a section. """

    @property
    def title(self):
        return f"{super().title}, section {self.sect}"


class RelReadHist(FieldReadHist):
    """ Bar graph of related data from one sample at each position in a
    sequence. """

    @classmethod
    def get_source(cls):
        return "Related"


class MaskReadHist(FieldReadHist, SectReadHist):
    """ Bar graph of masked data from one sample at each position in a
    sequence. """

    @classmethod
    def get_source(cls):
        return "Masked"