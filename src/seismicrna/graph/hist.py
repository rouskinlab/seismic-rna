from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
import pandas as pd

from .color import ColorMapGraph, RelColorMap
from .onetable import (OneTableRelClusterGroupGraph,
                       OneTableRelClusterGroupRunner,
                       OneTableRelClusterGroupWriter)
from .rel import MultiRelsGraph
from .trace import (HIST_COUNT_NAME,
                    HIST_LOWER_NAME,
                    HIST_UPPER_NAME,
                    iter_hist_traces)
from ..core.arg import opt_hist_bins, opt_hist_margin
from ..core.header import parse_header


def get_edges_index(edges: np.ndarray, use_ratio: bool):
    """ Generate an index for the edges of histogram bins.

    Parameters
    ----------
    edges: numpy.ndarray
        Edges of histogram bins.
    use_ratio: bool
        Assume the edges represent ratios rather than counts.

    Returns
    -------
    pandas.Index | pandas.MultiIndex
        Index for the edges.
    """
    # Determine the lower and upper edge of each bin.
    lower = edges[:-1]
    upper = edges[1:]
    if use_ratio:
        # Make a MultiIndex of the lower and upper edge of each bin.
        return pd.MultiIndex.from_arrays([lower, upper],
                                         names=[HIST_LOWER_NAME,
                                                HIST_UPPER_NAME])
    # Make an Index of the count for each bin.
    if lower.size > 0:
        min_count = round(lower[0] + 0.5)
        counts = np.arange(min_count, min_count + lower.size)
        if not np.allclose(lower + 0.5, counts):
            raise ValueError("If the edges represent counts, then every lower "
                             "edge must be a half-integer and 1 more than the "
                             f"previous lower edge, but got {lower}")
    else:
        counts = np.arange(0)
    return pd.Index(counts, name=HIST_COUNT_NAME)


class HistogramGraph(OneTableRelClusterGroupGraph,
                     MultiRelsGraph,
                     ColorMapGraph,
                     ABC):
    """ Histogram of relationship(s) in one table. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    def __init__(self, *, hist_bins: int, hist_margin: float, **kwargs):
        super().__init__(**kwargs)
        if hist_bins <= 0:
            raise ValueError(
                f"hist_bins must be a positive integer, but got {hist_bins}"
            )
        self.num_bins = hist_bins
        if hist_margin < 0.:
            raise ValueError(f"hist_margin must be ≥ 0, but got {hist_margin}")
        self.margin = hist_margin

    @property
    def x_title(self):
        return self.data_kind.capitalize()

    @cached_property
    def data(self):
        # Fetch the raw data from the table.
        data = self._fetch_data(self.table,
                                k=self.k,
                                clust=self.clust)
        # Determine the edges of the bins.
        edges = self.get_edges(data)
        index = get_edges_index(edges, self.use_ratio)
        # Construct a histogram for each column.
        hist = dict()
        for col_name in data.columns:
            col_hist, _ = np.histogram(data[col_name], bins=edges)
            hist[col_name] = pd.Series(col_hist, index=index)
        return pd.DataFrame.from_dict(hist).reindex(columns=data.columns)

    @cached_property
    def data_header(self):
        """ Header of the selected data (not of the entire table). """
        return parse_header(self.data.columns)

    def get_bounds(self, data: pd.DataFrame):
        """ Get the lower and upper bounds of the histogram. """
        try:
            lo = float(np.nanmin(data))
        except ValueError:
            # If data contains no non-NaN values, then default to 0.
            lo = 0.
        else:
            # If the minimum is slightly greater than 0, then make it 0.
            if self.use_ratio and 0. < lo < self.margin:
                lo = 0.
        try:
            up = float(np.nanmax(data))
        except ValueError:
            # If data contains no non-NaN values, then default to 1.
            up = 1.
        else:
            # If the maximum is slightly less than 1, then make it 1.
            if self.use_ratio and 0. < 1. - up < self.margin:
                up = 1.
        return lo, up

    def get_edges(self, data: pd.DataFrame):
        """ Get the edges of the histogram bins. """
        # Make bins of floating-point values.
        lower, upper = self.get_bounds(data)
        if self.use_ratio:
            # Use the specified number of bins.
            num_bins = self.num_bins
        else:
            # Make one bin per integer, with half-integer bounds.
            lower = round(lower) - 0.5
            upper = round(upper) + 0.5
            num_bins = round(upper - lower)
        return np.linspace(lower, upper, num_bins + 1)

    def get_traces(self):
        for row, index in zip(range(1, self.nrows + 1),
                              self.data_header.iter_clust_indexes(),
                              strict=True):
            for trace in iter_hist_traces(self.data.loc[:, index], self.cmap):
                yield (row, 1), trace


class HistogramWriter(OneTableRelClusterGroupWriter, ABC):

    @classmethod
    @abstractmethod
    def get_graph_type(cls) -> type[HistogramGraph]:
        """ Type of graph. """

    def get_graph(self, rels_group: str, **kwargs):
        graph_type = self.get_graph_type()
        return graph_type(table=self.table, rels=rels_group, **kwargs)


class HistogramRunner(OneTableRelClusterGroupRunner, ABC):

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_hist_bins, opt_hist_margin]
