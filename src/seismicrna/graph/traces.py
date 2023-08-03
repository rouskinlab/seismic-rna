from logging import getLogger

import pandas as pd
from plotly import graph_objects as go

from .base import PRECISION
from .color import ColorMap
from ..core.sect import BASE_NAME, POS_NAME
from ..core.seq import BASES

logger = getLogger(__name__)


def get_seq_base_scatter_trace(xdata: pd.Series, ydata: pd.Series,
                               cmap: ColorMap, base_int: int):
    # Validate the indexes.
    if not xdata.index.equals(ydata.index):
        raise ValueError("Indexes of x and y data must match, "
                         f"but got {xdata.index} and {ydata.index}")
    # Validate the base.
    base = chr(base_int)
    if base not in BASES.decode():
        raise ValueError(f"Invalid DNA base: '{base}'")
    # Find the position of every base of that type.
    seq_mask = xdata.index.get_level_values(BASE_NAME) == base
    # Get the values at those positions, excluding NaN values.
    xvals = xdata.loc[seq_mask].dropna()
    yvals = ydata.loc[seq_mask].dropna()
    # Join the x and y values into one DataFrame with only the positions
    # that are not missing in both series.
    xvals.name = 'x'
    yvals.name = 'y'
    vals = pd.concat([xvals, yvals], axis=1, join="inner")
    # Set the index of the values to the numerical positions.
    vals.index = vals.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{i}: {round(x, PRECISION), round(y, PRECISION)}"
                 for i, x, y in zip(vals.index, vals.x, vals.y)]
    # Create a trace comprising all bars for this base type.
    return go.Scatter(name=base, x=vals.x, y=vals.y,
                      mode="markers", marker_color=cmap[base_int],
                      hovertext=hovertext, hoverinfo="text")


def iter_seq_base_scatter_traces(xdata: pd.Series, ydata: pd.Series,
                                 cmap: ColorMap):
    for base in BASES:
        yield get_seq_base_scatter_trace(xdata, ydata, cmap, base)


def get_seq_base_bar_trace(data: pd.Series, cmap: ColorMap, base_int: int):
    # Validate the base.
    base = chr(base_int)
    if base not in BASES.decode():
        raise ValueError(f"Invalid DNA base: '{base}'")
    # Find the position of every base of that type.
    seq_mask = data.index.get_level_values(BASE_NAME) == base
    # Get the values at those positions, excluding NaN values.
    vals = data.loc[seq_mask].dropna()
    # Set the index of the values to the numerical positions.
    vals.index = vals.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x}: {round(y, PRECISION)}" for x, y in vals.items()]
    # Create a trace comprising all bars for this base type.
    return go.Bar(name=base, x=vals.index, y=vals,
                  marker_color=cmap[base_int],
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_seq_base_bar_traces(data: pd.Series, cmap: ColorMap):
    for base in BASES:
        yield get_seq_base_bar_trace(data, cmap, base)


def get_seq_stack_bar_trace(data: pd.Series, cmap: ColorMap):
    # Get the relationship from the name of the data series.
    rel = data.name
    # Get the sequence and positions.
    bases = data.index.get_level_values(BASE_NAME)
    pos = data.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x} {rel}: {round(y, PRECISION)}"
                 for base, x, y in zip(bases, pos, data, strict=True)]
    # Create a trace comprising all bars for this field.
    return go.Bar(name=rel, x=pos, y=data,
                  marker_color=cmap[rel],
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_seq_stack_bar_traces(data: pd.DataFrame, cmap: ColorMap):
    for rel, series in data.items():
        yield get_seq_stack_bar_trace(series, cmap)


def get_seq_line_trace(data: pd.Series, description: str):
    # Get the relationship from the name of the data series.
    rel = data.name
    # Get the positions.
    pos = data.index.get_level_values(POS_NAME)
    # Create a trace comprising all bars for this base type.
    return go.Scatter(name="correlation",
                      x=data.index.get_level_values(POS_NAME),
                      y=data)
                      #mode="lines",)
                      #hoverinfo="text")


def iter_seq_line_traces(data: pd.Series, description: str):
    yield get_seq_line_trace(data, description)
