import numpy as np
import pandas as pd
from plotly import graph_objects as go

from .color import ColorMap
from ..core.header import REL_NAME
from ..core.rna import compute_auc
from ..core.seq import BASE_NAME, POS_NAME, DNA

# Number of digits behind the decimal point to be kept.
PRECISION = 6
AUC_PRECISION = 3


def get_seq_base_scatter_trace(xdata: pd.Series,
                               ydata: pd.Series,
                               cmap: ColorMap,
                               base: str):
    # Validate the indexes.
    if not xdata.index.equals(ydata.index):
        raise ValueError("Indexes of x and y data must match, "
                         f"but got {xdata.index} and {ydata.index}")
    # Validate the base.
    if base not in DNA.alph():
        raise ValueError(f"Invalid DNA base: {repr(base)}")
    # Find the position of every base of that type.
    seq_mask = xdata.index.get_level_values(BASE_NAME) == base
    # Get the values at those positions, excluding NaN values.
    xvals = xdata.loc[seq_mask].dropna()
    yvals = ydata.loc[seq_mask].dropna()
    # Join the x and y values into one DataFrame with only the positions
    # that are not missing in both series.
    xvals.name = "x"
    yvals.name = "y"
    vals = pd.concat([xvals, yvals], axis=1, join="inner")
    # Set the index of the values to the numerical positions.
    vals.index = vals.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{i}: {round(x, PRECISION), round(y, PRECISION)}"
                 for i, x, y in zip(vals.index, vals.x, vals.y)]
    # Create a trace comprising all bars for this base type.
    return go.Scatter(name=base,
                      x=vals.x,
                      y=vals.y,
                      mode="markers",
                      marker_color=cmap.get(base),
                      hovertext=hovertext,
                      hoverinfo="text")


def iter_seq_base_scatter_traces(xdata: pd.Series,
                                 ydata: pd.Series,
                                 cmap: ColorMap):
    for base in DNA.alph():
        yield get_seq_base_scatter_trace(xdata, ydata, cmap, base)


def get_seq_base_bar_trace(data: pd.Series, cmap: ColorMap, base: str):
    # Validate the base.
    if base not in DNA.alph():
        raise ValueError(f"Invalid DNA base: {repr(base)}")
    # Find the position of every base of that type.
    seq_mask = data.index.get_level_values(BASE_NAME) == base
    # Get the values at those positions, excluding NaN values.
    vals = data.loc[seq_mask].dropna()
    # Set the index of the values to the numerical positions.
    vals.index = vals.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x}: {round(y, PRECISION)}" for x, y in vals.items()]
    # Create a trace comprising all bars for this base type.
    return go.Bar(name=base,
                  x=vals.index,
                  y=vals,
                  marker_color=cmap.get(base),
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_seq_base_bar_traces(data: pd.Series, cmap: ColorMap):
    for base in DNA.alph():
        yield get_seq_base_bar_trace(data, cmap, base)


def get_seq_stack_bar_trace(data: pd.Series, rel: str, cmap: ColorMap):
    # Get the sequence and positions.
    bases = data.index.get_level_values(BASE_NAME)
    pos = data.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [f"{base}{x} {rel}: {round(y, PRECISION)}"
                 for base, x, y in zip(bases, pos, data, strict=True)]
    # Create a trace comprising all bars for this field.
    return go.Bar(name=rel,
                  x=pos,
                  y=data,
                  marker_color=cmap.get(rel),
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_seqbar_stack_traces(data: pd.DataFrame, cmap: ColorMap):
    for (_, series), rel in zip(data.items(),
                                data.columns.get_level_values(REL_NAME),
                                strict=True):
        yield get_seq_stack_bar_trace(series, rel, cmap)


HIST_COUNT_NAME = "Count"
HIST_LOWER_NAME = "Lower"
HIST_UPPER_NAME = "Upper"


def get_hist_trace(data: pd.Series, rel: str, cmap: ColorMap):
    # Get the edges of the bins.
    if isinstance(data.index, pd.MultiIndex):
        lower = data.index.get_level_values(HIST_LOWER_NAME)
        upper = data.index.get_level_values(HIST_UPPER_NAME)
        center = (lower + upper) / 2.
        hovertext = [(f"[{round(lo, PRECISION)} - {round(up, PRECISION)}] "
                      f"{rel}: {round(value, PRECISION)}")
                     for lo, up, value in zip(lower, upper, data, strict=True)]
    else:
        if data.index.name != HIST_COUNT_NAME:
            raise ValueError(
                f"Expected index to be named {repr(HIST_COUNT_NAME)}, "
                f"but got {repr(data.index.name)}"
            )
        center = data.index.values
        hovertext = [f"{count} {rel}: {round(value, PRECISION)}"
                     for count, value in data.items()]
    # Create a trace comprising all bars for this field.
    return go.Bar(name=rel,
                  x=center,
                  y=data,
                  marker_color=cmap.get(rel),
                  hovertext=hovertext,
                  hoverinfo="text")


def iter_hist_traces(data: pd.DataFrame, cmap: ColorMap):
    for (_, series), rel in zip(data.items(),
                                data.columns.get_level_values(REL_NAME),
                                strict=True):
        yield get_hist_trace(series, rel, cmap)


def get_seq_line_trace(data: pd.Series):
    return go.Scatter(x=data.index.get_level_values(POS_NAME),
                      y=data)


def iter_seq_line_traces(data: pd.Series, *_, **__):
    yield get_seq_line_trace(data)


def _format_profile_struct(profile: str, struct: str, auc: float | None = None):
    text = f"{profile}, {struct}"
    if auc is not None:
        text = f"{text} (AUC = {round(auc, AUC_PRECISION)})"
    return text


def get_roc_trace(fpr: pd.Series, tpr: pd.Series, profile: str, struct: str):
    name = _format_profile_struct(profile,
                                  struct,
                                  compute_auc(fpr.values, tpr.values))
    return go.Scatter(x=fpr, y=tpr, name=name)


def iter_roc_traces(fprs: pd.DataFrame, tprs: pd.DataFrame, profile: str):
    for (sf, fpr), (st, tpr) in zip(fprs.items(), tprs.items(), strict=True):
        if sf != st:
            raise ValueError(f"Structure names differ: {repr(sf)} ≠ {repr(st)}")
        yield get_roc_trace(fpr, tpr, profile, str(sf))


def get_rolling_auc_trace(auc: pd.Series, profile: str, struct: str):
    return go.Scatter(x=auc.index.get_level_values(POS_NAME),
                      y=auc,
                      name=_format_profile_struct(profile, struct))


def iter_rolling_auc_traces(aucs: pd.DataFrame, profile: str):
    for struct, auc in aucs.items():
        yield get_rolling_auc_trace(auc, profile, str(struct))


def get_line_trace(gini: pd.Series, cluster: str):
    return go.Scatter(x=gini.index.get_level_values(POS_NAME),
                      y=gini,
                      name=cluster)


def iter_line_traces(lines: pd.DataFrame):
    for cluster, line in lines.items():
        yield get_line_trace(line, str(cluster))


def get_pairwise_position_trace(data: pd.Series, end5: int, end3: int):
    # The data must be a long-form Series with a two-level MultiIndex.
    # Convert the data to wide-form and make them symmetric.
    if not isinstance(data, pd.Series):
        raise TypeError("data must be a Series, "
                        f"but got {type(data).__name__}")
    if not isinstance(data.index, pd.MultiIndex):
        raise TypeError("data.index must be a MultiIndex, "
                        f"but got {type(data.index).__name__}")
    if data.index.nlevels != 2:
        raise ValueError("data.index must have 2 levels, "
                         f"but got {data.index.nlevels}")
    matrix_index = pd.RangeIndex(end5, end3 + 1)
    matrix = pd.DataFrame(np.nan, matrix_index, matrix_index)
    for (pos_x, pos_y), value in data.items():
        matrix.at[(pos_x, pos_y)] = value
        matrix.at[(pos_y, pos_x)] = value
    return go.Heatmap(x=matrix_index,
                      y=matrix_index,
                      z=matrix,
                      hoverongaps=False,
                      colorscale="rdbu_r",
                      zmid=0)


def iter_stack_bar_traces(data: pd.DataFrame):
    for column_label, column in data.items():
        yield go.Bar(name=f"{data.columns.name} {column_label}",
                     x=data.index,
                     y=column,
                     hovertext=[(f"{data.index.name} {index_label}, "
                                 f"{data.columns.name} {column_label}: "
                                 f"{round(value, PRECISION)}")
                                for index_label, value in column.items()],
                     hoverinfo="text")
