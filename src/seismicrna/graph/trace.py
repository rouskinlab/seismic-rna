from __future__ import annotations

from .color import ColorMap
from ..core.header import REL_NAME
from ..core.rna.roc import compute_auc
from ..core.seq.region import BASE_NAME, POS_NAME
from ..core.seq.xna import DNA

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


class _GraphObjects:
    """Lazy proxy for ``plotly.graph_objects``.

    Importing plotly is slow, and this module is imported whenever any graph
    subcommand is loaded (e.g. for ``--help``).  Deferring the plotly import
    until a trace is actually built keeps loading this module fast.
    """

    def __getattr__(self, name):
        from plotly import graph_objects as go

        return getattr(go, name)


go = _GraphObjects()

# Number of digits behind the decimal point to be kept.
PRECISION = 6
AUC_PRECISION = 3


def get_seq_base_scatter_trace(
    xdata: pd.Series, ydata: pd.Series, cmap: ColorMap, base: str
):
    """Build a scatter trace for one DNA base type.

    Parameters
    ----------
    xdata: pd.Series
        x-axis data indexed by a MultiIndex with levels including
        ``BASE_NAME`` and ``POS_NAME``.
    ydata: pd.Series
        y-axis data with the same index as ``xdata``.
    cmap: ColorMap
        Color map used to look up the color for ``base``.
    base: str
        Single-character DNA base code (A, C, G, T, or N).

    Returns
    -------
    plotly.graph_objects.Scatter
        A scatter trace for the positions of the given base type.
    """
    import pandas as pd

    # Validate the indexes.
    if not xdata.index.equals(ydata.index):
        raise ValueError(
            "Indexes of x and y data must match, "
            f"but got {xdata.index} and {ydata.index}"
        )
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
    hovertext = [
        f"{base}{i}: {round(x, PRECISION), round(y, PRECISION)}"
        for i, x, y in zip(vals.index, vals.x, vals.y)
    ]
    # Create a trace comprising all bars for this base type.
    return go.Scatter(
        name=base,
        x=vals.x,
        y=vals.y,
        mode="markers",
        marker_color=cmap.get(base),
        hovertext=hovertext,
        hoverinfo="text",
    )


def iter_seq_base_scatter_traces(xdata: pd.Series, ydata: pd.Series, cmap: ColorMap):
    """Yield one scatter trace per DNA base type.

    Parameters
    ----------
    xdata: pd.Series
        x-axis data indexed by a MultiIndex including ``BASE_NAME`` and
        ``POS_NAME``.
    ydata: pd.Series
        y-axis data with the same index as ``xdata``.
    cmap: ColorMap
        Color map for base coloring.

    Yields
    ------
    plotly.graph_objects.Scatter
        One trace per base in ``DNA.alph()``.
    """
    for base in DNA.alph():
        yield get_seq_base_scatter_trace(xdata, ydata, cmap, base)


def get_seq_base_bar_trace(data: pd.Series, cmap: ColorMap, base: str):
    """Build a bar trace for one DNA base type.

    Parameters
    ----------
    data: pd.Series
        Data indexed by a MultiIndex with levels including ``BASE_NAME``
        and ``POS_NAME``.
    cmap: ColorMap
        Color map used to look up the color for ``base``.
    base: str
        Single-character DNA base code (A, C, G, T, or N).

    Returns
    -------
    plotly.graph_objects.Bar
        A bar trace for the positions of the given base type.
    """
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
    return go.Bar(
        name=base,
        x=vals.index,
        y=vals,
        marker_color=cmap.get(base),
        hovertext=hovertext,
        hoverinfo="text",
    )


def iter_seq_base_bar_traces(data: pd.Series, cmap: ColorMap):
    for base in DNA.alph():
        yield get_seq_base_bar_trace(data, cmap, base)


def get_seq_stack_bar_trace(data: pd.Series, rel: str, cmap: ColorMap):
    """Build a stacked bar trace for one relationship type.

    Parameters
    ----------
    data: pd.Series
        Per-position data indexed by a MultiIndex with levels including
        ``BASE_NAME`` and ``POS_NAME``.
    rel: str
        Name of the relationship (used as the trace name and for color
        lookup).
    cmap: ColorMap
        Color map used to look up the color for ``rel``.

    Returns
    -------
    plotly.graph_objects.Bar
        A bar trace for all positions, colored by relationship type.
    """
    # Get the sequence and positions.
    bases = data.index.get_level_values(BASE_NAME)
    pos = data.index.get_level_values(POS_NAME)
    # Define the text shown on hovering over a bar.
    hovertext = [
        f"{base}{x} {rel}: {round(y, PRECISION)}"
        for base, x, y in zip(bases, pos, data, strict=True)
    ]
    # Create a trace comprising all bars for this field.
    return go.Bar(
        name=rel,
        x=pos,
        y=data,
        marker_color=cmap.get(rel),
        hovertext=hovertext,
        hoverinfo="text",
    )


def iter_seqbar_stack_traces(data: pd.DataFrame, cmap: ColorMap):
    for (_, series), rel in zip(
        data.items(), data.columns.get_level_values(REL_NAME), strict=True
    ):
        yield get_seq_stack_bar_trace(series, rel, cmap)


HIST_COUNT_NAME = "Count"
HIST_LOWER_NAME = "Lower"
HIST_UPPER_NAME = "Upper"


def get_hist_trace(data: pd.Series, rel: str, cmap: ColorMap):
    """Build a histogram bar trace for one relationship type.

    Parameters
    ----------
    data: pd.Series
        Histogram data.  If the index is a ``pd.MultiIndex`` its levels
        are ``HIST_LOWER_NAME`` and ``HIST_UPPER_NAME`` (ratio bins);
        otherwise it is a plain ``pd.Index`` named ``HIST_COUNT_NAME``
        (count bins).
    rel: str
        Name of the relationship (used as the trace name and for color
        lookup).
    cmap: ColorMap
        Color map used to look up the color for ``rel``.

    Returns
    -------
    plotly.graph_objects.Bar
        A bar trace representing the histogram for ``rel``.
    """
    import pandas as pd

    # Get the edges of the bins.
    if isinstance(data.index, pd.MultiIndex):
        lower = data.index.get_level_values(HIST_LOWER_NAME)
        upper = data.index.get_level_values(HIST_UPPER_NAME)
        center = (lower + upper) / 2.0
        hovertext = [
            (
                f"[{round(lo, PRECISION)} - {round(up, PRECISION)}] "
                f"{rel}: {round(value, PRECISION)}"
            )
            for lo, up, value in zip(lower, upper, data, strict=True)
        ]
    else:
        if data.index.name != HIST_COUNT_NAME:
            raise ValueError(
                f"Expected index to be named {repr(HIST_COUNT_NAME)}, "
                f"but got {repr(data.index.name)}"
            )
        center = data.index.values
        hovertext = [
            f"{count} {rel}: {round(value, PRECISION)}" for count, value in data.items()
        ]
    # Create a trace comprising all bars for this field.
    return go.Bar(
        name=rel,
        x=center,
        y=data,
        marker_color=cmap.get(rel),
        hovertext=hovertext,
        hoverinfo="text",
    )


def iter_hist_traces(data: pd.DataFrame, cmap: ColorMap):
    for (_, series), rel in zip(
        data.items(), data.columns.get_level_values(REL_NAME), strict=True
    ):
        yield get_hist_trace(series, rel, cmap)


def get_seq_line_trace(data: pd.Series):
    return go.Scatter(x=data.index.get_level_values(POS_NAME), y=data)


def iter_seq_line_traces(data: pd.Series, *_, **__):
    """Yield a single line trace for sequence data.

    Parameters
    ----------
    data: pd.Series
        Per-position data indexed by a MultiIndex including
        ``POS_NAME``.
    *_, **__
        Ignored extra arguments (for call-signature compatibility).

    Yields
    ------
    plotly.graph_objects.Scatter
        One line trace for the entire series.
    """
    yield get_seq_line_trace(data)


def _format_profile_struct(profile: str, struct: str, auc: float | None = None):
    """Format a trace name combining profile, structure, and AUC.

    Parameters
    ----------
    profile: str
        Name of the mutational profile.
    struct: str
        Name or identifier of the RNA structure.
    auc: float or None, optional
        Area under the ROC curve; if provided it is appended to the
        label.

    Returns
    -------
    str
        Formatted label string.
    """
    text = f"{profile}, {struct}"
    if auc is not None:
        text = f"{text} (AUC = {round(auc, AUC_PRECISION)})"
    return text


def get_roc_trace(fpr: pd.Series, tpr: pd.Series, profile: str, struct: str):
    """Build an ROC curve trace for one profile/structure combination.

    Parameters
    ----------
    fpr: pd.Series
        False positive rates.
    tpr: pd.Series
        True positive rates.
    profile: str
        Name of the mutational profile.
    struct: str
        Name or identifier of the RNA structure.

    Returns
    -------
    plotly.graph_objects.Scatter
        An ROC curve trace with AUC annotated in the trace name.
    """
    name = _format_profile_struct(profile, struct, compute_auc(fpr.values, tpr.values))
    return go.Scatter(x=fpr, y=tpr, name=name)


def iter_roc_traces(fprs: pd.DataFrame, tprs: pd.DataFrame, profile: str):
    """Yield one ROC trace per RNA structure.

    Parameters
    ----------
    fprs: pd.DataFrame
        DataFrame whose columns correspond to structures and whose
        values are false positive rates.
    tprs: pd.DataFrame
        DataFrame whose columns correspond to structures and whose
        values are true positive rates; columns must match ``fprs``.
    profile: str
        Name of the mutational profile.

    Yields
    ------
    plotly.graph_objects.Scatter
        One ROC trace per structure column.
    """
    for (sf, fpr), (st, tpr) in zip(fprs.items(), tprs.items(), strict=True):
        if sf != st:
            raise ValueError(f"Structure names differ: {repr(sf)} ≠ {repr(st)}")
        yield get_roc_trace(fpr, tpr, profile, str(sf))


def get_rolling_auc_trace(auc: pd.Series, profile: str, struct: str):
    return go.Scatter(
        x=auc.index.get_level_values(POS_NAME),
        y=auc,
        name=_format_profile_struct(profile, struct),
    )


def iter_rolling_auc_traces(aucs: pd.DataFrame, profile: str):
    for struct, auc in aucs.items():
        yield get_rolling_auc_trace(auc, profile, str(struct))


def get_line_trace(data: pd.Series, cluster: str):
    return go.Scatter(x=data.index.get_level_values(POS_NAME), y=data, name=cluster)


def iter_line_traces(lines: pd.DataFrame):
    for cluster, line in lines.items():
        yield get_line_trace(line, str(cluster))


def get_pairwise_position_trace(data: pd.Series, end5: int, end3: int):
    """Build a pairwise-position heatmap trace.

    Parameters
    ----------
    data: pd.Series
        Long-form data with a two-level ``pd.MultiIndex`` of
        ``(position_a, position_b)`` pairs and values representing the
        pairwise statistic.
    end5: int
        5' end position of the region (inclusive, 1-based).
    end3: int
        3' end position of the region (inclusive, 1-based).

    Returns
    -------
    plotly.graph_objects.Heatmap
        A symmetric heatmap trace covering positions ``end5`` to
        ``end3``.
    """
    import numpy as np
    import pandas as pd

    # The data must be a long-form Series with a two-level MultiIndex.
    # Convert the data to wide-form and make them symmetric.
    if not isinstance(data, pd.Series):
        raise TypeError(f"data must be a Series, but got {type(data).__name__}")
    if not isinstance(data.index, pd.MultiIndex):
        raise TypeError(
            f"data.index must be a MultiIndex, but got {type(data.index).__name__}"
        )
    if data.index.nlevels != 2:
        raise ValueError(f"data.index must have 2 levels, but got {data.index.nlevels}")
    matrix_index = pd.RangeIndex(end5, end3 + 1)
    matrix = pd.DataFrame(np.nan, matrix_index, matrix_index)
    for (pos_x, pos_y), value in data.items():
        matrix.at[(pos_x, pos_y)] = value
        matrix.at[(pos_y, pos_x)] = value
    return go.Heatmap(
        x=matrix_index,
        y=matrix_index,
        z=matrix,
        hoverongaps=False,
        colorscale="rdbu_r",
        zmid=0,
    )


def iter_stack_bar_traces(data: pd.DataFrame):
    for column_label, column in reversed(list(data.items())):
        yield go.Bar(
            name=f"{data.columns.name} {column_label}",
            x=data.index,
            y=column,
            hovertext=[
                (
                    f"{data.index.name} {index_label}, "
                    f"{data.columns.name} {column_label}: "
                    f"{round(value, PRECISION)}"
                )
                for index_label, value in column.items()
            ],
            hoverinfo="text",
        )
