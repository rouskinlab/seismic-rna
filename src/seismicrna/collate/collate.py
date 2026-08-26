import json
import os
import re

from functools import cache
from hashlib import sha1
from pathlib import Path
from jinja2 import Template
from typing import Iterable

from ..core import path
from ..core.logs import logger
from ..logo import draw_seismic_logo

ASSETS_DIR = Path(__file__).parent / "assets"

# Uniform tile height (in GridStack rows) for groups with more than one graph.
MULTI_TILE_ROWS = 6
# Default rows for a lone graph; the report JS grows it to fill the window.
SOLO_TILE_ROWS = 7


def _grid_width(n: int, svg_only: bool = False) -> int:
    """Uniform default tile width (in 12-column units) for a group of ``n``
    graphs, chosen so the tiles align in clean, evenly sized rows.

    Graphs default to at most two per row, because plots become too thin at
    three or more across (users can still make them narrower by hand). Groups
    that contain only SVG structure drawings are exempt: those stay legible
    when packed more densely, so they may go up to three per row.
    """
    if n <= 1:
        return 12
    if svg_only:
        if n == 2:
            return 6  # 2 across
        return 4  # 3 across
    return 6  # graphs: at most 2 across


def _place_band(entries: list, width: int, y0: int, complete_last: bool):
    """Lay a run of same-kind tiles into rows of a fixed width starting at row
    ``y0``. When ``complete_last`` is set, a partial final row is widened to
    span the full 12 columns so nothing below can rise into the leftover gap
    (which would produce a ragged mixed row)."""
    cols = max(1, 12 // width)
    n = len(entries)
    remainder = n % cols
    for i, entry in enumerate(entries):
        entry["h"] = MULTI_TILE_ROWS
        row, col = divmod(i, cols)
        entry["y"] = y0 + row * MULTI_TILE_ROWS
        in_last_partial = complete_last and remainder and i >= n - remainder
        if in_last_partial:
            fill = 12 // remainder
            entry["w"] = fill
            entry["x"] = (i - (n - remainder)) * fill
        else:
            entry["w"] = width
            entry["x"] = col * width
    return (n + cols - 1) // cols  # number of rows used


def _layout_group(files_data: list):
    """Assign default position and size to every tile in a group.

    A lone graph fills the width. A single-kind group tiles uniformly. A mixed
    group places the plots (at most two per row) in a top band and the SVG
    structure drawings (up to three per row) in a band beneath it, so the two
    kinds keep kind-appropriate widths without interleaving into ragged rows.
    """
    n = len(files_data)
    if n == 1:
        files_data[0].update(w=12, h=SOLO_TILE_ROWS, x=0, y=0)
        return
    svgs = [e for e in files_data if e["kind"] == "svg"]
    plots = [e for e in files_data if e["kind"] != "svg"]
    if svgs and plots:
        rows = _place_band(plots, 6, 0, complete_last=True)
        _place_band(
            svgs,
            _grid_width(len(svgs), svg_only=True),
            rows * MULTI_TILE_ROWS,
            complete_last=False,
        )
    elif svgs:
        _place_band(svgs, _grid_width(n, svg_only=True), 0, complete_last=False)
    else:
        _place_band(plots, _grid_width(n, svg_only=False), 0, complete_last=False)


@cache
def _load_asset(name: str) -> str:
    """Read a vendored/report asset file bundled with this package."""
    return (ASSETS_DIR / name).read_text(encoding="utf-8")


# Regexes for pulling apart a Plotly-generated standalone HTML file so that a
# single copy of the plotly.js library can be shared across many graphs.
_PLOTLY_DIV_RE = re.compile(
    r'<div id="(?P<id>[^"]+)" class="plotly-graph-div"[^>]*>\s*</div>'
)
_SCRIPT_RE = re.compile(
    r'<script type="text/javascript"[^>]*>(.*?)</script>', re.DOTALL
)
_PLOTLY_VERSION_RE = re.compile(r"plotly\.js v([\d.]+)")


def parse_plotly_html(html: str) -> dict | None:
    """Extract the reusable pieces of a Plotly standalone HTML document.

    Returns a dict with the plot ``<div>``, its id, the ``Plotly.newPlot``
    data script, the plotly.js library source, and the library version; or
    ``None`` if the document does not look like Plotly output (e.g. an
    arbitrary HTML graph), in which case the caller should embed it verbatim.
    """
    div_match = _PLOTLY_DIV_RE.search(html)
    if div_match is None:
        return None
    lib_js = config_js = data_js = version = None
    for script in _SCRIPT_RE.finditer(html):
        body = script.group(1)
        if "plotly.js v" in body:
            lib_js = body
            version_match = _PLOTLY_VERSION_RE.search(body)
            version = version_match.group(1) if version_match else None
        elif "PlotlyConfig" in body:
            config_js = body
        elif "Plotly.newPlot" in body:
            data_js = body
    if lib_js is None or data_js is None or version is None:
        return None
    return {
        "div_id": div_match.group("id"),
        "div_html": div_match.group(0),
        "lib_js": lib_js,
        "config_js": config_js or "",
        "data_js": data_js,
        "version": version,
    }


def _slug(text: str, fallback: str) -> str:
    """Make a DOM-id-safe slug from arbitrary text."""
    slug = "".join(ch if ch.isalnum() else "-" for ch in text.lower()).strip("-")
    return slug or fallback


# Accepted --group values mapped to the path field they group by.
GROUP_FIELDS = {
    "sample": "sample",
    "samp": "sample",
    "ref": "ref",
    "reference": "ref",
    "reg": "reg",
    "region": "reg",
    "graph": "graph",
    "branch": "branches",
    "branches": "branches",
    "top": "top",
    "out": "top",
}


def _group_files(input_graphs: Iterable[Path], group: str):
    """Group input files by the requested field, preserving graph ordering."""
    group = (group or "sample").lower()
    if group == "all":
        field = None
    else:
        field = GROUP_FIELDS.get(group)
        if field is None:
            logger.warning(
                f"Unknown --group value {group!r}; grouping by sample instead. "
                f"Valid values: all, {', '.join(sorted(set(GROUP_FIELDS)))}."
            )
            field = "sample"
    groups: dict[str, list] = {}
    for file_path in input_graphs:
        try:
            path_segs = path.parse(file_path, path.GRAPH_SEGS)
        except Exception:
            # SVG structure drawings follow the draw path convention.
            try:
                path_segs = path.parse(file_path, path.DRAW_SEGS)
            except Exception:
                path_segs = {}
        is_structure = (
            path_segs.get("step") == "fold" or file_path.suffix.lower() == ".svg"
        )
        if field is None:
            group_key = "All graphs"
        elif field == "branches":
            branches = path_segs.get("branches") or []
            group_key = f"[{', '.join(branches)}]" if branches else "(no branches)"
        elif field == "graph" and is_structure:
            # Structure drawings have heterogeneous names, so collect them into
            # one group rather than scattering them by individual name.
            group_key = "Structure drawings"
        else:
            group_key = path_segs.get(field) or "unknown"
        groups.setdefault(group_key, []).append((path_segs, file_path))
    return [
        (group_key, sorted(files, key=lambda x: x[0].get("graph", "")))
        for group_key, files in sorted(groups.items())
    ]


def _build_file_entry(
    metadata: dict,
    file_path: Path,
    file_id: str,
    out_path: Path,
    portable: bool,
    plotly_libs: dict[str, str],
) -> dict:
    """Build the template context for a single graph/drawing."""
    name = file_path.stem
    branches = metadata.get("branches")
    if isinstance(branches, list):
        branches = f"[{', '.join(branches)}]" if branches else "[]"
    else:
        branches = branches or "unknown"
    entry = {
        "name": name,
        "id": file_id,
        "sample": metadata.get("sample", "unknown"),
        "reference": metadata.get("ref", "unknown"),
        "region": metadata.get("reg", "unknown"),
        "branches": branches,
    }
    is_svg = file_path.suffix.lower() == ".svg"

    if is_svg:
        try:
            entry["svg"] = file_path.read_text(encoding="utf-8")
            entry["kind"] = "svg"
        except Exception:
            entry["kind"] = "error"
    elif portable:
        try:
            html = file_path.read_text(encoding="utf-8")
        except Exception:
            html = None
        parsed = parse_plotly_html(html) if html else None
        if parsed is not None:
            # Share one copy of plotly.js per unique version.
            plotly_libs.setdefault(parsed["version"], parsed["lib_js"])
            entry["kind"] = "plotly"
            entry["plot_div"] = parsed["div_html"]
            entry["div_id"] = parsed["div_id"]
            entry["data_js"] = parsed["data_js"]
            entry["plotly_version"] = parsed["version"]
            entry["_config_js"] = parsed["config_js"]
        elif html is not None:
            # Unknown HTML: embed verbatim in a sandboxed iframe.
            entry["kind"] = "srcdoc"
            entry["srcdoc"] = html.replace("&", "&amp;").replace('"', "&quot;")
        else:
            entry["kind"] = "error"
    else:
        # Non-portable: link to the graph file via an iframe.
        entry["kind"] = "iframe"
        entry["src"] = os.path.relpath(file_path, out_path.parent)

    # Position and size are assigned per-group by _layout_group().
    return entry


def collate_graphs(
    input_graphs: Iterable[Path],
    out_path: Path,
    group: str,
    portable: bool,
    title: str = "SEISMIC-RNA Report",
    **kwargs,
) -> Path:
    """Collate HTML graph files (and SVG drawings) into one HTML file.

    Parameters
    ----------
    input_graphs: Iterable[Path]
        Paths to the HTML and/or SVG files to collate.
    out_path: Path
        Path to write the output collated HTML file.
    group: str
        Field by which to group graphs in the output (e.g. "sample",
        "ref", "reg", "branches", or "all").
    portable: bool
        If True, embed all resources so the HTML file is self-contained.
        Portable Plotly graphs of the same version share one copy of
        plotly.js; graphs made with a different plotly.js version get
        their own copy.
    title: str
        Title shown in the report header and browser tab.
    """
    sorted_groups = _group_files(input_graphs, group)

    plotly_libs: dict[str, str] = {}
    groups_data = []
    file_ids: set[str] = set()
    plotly_config = ""
    for grp_idx, (group_key, files_list) in enumerate(sorted_groups):
        group_id = f"{_slug(group_key, 'group')}-{grp_idx}"
        files_data = []
        for idx, (metadata, file_path) in enumerate(files_list):
            file_id = f"{_slug(file_path.stem, 'file')}-{grp_idx}_{idx}"
            if file_id in file_ids:
                raise ValueError("Multiple files have the same file_id.")
            file_ids.add(file_id)
            entry = _build_file_entry(
                metadata, file_path, file_id, out_path, portable, plotly_libs
            )
            if entry.get("kind") == "plotly" and not plotly_config:
                plotly_config = entry.get("_config_js", "")
            entry.pop("_config_js", None)
            files_data.append(entry)
        _layout_group(files_data)
        groups_data.append(
            {"group_name": group_key, "group_id": group_id, "files": files_data}
        )

    contains_svg = any(f.get("kind") == "svg" for g in groups_data for f in g["files"])
    n_graphs = sum(len(g["files"]) for g in groups_data)
    subtitle = (
        f"{n_graphs} graph{'s' if n_graphs != 1 else ''}"
        f" across {len(groups_data)} group{'s' if len(groups_data) != 1 else ''}"
    )

    # A stable per-report id so layout preferences persist across reopenings
    # but do not collide between different reports.
    report_hash = sha1("|".join(sorted(file_ids)).encode("utf-8")).hexdigest()[:10]
    report_id = f"{_slug(out_path.stem, 'report')}-{report_hash}"
    report_config = json.dumps(
        {"reportId": report_id, "groupIds": [g["group_id"] for g in groups_data]}
    )

    template = Template(_load_asset("report.html.jinja"))
    out_path.write_text(
        template.render(
            title=title,
            subtitle=subtitle,
            groups=groups_data,
            plotly_libs=plotly_libs,
            plotly_config=plotly_config,
            contains_svg=contains_svg,
            logo_svg=draw_seismic_logo(report=True),
            report_config=report_config,
            gridstack_css=_load_asset("gridstack.css"),
            report_css=_load_asset("report.css"),
            gridstack_js=_load_asset("gridstack-all.js"),
            sortable_js=_load_asset("sortable.js") if len(groups_data) > 1 else "",
            svg_pan_zoom_js=_load_asset("svg-pan-zoom.js") if contains_svg else "",
            report_js=_load_asset("report.js"),
        ),
        encoding="utf-8",
    )
    return out_path
