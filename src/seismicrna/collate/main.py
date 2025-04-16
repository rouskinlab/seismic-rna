import os

from click import command
from pathlib import Path
from typing import Iterable

from .collate import collate_graphs
from ..core import path
from ..core.arg import (CMD_COLLATE,
                        arg_input_path,
                        opt_name,
                        opt_verbose_name,
                        opt_collate_out_dir,
                        opt_include_svg,
                        opt_include_graph,
                        opt_group,
                        opt_portable,
                        opt_force)
from ..core.logs import logger
from ..core.run import run_func
from ..core.write import need_write


def get_out_path(name: str,
                 verbose_name: bool,
                 input_files: list[Path],
                 collate_out_dir: str | Path | None):
    if collate_out_dir is not None:
        collate_out_dir = Path(collate_out_dir)
    if verbose_name:
        types = set()
        samples = set()
        refs = set()
        regs = set()
        for input_file in input_files:
            try:
                path_segs = path.parse(input_file, path.GRAPH_SEGS)
            except path.PathValueError:
                try:
                    path_segs = path.parse(input_file, path.DRAW_SEGS)
                except path.PathValueError:
                    continue
            assert path_segs.get("ext", None) in ('.html', '.svg')
            if path_segs.get("step") == "fold":
                graph = "structure-model"
            else:
                graph = path_segs.get("graph", None)
            if "collated" in graph:
                continue
            types.add(graph.split("_")[0])
            samples.add(path_segs.get("sample", None))
            refs.add(path_segs.get("ref", None))
            regs.add(path_segs.get("reg", None))

        types = sorted(list(types))
        samples = sorted(list(samples))
        refs = sorted(list(refs))
        regs = sorted(list(regs))

        info = f"{'-'.join(samples)}_{'-'.join(refs)}_{'-'.join(regs)}_{'-'.join(types)}"

    if not collate_out_dir:
        try:
            collate_out_dir = Path(os.path.commonpath(input_files))
        except ValueError:
            return None
        if collate_out_dir.is_file():
            collate_out_dir = collate_out_dir.parent
    if verbose_name:
        out_segs = [path.CollateInfoSeg]
        out_field_values = {path.COLLATE_NAME: name,
                            path.COLLATE_INFO: info,
                            path.EXT: path.HTML_EXT,
                            path.TOP: collate_out_dir}
    else:
        out_segs = [path.CollateSeg]
        out_field_values = {path.COLLATE_NAME: name,
                            path.EXT: path.HTML_EXT,
                            path.TOP: collate_out_dir}

    return path.build(out_segs, out_field_values)

@run_func(CMD_COLLATE)
def run(input_path: Iterable[str | Path],
        name: str,
        verbose_name: bool,
        include_svg: bool,
        include_graph: bool, *,
        group: str = 'sample',
        portable: bool = False,
        collate_out_dir: str | Path | None = None,
        force: bool,
        **kwargs) -> list[Path]:
    """ Collate HTML graphs into one HTML file. """
    # Collect all HTML graphs.
    input_files = list()
    if include_graph:
        for graph in path.find_files_chain(input_path, [path.HtmlSeg]):
            try:
                graph_segs = path.parse(graph, path.GRAPH_SEGS)
                if graph_segs.get("step", None) == "graph" and "collated" not in graph.name:
                    input_files.append(graph)
            except path.PathValueError:
                continue
    if include_svg:
        for svg in path.find_files_chain(input_path, [path.SvgSeg]):
            try:
                svg_segs = path.parse(svg, path.DRAW_SEGS)
                if svg_segs.get("step", None) == "fold":
                    input_files.append(svg)
            except path.PathValueError:
                    continue

    if len(input_files) == 0:
        logger.warning("No files found to collate.")
        return None

    out_path = get_out_path(name, verbose_name, input_files, collate_out_dir=collate_out_dir)

    if need_write(out_path, force=force):
        collate_graphs(input_files,
                       out_path,
                       group,
                       portable,
                       **kwargs)

    return out_path


params = [
    arg_input_path,
    opt_name,
    opt_verbose_name,
    opt_collate_out_dir,
    opt_include_svg,
    opt_include_graph,
    opt_group,
    opt_portable,
    opt_force
]


@command(CMD_COLLATE, params=params)
def cli(*args, **kwargs):
    """ Collate HTML graphs and SVG drawings into an HTML report file. """
    return run(*args, **kwargs)
