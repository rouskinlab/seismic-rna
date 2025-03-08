from typing import Iterable
from click import command
from pathlib import Path
from ..core import path
from ..core.arg import (CMD_DRAW,
                        arg_input_path,
                        opt_force,
                        opt_keep_tmp,
                        opt_struct_num,
                        opt_color,
                        opt_verify_times,
                        opt_draw_svg,
                        opt_draw_png,
                        opt_update_rnartistcore,
                        opt_num_cpus)
from ..core.run import run_func
from ..core.task import dispatch, as_list_of_tuples
from .draw import draw


@run_func(CMD_DRAW, with_tmp=True, pass_keep_tmp=True)
def run(input_path: Iterable[str | Path], *,
        struct_num: Iterable[int],
        color: bool,
        draw_svg: bool,
        draw_png: bool,
        update_rnartistcore: bool,
        force: bool,
        num_cpus: int,
        tmp_dir: Path,
        verify_times: bool,
        keep_tmp: bool) -> list[Path]:
    """ Draw RNA structures with reactivities using RNArtistCore. """
    # Generate the positional arguments for draw.
    args = as_list_of_tuples(path.find_files_chain(input_path,
                                                   [path.FoldRepSeg]))
    # Draw the files.
    return dispatch(draw,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=args,
                    kwargs=dict(struct_num=struct_num,
                                color=color,
                                draw_svg=draw_svg,
                                draw_png=draw_png,
                                update=update_rnartistcore,
                                tmp_dir=tmp_dir,
                                keep_tmp=keep_tmp,
                                verify_times=verify_times,
                                force=force))


params = [
    arg_input_path,
    opt_struct_num,
    opt_color,
    opt_draw_svg,
    opt_draw_png,
    opt_update_rnartistcore,
    opt_force,
    opt_keep_tmp,
    opt_verify_times,
    opt_num_cpus
]


@command(CMD_DRAW, params=params)
def cli(*args, **kwargs):
    """ Draw RNA structures with reactivities using RNArtistCore. """
    return run(*args, **kwargs)
