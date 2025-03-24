import os
import re
from functools import cached_property
from pathlib import Path
from typing import Iterable

from jinja2 import Template

from ..cluster.data import ClusterPositionTable, ClusterPositionTableLoader
from ..core import path
from ..core.extern.shell import (args_to_cmd,
                                 run_cmd,
                                 JAVA_CMD,
                                 JAR_CMD,
                                 JGO_CMD,
                                 ShellCommandFailedError)
from ..core.header import AVERAGE_PREFIX
from ..core.logs import logger
from ..core.report import SampleF, BranchesF, RefF, RegF, ProfileF
from ..core.rna.io import from_ct
from ..core.rna.state import RNAState
from ..core.task import dispatch
from ..core.write import need_write
from ..fold.report import FoldReport
from ..mask.table import MaskPositionTable, MaskPositionTableLoader

TEMPLATE_STRING = """
rnartist {
    {% if draw_svg %}
    svg {
        path = "{{ path }}"
    }
    {% endif %}
    {% if draw_png %}
    png {
        path = "{{ path }}"
    }
    {% endif %}
    ss {
        bn {
            seq = "{{ seq }}"
            value = "{{ value }}"
            name = "{{ name }}"
        }
    }
    {% if color_dict %}
    data {
        {% for key, val in color_dict.items() %}
        {{ key }} to {{ val }}
        {% endfor %}
    }
    {% endif %}
    theme {
        details {
            value = {{ details_value }}
        }
        {% for block in color_blocks %}
        {% if color_dict or not block.color_to %}
        color {
            type = "{{ block.color_type }}"
            value = "{{ block.color_value }}"
            {% if block.color_to %}
            to = "{{ block.color_to }}"
            {% endif %}
            {% if block.color_filter %}
            {% for filter_type, filter_value in block.color_filter.items() %}
            {% if filter_type == 'between' %}
            data between {{ filter_value[0] }}..{{ filter_value[1] }}
            {% else %}
            data {{filter_type}} {{ filter_value }}
            {% endif %}
            {% endfor %}
            {% endif %}
            {% if block.location %}
            location {
                {{ block.location[0] }} to {{ block.location[1] }}
            }
            {% endif %}
        }
        {% endif %}
        {% endfor %}
    }
}
"""

TEMPLATE = Template(TEMPLATE_STRING)

RNARTIST_REPO = "maven-snapshots=https://oss.sonatype.org/content/repositories/snapshots"
RNARTIST_GROUP_ID = "io.github.fjossinet.rnartist"
RNARTIST_ARTIFACT_ID = "rnartistcore"
RNARTIST_FALLBACK_VERSION = "0.4.7-SNAPSHOT"
RNARTIST_VERSION = os.getenv("RNARTISTCORE_VERSION", RNARTIST_FALLBACK_VERSION)
RNARTIST_MAIN_CLASS = "io.github.fjossinet.rnartist.core.MainKt"
RNARTIST_ARTIFACT_NO_VERSION = f"{RNARTIST_GROUP_ID}:{RNARTIST_ARTIFACT_ID}:{RNARTIST_MAIN_CLASS}"
RNARTIST_ARTIFACT_WITH_VERSION = f"{RNARTIST_GROUP_ID}:{RNARTIST_ARTIFACT_ID}:{RNARTIST_VERSION}:{RNARTIST_MAIN_CLASS}"
RNARTIST_ARTIFACT_FALLBACK = f"{RNARTIST_GROUP_ID}:{RNARTIST_ARTIFACT_ID}:{RNARTIST_FALLBACK_VERSION}:{RNARTIST_MAIN_CLASS}"

TABLES = {AVERAGE_PREFIX: (MaskPositionTable,
                           MaskPositionTableLoader),
          path.CLUSTER_STEP: (ClusterPositionTable,
                              ClusterPositionTableLoader)}


class ColorBlock:
    def __init__(self,
                 color_type: str,
                 color_value: str,
                 color_to: str = None,
                 color_filter: str = None,
                 location: tuple[int] = None):
        self.color_type = color_type
        self.color_value = color_value
        self.color_to = color_to
        self.color_filter = color_filter
        self.location = location

    def to_dict(self):
        return {
            'color_type': self.color_type,
            'color_value': self.color_value,
            'color_to': self.color_to,
            'color_filter': self.color_filter,
            'location': self.location
        }


class JinjaData:
    def __init__(self,
                 path_: Path,
                 seq: str,
                 value: str,
                 name: str,
                 color_dict: dict,
                 details_value: int,
                 color_blocks: list[ColorBlock],
                 draw_svg: bool,
                 draw_png: bool):
        self.path = path_
        self.seq = seq
        self.value = value
        self.name = name
        self.color_dict = color_dict
        self.details_value = details_value
        self.color_blocks = color_blocks
        self.draw_svg = draw_svg
        self.draw_png = draw_png

    def to_dict(self):
        return {
            'path': self.path,
            'draw_svg': self.draw_svg,
            'draw_png': self.draw_png,
            'seq': self.seq,
            'value': self.value,
            'name': self.name,
            'color_dict': self.color_dict,
            'details_value': self.details_value,
            'color_blocks': [block.to_dict()
                             if hasattr(block, 'to_dict')
                             else block
                             for block in self.color_blocks]
        }


def parse_color_file(file_path):
    color_dict = dict()
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line_num, line in enumerate(lines, start=1):
            pos = line_num
            value = float(line)
            color_dict[pos] = value
    return color_dict


def build_jinja_data(struct: str,
                     color_dict: dict,
                     name: str,
                     out_dir: Path,
                     draw_svg: bool,
                     draw_png: bool,
                     highlight_pos: Iterable[int] = None):
    color_blocks = [ColorBlock("N", "#caccce"),
                    ColorBlock("n", "black"),
                    ColorBlock("N", "#88CCEE", "#7287D9",
                               {"between": (0.0, 0.2)}),
                    ColorBlock("N", "#7287D9", "#8750D1",
                               {"between": (0.2, 0.4)}),
                    ColorBlock("N", "#8750D1", "#852075",
                               {"between": (0.4, 0.8)}),
                    ColorBlock("N", "#852075", "#661100",
                               {"between": (0.8, 1.0)}),
                    ColorBlock("n", "white", "white", {"between": (0.2, 1.0)}),
                    ColorBlock("n", "black", "black", {"between": (0.0, 0.2)})]

    if highlight_pos:
        for pos in highlight_pos:
            color_blocks.append(ColorBlock("N", "#00c90a",
                                           location=(pos, pos)))
            color_blocks.append(ColorBlock("n", "black",
                                           location=(pos, pos)))

    jinja_data = JinjaData(path_=out_dir,
                           seq=struct["seq"],
                           value=struct["value"],
                           name=name,
                           color_dict=color_dict,
                           details_value=5,
                           color_blocks=color_blocks,
                           draw_svg=draw_svg,
                           draw_png=draw_png)
    return jinja_data


class RNArtistRun(object):

    def _parse_profile(self):
        match = re.search(r'(.+?)__(?:(average|.+?(?=-)))', self.profile)
        self.data_reg = match.group(1) if match else None
        self.table_type = match.group(2) if match else None
        if not match:
            logger.warning(f"Could not parse profile: {self.profile}.")

    def __init__(self,
                 report_file: Path,
                 tmp_dir: Path,
                 struct_num: Iterable[int],
                 color: bool,
                 verify_times: bool,
                 draw_svg: bool,
                 draw_png: bool,
                 update: bool,
                 num_cpus: int):
        self.top, _ = FoldReport.parse_path(report_file)
        report = FoldReport.load(report_file)
        self.sample = report.get_field(SampleF)
        self.branches = report.get_field(BranchesF)
        self.ref = report.get_field(RefF)
        self.reg = report.get_field(RegF)
        self.profile = report.get_field(ProfileF)
        self.tmp_dir = tmp_dir
        self.struct_num = list(struct_num)
        self.color = color
        self.draw_svg = draw_svg
        self.draw_png = draw_png
        self.update = update
        self.num_cpus = num_cpus
        self.verify_times = verify_times
        self._parse_profile()

    def _get_dir_fields(self, top: Path):
        """ Get the path fields for the directory of this RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        dict[str, str | pathlib.Path]
            Path fields.
        """
        return {path.TOP: top,
                path.STEP: path.FOLD_STEP,
                path.SAMPLE: self.sample,
                path.BRANCHES: self.branches,
                path.REF: self.ref,
                path.REG: self.reg}

    def _get_dir(self, top: Path):
        """ Get the directory in which to write files of this RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Parent directory of files for this RNA.
        """
        return path.builddir(path.REG_DIR_SEGS, self._get_dir_fields(top))

    def _get_file(self, top: Path, file_seg: path.PathSegment, file_fields):
        """ Get the path to a file of the RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.
        file_seg: path.PathSegment
            Segment of the file component of the path.
        file_fields
            Fields for the file segment.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_dir(top).joinpath(file_seg.build(file_fields))

    @cached_property
    def table_classes(self):
        if self.table_type and self.table_type in TABLES:
            return TABLES.get(self.table_type)
        else:
            logger.warning(f"Cannot use table of type {self.table_type} to "
                           "calculate AUROC.")
        return None, None

    @property
    def table_class(self):
        return self.table_classes[0]

    @property
    def table_loader(self):
        return self.table_classes[1]

    @property
    def table_file(self):
        table_branches = self.branches.copy()
        table_branches.pop("fold", None)
        return (self.table_class.build_path({path.TOP: self.top,
                                             path.SAMPLE: self.sample,
                                             path.BRANCHES: table_branches,
                                             path.REF: self.ref,
                                             path.REG: self.data_reg})
                if self.table_class else None)

    @cached_property
    def table(self):
        if self.table_file.exists():
            return (self.table_loader(self.table_file,
                                      verify_times=self.verify_times)
                    if self.table_loader else None)
        else:
            logger.warning(f"{self.table_file} does not exist.")
            return None

    def get_ct_file(self, top: Path):
        """ Get the path to the connectivity table (CT) file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.ConnectTableSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.CT_EXT})

    def get_db_file(self, top: Path):
        """ Get the path to the dot-bracket (DB) file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.DotBracketSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.DOT_EXTS[0]})

    def get_svg_file(self, top: Path, struct: int):
        """ Get the path to the SVG file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.SvgSeg,
                              {path.PROFILE: self.profile,
                               path.STRUCT: struct,
                               path.EXT: path.SVG_EXT})

    def get_png_file(self, top: Path, struct: int):
        """ Get the path to the PNG file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.PngSeg,
                              {path.PROFILE: self.profile,
                               path.STRUCT: struct,
                               path.EXT: path.PNG_EXT})

    def get_varna_color_file(self, top: Path):
        """ Get the path to the VARNA color file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.VarnaColorSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.TXT_EXT})

    def get_script_file(self, top: Path, struct: int):
        """ Get the path to the RNArtist script (.kts) file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.KtsSeg,
                              {path.PROFILE: self.profile,
                               path.STRUCT: struct,
                               path.EXT: path.KTS_EXT})

    @cached_property
    def best_struct(self):
        if not self.table:
            logger.warning(f"Could not find best structure for {self.profile}. "
                           f"Drawing the MFE structure by default.")
            return 0
        max_auc = 0
        best_auc = 0
        structs = [struct for struct in from_ct(self.get_ct_file(self.top))]
        regions = [structs[0].region]
        for profile in self.table.iter_profiles(regions=regions):
            if self.profile == profile.profile:
                for struct_num, struct in enumerate(structs):
                    state = RNAState.from_struct_profile(struct, profile)
                    if state.auc > max_auc:
                        max_auc = state.auc
                        best_auc = struct_num
        return best_auc

    @cached_property
    def edited_numbers(self):
        edited_pattern = re.compile(r'edited_(\d+(?:_\d+)*)')
        edited_numbers = list()
        match = edited_pattern.search(self.profile)
        if match:
            numbers = list(map(int, match.group(1).split('_')))
            edited_numbers.extend(numbers)
        return edited_numbers if edited_numbers else None

    @cached_property
    def color_dict(self):
        color_file = self.get_varna_color_file(self.top)
        if color_file.exists():
            if self.color:
                return parse_color_file(color_file)
        else:
            logger.warning(f"{color_file} does not exist. "
                           "Defaulting to --no-color")
            self.color = False
        return dict()

    def process_struct(self,
                       struct_name: str,
                       struct: str,
                       svg_path: Path,
                       png_path: Path,
                       script_file: Path,
                       keep_tmp: bool,
                       force: bool):
        if not (self.draw_svg or self.draw_png):
            logger.warning("Both --no-draw-svg and --no-draw-png are set, defaulting to --svg")
            self.draw_svg = True
        if (self.draw_svg and need_write(svg_path, force=force)) or (
                self.draw_png and need_write(png_path, force=force)):
            jinja_data = build_jinja_data(struct=struct,
                                          color_dict=self.color_dict,
                                          name=struct_name,
                                          out_dir=svg_path.parent,
                                          highlight_pos=self.edited_numbers,
                                          draw_svg=self.draw_svg,
                                          draw_png=self.draw_png)

            rnartist_script = TEMPLATE.render(jinja_data.to_dict())
            with open(script_file, 'w') as f:
                f.write(rnartist_script)
            if self.update:
                args_to_try = [[JGO_CMD,
                                "-U",
                                "-r",
                                RNARTIST_REPO,
                                RNARTIST_ARTIFACT_NO_VERSION,
                                script_file]]
            else:
                args_to_try = []
            args_to_try.extend([[JGO_CMD,
                                 "-r",
                                 RNARTIST_REPO,
                                 RNARTIST_ARTIFACT_NO_VERSION,
                                 script_file],
                                [JGO_CMD,
                                 "-r",
                                 RNARTIST_REPO,
                                 RNARTIST_ARTIFACT_WITH_VERSION,
                                 script_file],
                                [JGO_CMD,
                                 "-r",
                                 RNARTIST_REPO,
                                 RNARTIST_ARTIFACT_FALLBACK,
                                 script_file]])
            if RNARTIST_VERSION != RNARTIST_FALLBACK_VERSION:
                logger.action(f"Attempting to load RNArtistCore version {RNARTIST_VERSION}")
            for args in args_to_try:
                try:
                    rnartist_cmd = args_to_cmd(args)
                    run_cmd(rnartist_cmd)
                    break
                except ShellCommandFailedError:
                    continue
            else:
                logger.warning("Running RNArtistCore with jgo failed. Falling back to manual installation.")
                from ..core.arg import CMD_DRAW
                from ..core.extern import require_env_var
                require_env_var("RNARTISTCORE", CMD_DRAW)
                RNARTIST_PATH = os.environ.get("RNARTISTCORE")
                rnartist_cmd = args_to_cmd([JAVA_CMD,
                                            JAR_CMD,
                                            RNARTIST_PATH,
                                            script_file])
                run_cmd(rnartist_cmd)
            if not keep_tmp:
                script_file.unlink(missing_ok=True)
        out_paths = [path for path, write in zip((svg_path, png_path), (self.draw_svg, self.draw_png)) if write]
        return out_paths

    def run(self, keep_tmp: bool, force: bool):
        structs = dict()
        if not self.struct_num:
            self.struct_num = (self.best_struct,)
        for struct_num, struct in enumerate(
                from_ct(self.get_ct_file(self.top))):
            if struct_num in self.struct_num or -1 in self.struct_num:
                structs[struct_num] = dict(seq=struct.seq,
                                           value=struct.db_structure)
        args = [(f"{self.profile}-{struct_num}",
                 struct,
                 self.get_svg_file(self.top, struct=struct_num),
                 self.get_png_file(self.top, struct=struct_num),
                 self.get_script_file(top=self.tmp_dir, struct=struct_num))
                for struct_num, struct in structs.items()]
        return dispatch(self.process_struct,
                        num_cpus=self.num_cpus,
                        pass_num_cpus=False,
                        as_list=True,
                        ordered=False,
                        raise_on_error=False,
                        args=args,
                        kwargs=dict(force=force,
                                    keep_tmp=keep_tmp))


def draw(report_path: Path, *,
         struct_num: Iterable[int],
         color: bool,
         draw_svg: bool,
         draw_png: bool,
         update: bool,
         tmp_dir: Path,
         keep_tmp: bool,
         verify_times: bool,
         num_cpus: int,
         force: bool = False):
    """ Draw RNA structure(s) from a FoldReport. """
    rnartist = RNArtistRun(report_path,
                           tmp_dir,
                           struct_num,
                           color,
                           verify_times,
                           draw_svg,
                           draw_png,
                           update,
                           num_cpus)
    # By convention, a function must return a Path for dispatch to deem
    # that it has completed successfully.
    return rnartist.run(keep_tmp, force)
