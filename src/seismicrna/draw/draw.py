import re
import os
from functools import cached_property
from typing import Iterable
from pathlib import Path

from ..core import path
from ..core.logs import logger
from ..core.task import dispatch
from ..core.write import need_write
from ..core.extern.shell import args_to_cmd, run_cmd, JAVA_CMD, JAR_CMD
from ..core.rna.io import from_ct
from ..core.rna.state import RNAState
from ..fold.report import FoldReport
from ..core.header import AVERAGE_PREFIX

from ..relate.table import RelatePositionTable, RelatePositionTableLoader
from ..mask.table import MaskPositionTable, MaskPositionTableLoader
from ..cluster.table import ClusterPositionTable, ClusterPositionTableLoader


from jinja2 import Template

TEMPLATE_STRING = """
rnartist {
    svg {
        path = "{{ path }}"
    }
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
RNARTIST_PATH = os.environ.get("RNARTISTCORE")

TABLES = {AVERAGE_PREFIX: (MaskPositionTable,
                           MaskPositionTableLoader),
          path.CMD_CLUST_DIR: (ClusterPositionTable,
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
                 path: Path,
                 seq: str,
                 value: str,
                 name: str,
                 color_dict: dict,
                 details_value: int,
                 color_blocks: list[ColorBlock]):
        self.path = path
        self.seq = seq
        self.value = value
        self.name = name
        self.color_dict = color_dict
        self.details_value = details_value
        self.color_blocks = color_blocks

    def to_dict(self):
        return {
            'path': self.path,
            'seq': self.seq,
            'value': self.value,
            'name': self.name,
            'color_dict': self.color_dict,
            'details_value': self.details_value,
            'color_blocks': [block.to_dict() if hasattr(block, 'to_dict') else block for block in self.color_blocks]
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

    jinja_data = JinjaData(path=out_dir,
                           seq=struct["seq"],
                           value=struct["value"],
                           name=name,
                           color_dict=color_dict,
                           details_value=5,
                           color_blocks=color_blocks)
    return jinja_data


class RNArtistRun(object):

    def _parse_profile(self):
        match = re.search(r'(.+?)__(?:(average|.+?(?=-)))', self.profile)
        self.data_reg = match.group(1) if match else None
        self.table_type = match.group(2) if match else None
        if not match:
            logger.warning(f"Could not parse profile: {self.profile}.")

    def __init__(self,
                 report: Path,
                 tmp_dir: Path,
                 struct_num: Iterable[int],
                 color: bool,
                 n_procs: int):
        self.top, self.fields = FoldReport.parse_path(report)
        self.sample = self.fields.get(path.SAMP)
        self.ref = self.fields.get(path.REF)
        self.reg = self.fields.get(path.REG)
        self.profile = self.fields.get(path.PROFILE)
        self.report = FoldReport.load(report)
        self.tmp_dir = tmp_dir
        self.struct_num = list(struct_num)
        self.color = color
        self.n_procs = n_procs
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
                path.CMD: path.CMD_FOLD_DIR,
                path.SAMP: self.sample,
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
        return path.builddir(*path.REG_DIR_SEGS, **self._get_dir_fields(top))

    def _get_file(self, top: Path, file_seg: path.Segment, **file_fields):
        """ Get the path to a file of the RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.
        file_seg: path.Segment
            Segment of the file component of the path.
        **file_fields
            Fields for the file segment.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_dir(top).joinpath(file_seg.build(**file_fields))

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
        return (self.table_class.build_path(top=self.top,
                                            sample=self.sample,
                                            ref=self.ref,
                                            reg=self.data_reg)
                if self.table_class else None)

    @cached_property
    def table(self):
        if self.table_file.exists():
            return (self.table_loader(self.table_file)
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
                              profile=self.profile,
                              ext=path.CT_EXT)

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
                              profile=self.profile,
                              ext=path.DOT_EXTS[0])

    def get_svg_file(self, top: Path, struct: int):
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
                              path.SvgSeg,
                              profile=self.profile,
                              struct=struct,
                              ext=path.SVG_EXT)

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
                              profile=self.profile,
                              ext=path.TXT_EXT)

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
                              profile=self.profile,
                              struct=struct,
                              ext=path.KTS_EXT)

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
                       out_path: Path,
                       script_file: Path,
                       keep_tmp: bool,
                       force: bool):
        if need_write(out_path, force=force):
            jinja_data = build_jinja_data(struct=struct,
                                          color_dict=self.color_dict,
                                          name=struct_name,
                                          out_dir=out_path.parent,
                                          highlight_pos=self.edited_numbers)

            rnartist_script = TEMPLATE.render(jinja_data.to_dict())
            with open(script_file, 'w') as f:
                f.write(rnartist_script)
            rnartist_cmd = args_to_cmd([JAVA_CMD,
                                        JAR_CMD,
                                        RNARTIST_PATH,
                                        script_file])
            run_cmd(rnartist_cmd)
            if not keep_tmp:
                script_file.unlink(missing_ok=True)
        return out_path

    def run(self, keep_tmp: bool, force: bool):
        structs = dict()
        if not self.struct_num:
            self.struct_num = (self.best_struct,)
        for struct_num, struct in enumerate(
                from_ct(self.get_ct_file(self.top))):
            if struct_num in self.struct_num:
                structs[struct_num] = dict(seq=struct.seq,
                                           value=struct.db_structure)
        args = [
            (f"{self.profile}-{struct_num}",
             struct,
             self.get_svg_file(self.top, struct=struct_num),
             self.get_script_file(top=self.tmp_dir, struct=struct_num))
            for struct_num, struct in structs.items()]
        results = dispatch(self.process_struct,
                           self.n_procs,
                           args=args,
                           pass_n_procs=False,
                           kwargs=dict(force=force,
                                       keep_tmp=keep_tmp))
        return results


def draw(report_path: Path, *,
         struct_num: Iterable[int],
         color: bool,
         tmp_dir: Path,
         keep_tmp: bool,
         n_procs: int,
         force: bool = False):
    """ Draw RNA structure(s) from a FoldReport. """
    rnartist = RNArtistRun(report_path,
                           tmp_dir,
                           struct_num,
                           color,
                           n_procs)
    # By convention, a function must return a Path for dispatch to deem
    # that it has completed successfully.
    return rnartist.run(keep_tmp, force)

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
