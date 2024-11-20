import re
import os
from functools import cached_property
from typing import Iterable
from pathlib import Path

from ..core import path
from ..core.extern.shell import args_to_cmd, run_cmd, JAVA_CMD, JAR_CMD
from ..core.task import dispatch
from ..fold.report import FoldReport
from ..core.write import need_write

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
    data{
    {% for key, val in color_dict.items() %}
        {{ key }} to {{ val }}
    {% endfor %}
    }

    theme {
        details {
            value = {{ details_value }}
        }

        {% for block in color_blocks %}
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
        {% endfor %}
    }
}
"""

TEMPLATE = Template(TEMPLATE_STRING)
RNARTIST_PATH = os.environ["RNARTISTCORE"]

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

def parse_db_file(file_path, struct_nums):
    # List to hold all dictionaries
    entries = dict()
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # The sequence line is after the first energy line
        seq_line = lines[1].strip()
        # Start reading structures from the second line (0-based index 2)
        struct_num = 0
        for i in range(2, len(lines)):
            line = lines[i].strip()
            # Check if line starts with '>', if so, skip
            if line.startswith(">"):
                continue
            # The lines that remain are structure lines
            structure_line = line
            # Create dictionary for this block
            entry = {
                'seq': seq_line,
                'value': structure_line
            }
            if struct_num in struct_nums or -1 in struct_nums:
                entries[struct_num] = entry
            struct_num += 1
    return entries

def parse_color_file(file_path):
    color_dict = {}
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
    color_blocks = list()

    block1 = ColorBlock("N", "#caccce")
    block2 = ColorBlock("n", "black")
    block3 = ColorBlock("N", "#88CCEE", "#7287D9", {"between": (0.0, 0.2)})
    block4 = ColorBlock("N", "#7287D9", "#8750D1", {"between": (0.2, 0.4)})
    block5 = ColorBlock("N", "#8750D1", "#852075", {"between": (0.4, 0.8)})
    block6 = ColorBlock("N", "#852075", "#661100", {"between": (0.8, 1.0)})
    block7 = ColorBlock("G", "#caccce")
    block8 = ColorBlock("U", "#caccce")
    block9 = ColorBlock("n", "white", "white", {"between": (0.2, 1.0)})
    block10 = ColorBlock("n", "black", "black", {"between": (0.0, 0.2)})
    block11 = ColorBlock("g", "black")
    block12 = ColorBlock("u", "black")

    color_blocks = [block1, block2, block3, block4, block5, block6,
                    block7, block8, block9, block10, block11, block12]

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

    def __init__(self,
                 report: Path,
                 tmp_dir: Path,
                 struct_nums: Iterable[int],
                 n_procs: int):
        self.top, self.fields = FoldReport.parse_path(report)
        self.report = FoldReport.load(report)
        self.tmp_dir = tmp_dir
        self.struct_nums = struct_nums
        self.n_procs = n_procs

    @property
    def sample(self):
        return self.fields.get(path.SAMP)

    @property
    def ref(self):
        return self.fields.get(path.REF)

    @property
    def sect(self):
        return self.fields.get(path.SECT)

    @property
    def profile(self):
        return self.fields.get(path.PROFILE)

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
                path.SECT: self.sect}

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
        return path.builddir(*path.SECT_DIR_SEGS, **self._get_dir_fields(top))

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
        return parse_color_file(self.get_varna_color_file(self.top))

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
        db_file = self.get_db_file(self.top)
        structs = parse_db_file(db_file, struct_nums=self.struct_nums)
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
         struct_nums: Iterable[int],
         tmp_dir: Path,
         keep_tmp: bool,
         n_procs: int,
         force: bool = False):
    """ Draw RNA structure(s) from a FoldReport. """
    rnartist = RNArtistRun(report_path,
                           tmp_dir,
                           struct_nums,
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
