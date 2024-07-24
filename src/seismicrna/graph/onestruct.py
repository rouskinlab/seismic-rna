from abc import ABC
from functools import cached_property
from logging import getLogger
from pathlib import Path

from .base import cgroup_table
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import OneRelGraph
from ..core import path
from ..core.arg import (opt_struct_file,
                        opt_fold_sections_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_full)
from ..core.rna import RNAState, from_ct
from ..core.seq import DNA, RefSections
from ..fold.main import find_foldable_tables

logger = getLogger(__name__)


class StructOneTableGraph(OneTableGraph, OneRelGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *,
                 struct_file: Path | None,
                 struct_sect: str | None,
                 **kwargs):
        super().__init__(**kwargs)
        self._struct_file = struct_file
        self._struct_sect = struct_sect

    @property
    def struct_sect(self):
        """ Name of the section from which the structure comes. """
        if self._struct_file is not None:
            # Use the section from the given structure file.
            fields = path.parse(self._struct_file,
                                path.RefSeg,
                                path.SectSeg,
                                path.ConnectTableSeg)
            return fields[path.SECT]
        if self._struct_sect is None:
            raise ValueError("A structure section is required if no structure "
                             "file is given")
        return self._struct_sect

    def get_path_fields(self):
        fields = super().get_path_fields()
        # Replace the section with the structure section.
        fields[path.SECT] = self.struct_sect
        return fields

    def _get_struct_files(self, profile: str):
        """ Get the path of the CT file of RNA structures for a specific
        RNA profile. """
        if self._struct_file is not None:
            # Use the given structure file for every profile.
            return self._struct_file
        # Determine the path of the structure file from the profile.
        return path.build(*path.CT_FILE_SEGS,
                          top=self.top,
                          sample=self.sample,
                          cmd=path.CMD_FOLD_DIR,
                          ref=self.ref,
                          sect=self.struct_sect,
                          profile=profile,
                          ext=path.CT_EXT)

    @cached_property
    def path_subject(self):
        return f"{self.table.sect}__{super().path_subject}"

    @cached_property
    def _title_main(self):
        return [f"{self.what()} of structures "
                f"for reference {repr(self.ref)} "
                f"section {repr(self._struct_sect)} "
                f"vs. {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over section {repr(self.sect)}"]

    def iter_profiles(self):
        """ Yield each RNAProfile from the table. """
        yield from self.table.iter_profiles(quantile=self.quantile,
                                            rel=self.rel_name,
                                            k=self.k,
                                            clust=self.clust)

    def iter_states(self):
        """ Yield each RNAState. """
        for profile in self.iter_profiles():
            ct_file = self._get_struct_files(profile.profile)
            try:
                for struct in from_ct(ct_file):
                    yield RNAState.from_struct_profile(struct, profile)
            except FileNotFoundError:
                logger.error(f"Structure file {ct_file} does not exist; please "
                             f"obtain the file, e.g. using seismic fold")


class StructOneTableWriter(OneTableWriter, ABC):

    def iter_graphs(self,
                    rels: tuple[str, ...],
                    cgroup: str,
                    struct_file: tuple[str, ...] = (),
                    fold_coords: tuple[tuple[str, int, int], ...] = (),
                    fold_primers: tuple[tuple[str, DNA, DNA], ...] = (),
                    fold_sections_file: str | None = None,
                    fold_full: bool = opt_fold_full.default,
                    **kwargs):
        struct_files = list()
        for file in struct_file:
            # Use a given CT file of an RNA structure, and determine the
            # structure section name from the file path.
            fields = path.parse(file,
                                path.RefSeg,
                                path.SectSeg,
                                path.ConnectTableSeg)
            ref = fields[path.REF]
            if ref == self.table.ref:
                struct_files.append(file)
            else:
                logger.warning(f"Skipped CT file {file} in section directory "
                               f"{path.SECT} in reference directory {ref}, "
                               f"which differs from the reference name of "
                               f"the table file ({self.table.ref})")
        # Add the sections from the given coordinates/primers.
        ref_sections = RefSections([(self.table.ref, self.table.refseq)],
                                   sects_file=(Path(fold_sections_file)
                                               if fold_sections_file
                                               else None),
                                   coords=fold_coords,
                                   primers=fold_primers,
                                   default_full=fold_full)
        fold_sects = [section.name
                      for section in ref_sections.list(self.table.ref)]
        if not fold_sects:
            # Add the table's section if no other sections were defined.
            fold_sects.append(self.table.sect)
        # Generate a graph for each cluster, relationship, and section.
        for cparams in cgroup_table(self.table, cgroup):
            kwparams = kwargs | cparams
            for rels_group in rels:
                for file in struct_files:
                    yield self.get_graph(rels_group,
                                         struct_file=file,
                                         struct_sect=None,
                                         **kwparams)
                for sect in fold_sects:
                    yield self.get_graph(rels_group,
                                         struct_file=None,
                                         struct_sect=sect,
                                         **kwparams)


class StructOneTableRunner(OneTableRunner, ABC):

    @classmethod
    def get_table_finder(cls):
        def table_finder(*args, **kwargs):
            for table_file, _ in find_foldable_tables(*args, **kwargs):
                yield table_file

        return table_finder

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_struct_file,
                                       opt_fold_sections_file,
                                       opt_fold_coords,
                                       opt_fold_primers,
                                       opt_fold_full]

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
