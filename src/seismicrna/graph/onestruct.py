from abc import ABC
from functools import cached_property
from pathlib import Path

from .base import cgroup_table
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import OneRelGraph
from ..core import path
from ..core.arg import (opt_struct_file,
                        opt_fold_regions_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_full)
from ..core.logs import logger
from ..core.rna import RNAState, from_ct
from ..core.seq import DNA, RefRegions
from ..fold.main import find_foldable_tables


class StructOneTableGraph(OneTableGraph, OneRelGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *,
                 struct_file: Path | None,
                 struct_reg: str | None,
                 **kwargs):
        super().__init__(**kwargs)
        self._struct_file = struct_file
        self._struct_reg = struct_reg

    @property
    def struct_reg(self):
        """ Name of the region from which the structure comes. """
        if self._struct_file is not None:
            # Use the region from the given structure file.
            fields = path.parse(self._struct_file,
                                path.RefSeg,
                                path.RegSeg,
                                path.ConnectTableSeg)
            return fields[path.REG]
        if self._struct_reg is None:
            raise ValueError("A structure region is required if no structure "
                             "file is given")
        return self._struct_reg

    def get_path_fields(self):
        fields = super().get_path_fields()
        # Replace the region with the structure region.
        fields[path.REG] = self.struct_reg
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
                          reg=self.struct_reg,
                          profile=profile,
                          ext=path.CT_EXT)

    @cached_property
    def path_subject(self):
        return f"{self.table.reg}__{super().path_subject}"

    @cached_property
    def _title_main(self):
        return [f"{self.what()} of structures "
                f"for reference {repr(self.ref)} "
                f"region {repr(self._struct_reg)} "
                f"vs. {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over region {repr(self.reg)}"]

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
                logger.error(f"Structure file {ct_file} does not exist; "
                             "please obtain the file, e.g. using seismic fold")


class StructOneTableWriter(OneTableWriter, ABC):

    def iter_graphs(self,
                    rels: tuple[str, ...],
                    cgroup: str,
                    struct_file: tuple[str, ...] = (),
                    fold_coords: tuple[tuple[str, int, int], ...] = (),
                    fold_primers: tuple[tuple[str, DNA, DNA], ...] = (),
                    fold_regions_file: str | None = None,
                    fold_full: bool = opt_fold_full.default,
                    **kwargs):
        struct_files = list()
        for file in struct_file:
            # Use a given CT file of an RNA structure, and determine the
            # structure region name from the file path.
            fields = path.parse(file,
                                path.RefSeg,
                                path.RegSeg,
                                path.ConnectTableSeg)
            ref = fields[path.REF]
            if ref == self.table.ref:
                struct_files.append(file)
            else:
                logger.warning(f"Skipped CT file {file} in region directory "
                               f"{repr(path.REG)} in reference directory "
                               f"{repr(ref)}, which differs from the reference "
                               f"name of the table file {repr(self.table.ref)}")
        # Add the regions from the given coordinates/primers.
        ref_regions = RefRegions([(self.table.ref, self.table.refseq)],
                                 regs_file=(Path(fold_regions_file)
                                             if fold_regions_file
                                             else None),
                                 coords=fold_coords,
                                 primers=fold_primers,
                                 default_full=fold_full)
        fold_regs = [region.name
                      for region in ref_regions.list(self.table.ref)]
        if not fold_regs:
            # Add the table's region if no other regions were defined.
            fold_regs.append(self.table.reg)
        # Generate a graph for each cluster, relationship, and region.
        for cparams in cgroup_table(self.table, cgroup):
            kwparams = kwargs | cparams
            for rels_group in rels:
                for file in struct_files:
                    yield self.get_graph(rels_group,
                                         struct_file=file,
                                         struct_reg=None,
                                         **kwparams)
                for reg in fold_regs:
                    yield self.get_graph(rels_group,
                                         struct_file=None,
                                         struct_reg=reg,
                                         **kwparams)


class StructOneTableRunner(OneTableRunner, ABC):

    @classmethod
    def get_table_loader(cls):
        return find_foldable_tables

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_struct_file,
                                       opt_fold_regions_file,
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
