from abc import ABC
from functools import cached_property
from logging import getLogger
from pathlib import Path

from .onetable import OneTableGraph, OneTableRunner
from .rel import OneRelGraph
from ..core import path
from ..core.arg import opt_struct_file, opt_struct_sect
from ..core.rna import RNAState, from_ct

logger = getLogger(__name__)


class StructOneTableGraph(OneTableGraph, OneRelGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *, struct_file: Path | None, struct_sect: str, **kwargs):
        super().__init__(**kwargs)
        if struct_file is not None:
            # Use a given CT file of an RNA structure, and determine the
            # structure section name from the file path.
            self._struct_file = struct_file
            fields = path.parse(self._struct_file,
                                path.RefSeg,
                                path.SectSeg,
                                path.ConnectTableSeg)
            ref = fields[path.REF]
            if ref != self.ref:
                raise ValueError(f"Reference names differ in CT file ({ref}) "
                                 f"and table file ({self.ref})")
            self.struct_sect = fields[path.SECT]
            if struct_sect and struct_sect != self.struct_sect:
                logger.warning(f"Structure section names differ in CT file "
                               f"({self.struct_sect}) and given struct_sect "
                               f"argument ({struct_sect}); using the former")
        else:
            # Use the given structure section name if one was given,
            # otherwise default to the section from which the data came.
            self.struct_sect = struct_sect if struct_sect else self.sect
            # The structure file will vary by RNA profile, so leave it
            # blank for now.
            self._struct_file = None

    def get_path_fields(self):
        fields = super().get_path_fields()
        # Replace the section with the structure section.
        fields[path.SECT] = self.struct_sect
        return fields

    def _get_struct_file(self, profile: str):
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
                f"section {repr(self.struct_sect)} "
                f"vs. {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over section {repr(self.sect)}"]

    def iter_profiles(self):
        """ Yield each RNAProfile from the table. """
        yield from self.table.iter_profiles(quantile=self.quantile,
                                            rel=self.rel_name,
                                            order=self.order,
                                            clust=self.clust)

    def iter_states(self):
        """ Yield each RNAState. """
        for profile in self.iter_profiles():
            ct_file = self._get_struct_file(profile.profile)
            try:
                for struct in from_ct(ct_file):
                    yield RNAState.from_struct_profile(struct, profile)
            except FileNotFoundError:
                logger.error(f"Structure file {ct_file} does not exist; please "
                             f"obtain the file, e.g. using seismic fold")


class StructOneTableRunner(OneTableRunner, ABC):

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_struct_file, opt_struct_sect]

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
