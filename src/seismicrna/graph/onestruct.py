from abc import ABC
from functools import cached_property
from pathlib import Path

from .onetable import OneTableGraph, OneTableRunner
from .rel import OneRelGraph
from ..core import path
from ..core.arg import opt_structs
from ..core.rna import RNAState, from_ct


class StructOneTableGraph(OneTableGraph, OneRelGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *, structs: Path | None = None, **kwargs):
        super().__init__(**kwargs)
        self.structs_file = structs

    @cached_property
    def _struct_fields(self):
        """ Get the fields of the structure. """
        if self.structs_file is not None:
            fields = path.parse(self.structs_file,
                                path.RefSeg,
                                path.SectSeg,
                                path.ConnectTableSeg)
            return fields[path.REF], fields[path.SECT]
        return super().ref, self.sect

    @property
    def ref(self):
        ref, _ = self._struct_fields
        if ref != super().ref:
            raise ValueError(f"Reference names differ between CT file ({ref}) "
                             f"and table file ({super().ref})")
        return ref

    @property
    def struct_sect(self):
        """ Section of the reference that the structure model spans. """
        _, sect = self._struct_fields
        return sect

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
            ct_file = (self.structs_file
                       if self.structs_file is not None
                       else profile.get_ct_file(self.top))
            for struct in from_ct(ct_file):
                yield RNAState.from_struct_profile(struct, profile)


class StructOneTableRunner(OneTableRunner, ABC):

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_structs]

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
