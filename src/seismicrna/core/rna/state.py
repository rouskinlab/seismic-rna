from functools import cached_property

from .profile import RNAProfile
from .roc import compute_auc, compute_roc_curve, compute_rolling_auc
from .struct import RNAStructure


class RNAState(RNAStructure, RNAProfile):
    """ RNA secondary structure with mutation rates. """

    @classmethod
    def from_struct_profile(cls, struct: RNAStructure, profile: RNAProfile):
        """ Make an RNAState from an RNAStructure and an RNAProfile. """
        if struct.section.ref != profile.section.ref:
            raise ValueError("Reference names differ between "
                             f"structure ({repr(struct.section.ref)}) "
                             f"and profile ({repr(profile.section.ref)})")
        return cls(section=struct.section,
                   title=struct.title,
                   pairs=struct.pairs,
                   sample=profile.sample,
                   data_sect=profile.data_sect,
                   data_name=profile.data_name,
                   data=profile.data)

    @cached_property
    def roc(self):
        return compute_roc_curve(self.is_paired, self.data)

    @cached_property
    def auc(self):
        return compute_auc(*self.roc)

    def rolling_auc(self, size: int, min_data: int = 2):
        return compute_rolling_auc(self.is_paired, self.data, size, min_data)


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
