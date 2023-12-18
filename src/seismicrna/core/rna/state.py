from functools import cached_property

from .profile import RNAProfile
from .roc import compute_auc_roc, compute_roc_curve
from .struct import RNAStructure


class RNAState(RNAStructure, RNAProfile):
    """ RNA secondary structure with mutation rates. """

    @classmethod
    def from_struct_profile(cls, struct: RNAStructure, profile: RNAProfile):
        """ Make an RNAState from an RNAStructure and an RNAProfile. """
        return cls(section=struct.section,
                   title=struct.title,
                   pairs=struct.pairs,
                   sample=profile.sample,
                   data_sect=profile.data_sect,
                   data_name=profile.data_name,
                   data=profile.data)

    @cached_property
    def roc(self):
        return compute_roc_curve(self.table != 0, self.data)

    @cached_property
    def auc(self):
        return compute_auc_roc(*self.roc)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
