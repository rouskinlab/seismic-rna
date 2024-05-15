from .report import PoolReport
from ..core.data import LoadFunction, TallMutsDataset
from ..relate.data import RelateDataset


class PoolDataset(TallMutsDataset):
    """ Load pooled batches of relation vectors. """

    @classmethod
    def get_report_type(cls):
        return PoolReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_relate_dataset


load_relate_dataset = LoadFunction(RelateDataset, PoolDataset)

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
