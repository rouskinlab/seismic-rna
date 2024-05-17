from .batch import RelateBatch
from .io import RelateBatchIO
from .report import RelateReport
from ..core.data import LoadedMutsDataset


class RelateDataset(LoadedMutsDataset):
    """ Dataset from the Relate step. """

    @classmethod
    def get_batch_type(cls):
        return RelateBatchIO

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @property
    def pattern(self):
        return None

    def get_batch(self, batch: int):
        relate_batch = super().get_batch(batch)
        return RelateBatch(batch=relate_batch.batch,
                           seg_end5s=relate_batch.seg_end5s,
                           seg_end3s=relate_batch.seg_end3s,
                           muts=relate_batch.muts,
                           section=self.section,
                           sanitize=False)

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
