from .io import MaskReadBatchIO
from .report import MaskReport
from ..core.io import BatchedDatasetLoader, BatchedDatasetLinker
from ..relate.io import RelateBatchIO
from ..relate.data import RelateLoader


class MaskLoader(BatchedDatasetLoader):
    """ Load batches of masked relation vectors. """

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_batch_type(cls):
        return MaskReadBatchIO


class MaskLinker(BatchedDatasetLinker):
    """ Link relation """

    @classmethod
    def get_batch_type(cls):
        return MaskReadBatchIO

    @classmethod
    def get_data1_type(cls):
        return RelateLoader

    @classmethod
    def get_data2_type(cls):
        return MaskLoader

    def _link(self,
              batch1: RelateBatchIO,
              batch2: MaskReadBatchIO,
              *args,
              **kwargs):
        return batch1.mask(reads=batch2.read_nums)


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
