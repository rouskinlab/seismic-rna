from logging import getLogger

from .batch import MaskRelsBatch
from .files import SavedMaskBatch
from .report import MaskReport
from ..core.iodata import BatchedDatasetLoader, BatchedDatasetLinker
from ..relate.files import SavedRelateBatch
from ..relate.load import RelateLoader

logger = getLogger(__name__)

MASK_KEY = "mask-load"


class MaskLoader(BatchedDatasetLoader):
    """ Load batches of masked relation vectors. """

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_batch_type(cls):
        return SavedMaskBatch


class MaskLinker(BatchedDatasetLinker):
    """ Link relation """

    @classmethod
    def get_batch_type(cls):
        return SavedMaskBatch

    @classmethod
    def data1_type(cls):
        return RelateLoader

    @classmethod
    def data2_type(cls):
        return MaskLoader

    def _link(self,
              batch1: SavedRelateBatch,
              batch2: SavedMaskBatch,
              *args,
              **kwargs):
        return MaskRelsBatch.from_batch(batch1,
                                        reads=batch2.read_nums)


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
