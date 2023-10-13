from logging import getLogger

from .files import SavedRelateBatch
from .report import RelateReport
from ..core.iodata import BatchedDatasetLoader

logger = getLogger(__name__)

POS = "by_position"
VEC = "by_vector"
READ = "Read Name"


class RelateLoader(BatchedDatasetLoader):
    """ Load batches of relation vectors. """

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @classmethod
    def get_batch_type(cls):
        return SavedRelateBatch

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
