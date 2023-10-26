import re
from abc import ABC

from .file import BrickleIO
from ..batch import MutsBatch, ReadBatch


class ReadBatchIO(ReadBatch, BrickleIO, ABC):
    """ Pickled file of a batch of data. """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z_]*)batchio", cls.__name__.lower()).groups()
        return btype


class MutsBatchIO(MutsBatch, ReadBatchIO, ABC):
    """ Pickled file of a batch of mutational data. """

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
