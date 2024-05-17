import re
from abc import ABC

from .file import BrickleIO
from ..batch import MutsBatch, ReadBatch
from ..seq import Section


class ReadBatchIO(ReadBatch, BrickleIO, ABC):
    """ Pickled file of a batch of data. """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z_]*)batchio", cls.__name__.lower()).groups()
        return btype


class MutsBatchIO(MutsBatch, ReadBatchIO, ABC):
    """ Pickled file of a batch of mutational data. """

    def __init__(self, *args, section: Section, **kwargs):
        super().__init__(*args, **kwargs, section=section, ref=section.ref)

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
