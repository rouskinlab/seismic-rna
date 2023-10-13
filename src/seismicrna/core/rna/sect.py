from functools import cached_property
from logging import getLogger

from .. import path
from ..sect import Section

logger = getLogger(__name__)

IDX_FIELD = "Index"
BASE_FIELD = "Base"
PREV_FIELD = "Prev"
NEXT_FIELD = "Next"
PAIR_FIELD = "Pair"


# RNA Structure classes ################################################


class RnaSection(object):
    """ RNA sequence or section thereof. """

    def __init__(self, title: str, section: Section):
        self.title = path.fill_whitespace(title)
        self.section = section

    @property
    def ref(self):
        return self.section.ref

    @property
    def sect(self):
        return self.section.name

    @cached_property
    def seq(self):
        """ Sequence as RNA. """
        return self.section.seq.tr()

    @property
    def seq_record(self):
        return self.section.ref_sect, self.seq

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
