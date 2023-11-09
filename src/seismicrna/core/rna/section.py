from functools import cached_property
from logging import getLogger

from .. import path
from ..seq import Section

logger = getLogger(__name__)

IDX_FIELD = "Index"
BASE_FIELD = "Base"
PREV_FIELD = "Prev"
NEXT_FIELD = "Next"
PAIR_FIELD = "Pair"


# RNA Structure classes ################################################


class RnaSection(object):
    """ RNA sequence or section thereof. """

    def __init__(self, *, title: str, section: Section, **kwargs):
        super().__init__(**kwargs)
        self.title = path.fill_whitespace(title)
        self.section = section

    @property
    def ref(self):
        return self.section.ref

    @property
    def end5(self):
        return self.section.end5

    @property
    def end3(self):
        return self.section.end3

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

    def _subsect_kwargs(self, end5: int, end3: int, title: str | None = None):
        return dict(title=(title if title is not None
                           else f"{self.title}__{end5}-{end3}"),
                    section=self.section.subsection(end5, end3))

    def subsection(self, end5: int, end3: int, title: str | None = None):
        return self.__class__(**self._subsect_kwargs(end5, end3, title))

    def __str__(self):
        return f"{type(self).__name__} {repr(self.title)} over {self.section}"

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
