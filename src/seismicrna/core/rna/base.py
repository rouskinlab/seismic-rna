from functools import cached_property
from logging import getLogger

from ..seq import Section

logger = getLogger(__name__)


class RNASection(object):
    """ Section of an RNA sequence. """

    def __init__(self, *, section: Section, **kwargs):
        super().__init__(**kwargs)
        self.section = section

    @property
    def init_args(self):
        """ Arguments needed to initialize a new instance. """
        return dict(section=self.section)

    @property
    def ref(self):
        """ Name of the reference sequence. """
        return self.section.ref

    @property
    def end5(self):
        """ Position of the 5' end of the section. """
        return self.section.end5

    @property
    def end3(self):
        """ Position of the 3' end of the section. """
        return self.section.end3

    @property
    def sect(self):
        """ Name of the section. """
        return self.section.name

    @cached_property
    def seq(self):
        """ Sequence of the section as RNA. """
        return self.section.seq.tr()

    @property
    def seq_record(self):
        return self.section.ref_sect, self.seq

    def _subsection_kwargs(self, end5: int, end3: int):
        """ Keyword arguments used by self.subsection(). """
        return dict(section=self.section.subsection(end5, end3))

    def subsection(self, end5: int, end3: int):
        return self.__class__(**self._subsection_kwargs(end5, end3))

    def _renumber_from_args(self, seq5: int):
        """ Arguments needed to initialize a renumbered instance. """
        return self.init_args | dict(section=self.section.renumber_from(seq5))

    def renumber_from(self, seq5: int):
        """ Return a new RNASection renumbered starting from a position.

        Parameters
        ----------
        seq5: int
            Position from which to start the new numbering system.

        Returns
        -------
        RNASection
            Section with renumbered positions.
        """
        return self.__class__(**self._renumber_from_args(seq5))

    def __str__(self):
        return f"{type(self).__name__} over {self.section}"

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
