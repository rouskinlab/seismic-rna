from abc import ABC
from functools import cached_property

from .base import GraphBase
from ..table.base import get_rel_name


class OneRelGraph(GraphBase, ABC):
    """ Graph of exactly one type of relationship. """

    def __init__(self, *, rel: str, **kwargs):
        """
        Parameters
        ----------
        rel: str
            Relationship(s) whose data will be pulled from the table(s).
            It is given as a one-letter code, the definitions of which
            are given in seismicrna.table.base.REL_CODES, as follows:

            - Covered: v
            - Informed: n
            - Matched: r
            - Mutated: m
            - Subbed: s
            - Subbed to A: a
            - Subbed to C: c
            - Subbed to G: g
            - Subbed to T: t
            - Deleted: d
            - Inserted: i

            This type of graph requires that `rel` be a one-letter code.
        """
        super().__init__(**kwargs)
        if len(rel) != 1:
            raise ValueError(f"{type(self).__name__} expected exactly one "
                             f"relationship code, but got {repr(rel)}")
        self._rel = rel

    @property
    def codestring(self):
        return self._rel

    @cached_property
    def rel_name(self):
        """ Name of the relationship to graph. """
        return get_rel_name(self.codestring)

    @cached_property
    def rel_names(self):
        return [self.rel_name]


class MultiRelsGraph(GraphBase, ABC):
    """ Graph of one or more relationships. """

    def __init__(self, *, rels: str, **kwargs):
        """
        Parameters
        ----------
        rels: str
            Relationships whose data will be pulled from the table(s).
            Each is given as a one-letter code, the definitions of which
            are given in seismicrna.table.base.REL_CODES, as follows:

            - Covered: v
            - Informed: n
            - Matched: r
            - Mutated: m
            - Subbed: s
            - Subbed to A: a
            - Subbed to C: c
            - Subbed to G: g
            - Subbed to T: t
            - Deleted: d
            - Inserted: i

            More than one relationship can be specified by passing a
            string longer than one character; for example, `acgt` would
            mean to use all four types of substitutions.
        """
        super().__init__(**kwargs)
        if len(rels) == 0:
            raise ValueError(f"{type(self).__name__} expected one or more "
                             f"relationship codes, but got {repr(rels)}")
        self._rels = rels

    @property
    def codestring(self):
        return self._rels

    @cached_property
    def rel_names(self):
        return list(map(get_rel_name, self.codestring))

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
