from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .pair import (find_root_pairs,
                   pairs_to_dict,
                   pairs_to_table,
                   renumber_pairs,
                   table_to_pairs)
from .base import RNASection
from ..seq import POS_NAME, intersection

IDX_FIELD = "Index"
BASE_FIELD = "Base"
PREV_FIELD = "Prev"
NEXT_FIELD = "Next"
PAIR_FIELD = "Pair"


class Rna2dPart(object):
    """ Part of an RNA secondary structure. """

    def __init__(self, *sections: RNASection, **kwargs):
        super().__init__(**kwargs)
        if inter := intersection(*[section.section for section in sections]):
            raise ValueError(f"Sections intersect: {inter}")
        self.sections = sections


class Rna2dStem(Rna2dPart):
    """ An RNA stem (contiguous double helix). """

    def __init__(self, side1: RNASection, side2: RNASection, **kwargs):
        # The lengths of the two sections must equal.
        if side1.section.length != side2.section.length:
            raise ValueError(f"The lengths of side 1 ({side1.section.length}) "
                             f"and side 2 ({side2.section.length}) are unequal")
        # Determine the order of the sides.
        if side1.section.end3 < side2.section.end5:
            section5, section3 = side1, side2
        else:
            section5, section3 = side2, side1
        super().__init__(section5, section3, **kwargs)

    @property
    def section5(self):
        section5, _ = self.sections
        return section5

    @property
    def section3(self):
        _, section3 = self.sections
        return section3


class RnaJunction(Rna2dPart):
    """ A junction between stems in an RNA structure. """


class Rna2dStemLoop(RnaJunction):
    """ An RNA loop at the end of a stem. """

    def __init__(self, section: RNASection, **kwargs):
        super().__init__(section, **kwargs)

    @property
    def section(self):
        section, = self.sections
        return section


class RNAStructure(RNASection):
    """ Secondary structure of an RNA. """

    def __init__(self, *,
                 title: str,
                 pairs: Iterable[tuple[int, int]],
                 **kwargs):
        """
        Parameters
        ----------
        title: str
            Title of the RNA structure, as written in the CT/DB file.
        pairs: Iterable[tuple[int, int]]
            Base pairs in the structure.
        """
        super().__init__(**kwargs)
        self.title = title
        self.table = pairs_to_table(pairs, self.section)

    @cached_property
    def init_args(self):
        return super().init_args | dict(title=self.title, pairs=self.pairs)

    def _renumber_from_args(self, seq5: int):
        return super()._renumber_from_args(seq5) | dict(
            pairs=renumber_pairs(self.pairs, seq5 - self.section.end5)
        )

    @cached_property
    def pairs(self):
        """ Base pairs in the structure. """
        return tuple(table_to_pairs(self.table))

    @cached_property
    def dict(self):
        return pairs_to_dict(self.pairs)

    @cached_property
    def roots(self):
        return find_root_pairs(self.pairs)

    @cached_property
    def is_paired(self):
        """ Series where each index is a position and each value is True
        if the corresponding base is paired, otherwise False. """
        return self.table != 0

    def _subsection_kwargs(self,
                           end5: int,
                           end3: int,
                           title: str | None = None):
        return super()._subsection_kwargs(end5, end3) | dict(
            title=(title if title is not None
                   else f"{self.title}_{end5}-{end3}"),
            pairs=table_to_pairs(
                self.table[np.logical_and(
                    self.table.index.get_level_values(POS_NAME) >= end5,
                    self.table.index.get_level_values(POS_NAME) <= end3
                )]
            )
        )

    def iter_root_modules(self):
        for end5, end3 in self.roots:
            yield self.subsection(end5, end3)

    @property
    def ct_title(self):
        """ Header line for the CT file."""
        return f"{self.section.length}\t{self.title}"

    @cached_property
    def ct_data(self):
        """ Return the connectivity table as a DataFrame. """
        # Make an index the same length as the section and starting
        # from 1 (CT files must start at index 1).
        index = pd.Index(self.section.range_one, name=IDX_FIELD)
        # Adjust the numbers of the paired bases (i.e. pairs > 0) such
        # that they also index from 1.
        pairs = self.table.values.copy()
        pairs[pairs > 0] -= self.section.end5 - 1
        # Generate the data for the connectivity table.
        data = {BASE_FIELD: self.seq.array,
                PREV_FIELD: index.values - 1,
                NEXT_FIELD: index.values + 1,
                PAIR_FIELD: pairs,
                POS_NAME: self.section.range_int}
        # The last value of the next field must be 0.
        if index.size > 0:
            data[NEXT_FIELD][-1] = 0
        # Assemble the data into a DataFrame.
        return pd.DataFrame.from_dict({field: pd.Series(values, index=index)
                                       for field, values in data.items()})

    @property
    def ct_text(self):
        """ Return the connectivity table as text. """
        data = self.ct_data.reset_index()
        return f"{self.ct_title}\n{data.to_string(index=False, header=False)}\n"

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
