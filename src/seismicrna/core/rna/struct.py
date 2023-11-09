from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .ct import parse_ct
from .pair import find_root_pairs, pairs_to_dict, pairs_to_table, table_to_pairs
from .section import (BASE_FIELD,
                      IDX_FIELD,
                      NEXT_FIELD,
                      PAIR_FIELD,
                      PREV_FIELD,
                      RnaSection)
from ..seq import intersection, POS_NAME, Section

logger = getLogger(__name__)


class Rna2dPart(object):
    """ Part of an RNA secondary structure. """

    def __init__(self, *sections: RnaSection, **kwargs):
        super().__init__(**kwargs)
        if inter := intersection(*[section.section for section in sections]):
            raise ValueError(f"Sections intersect: {inter}")
        self.sections = sections


class Rna2dStem(Rna2dPart):
    """ An RNA stem (contiguous double helix). """

    def __init__(self, side1: RnaSection, side2: RnaSection, **kwargs):
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

    def __init__(self, section: RnaSection, **kwargs):
        super().__init__(section, **kwargs)

    @property
    def section(self):
        section, = self.sections
        return section


class Rna2dStructure(RnaSection):
    """ RNA secondary structure. """

    def __init__(self, *, pairs: Iterable[tuple[int, int]], **kwargs):
        super().__init__(**kwargs)
        self.table = pairs_to_table(pairs, self.section)

    @cached_property
    def pairs(self):
        return tuple(table_to_pairs(self.table))

    @cached_property
    def dict(self):
        return pairs_to_dict(self.pairs)

    @cached_property
    def roots(self):
        return find_root_pairs(self.pairs)

    def _subsect_kwargs(self, end5: int, end3: int, title: str | None = None):
        return super()._subsect_kwargs(end5, end3, title) | dict(
            pairs=table_to_pairs(
                self.table[np.logical_and(self.table.index.values >= end5,
                                          self.table.index.values <= end3)]
            )
        )

    def iter_root_modules(self):
        for end5, end3 in self.roots:
            yield self.subsection(end5, end3)

    @property
    def header(self):
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
        data = {
            BASE_FIELD: self.seq.to_str_array(),
            PREV_FIELD: index.values - 1,
            NEXT_FIELD: index.values + 1,
            PAIR_FIELD: pairs,
            POS_NAME: self.section.range_int,
        }
        # Assemble the data into a DataFrame.
        return pd.DataFrame.from_dict({field: pd.Series(values, index=index)
                                       for field, values in data.items()})

    @property
    def ct_text(self):
        """ Return the connectivity table as text. """
        data = self.ct_data.reset_index()
        return f"{self.header}\n{data.to_string(index=False, header=False)}\n"


def from_ct(ct_file: Path, section: Section):
    section_rna_seq = section.seq.tr()
    for title, seq, pairs in parse_ct(ct_file, section.end5):
        if seq == section_rna_seq:
            yield Rna2dStructure(title=title, section=section, pairs=pairs)
        else:
            logger.error(f"Expected {section_rna_seq}, but got {seq}")

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
