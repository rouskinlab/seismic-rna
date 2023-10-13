from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .ct import parse_ct_pairs
from .pair import find_root_pairs, pairs_to_dict, pairs_to_table, table_to_pairs
from .sect import (BASE_FIELD,
                   IDX_FIELD,
                   NEXT_FIELD,
                   PAIR_FIELD,
                   PREV_FIELD,
                   RnaSection)
from ..sect import intersection, POS_NAME, Section


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

    def __init__(self, pairs: Iterable[tuple[int, int]], **kwargs):
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


def parse_ct_structures(ct_file: Path, section: Section):
    section_rna_seq = section.seq.tr()
    for title, seq, pairs in parse_ct_pairs(ct_file, section.end5):
        if seq != section_rna_seq:
            raise ValueError(f"Expected {section_rna_seq}, but got {seq}")
        yield Rna2dStructure(title=title, section=section, pairs=pairs)
