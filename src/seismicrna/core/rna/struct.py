from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .base import RNARegion
from .db import DB_NAME_MARK, format_db_structure
from .pair import (UNPAIRED,
                   find_root_pairs,
                   pairs_to_dict,
                   pairs_to_table,
                   renumber_pairs,
                   table_to_pairs)
from ..seq import POS_NAME, intersect

IDX_FIELD = "Index"
BASE_FIELD = "Base"
PREV_FIELD = "Prev"
NEXT_FIELD = "Next"
PAIR_FIELD = "Pair"


class Rna2dPart(object):
    """ Part of an RNA secondary structure. """

    def __init__(self, *regions: RNARegion, **kwargs):
        super().__init__(**kwargs)
        if inter := intersect(*[region.region for region in regions]):
            raise ValueError(f"Regions intersect: {inter}")
        self.regions = regions


class Rna2dStem(Rna2dPart):
    """ An RNA stem (contiguous double helix). """

    def __init__(self, side1: RNARegion, side2: RNARegion, **kwargs):
        # The lengths of the two regions must equal.
        if side1.region.length != side2.region.length:
            raise ValueError(f"The lengths of side 1 ({side1.region.length}) "
                             f"and side 2 ({side2.region.length}) are unequal")
        # Determine the order of the sides.
        if side1.region.end3 < side2.region.end5:
            region5, region3 = side1, side2
        else:
            region5, region3 = side2, side1
        super().__init__(region5, region3, **kwargs)

    @property
    def region5(self):
        region5, _ = self.regions
        return region5

    @property
    def region3(self):
        _, region3 = self.regions
        return region3


class RnaJunction(Rna2dPart):
    """ A junction between stems in an RNA structure. """


class Rna2dStemLoop(RnaJunction):
    """ An RNA loop at the end of a stem. """

    def __init__(self, region: RNARegion, **kwargs):
        super().__init__(region, **kwargs)

    @property
    def region(self):
        region, = self.regions
        return region


class RNAStructure(RNARegion):
    """ Secondary structure of an RNA. """

    def __init__(self, *,
                 title: str,
                 pairs: Iterable[tuple[int, int]],
                 branch: str = "",
                 **kwargs):
        """
        Parameters
        ----------
        title: str
            Title of the RNA structure, as written in the CT/DB file.
        pairs: Iterable[tuple[int, int]]
            Base pairs in the structure.
        branch: str
            Branch of the workflow for folding (optional).
        """
        super().__init__(**kwargs)
        self.title = title
        self.table = pairs_to_table(pairs, self.region)
        self.branch = branch

    @cached_property
    def init_args(self):
        return super().init_args | dict(title=self.title,
                                        pairs=self.pairs,
                                        branch=self.branch)

    def _renumber_from_args(self, seq5: int):
        return super()._renumber_from_args(seq5) | dict(
            pairs=renumber_pairs(self.pairs, seq5 - self.region.end5)
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
        return self.table != UNPAIRED

    def _subregion_kwargs(self,
                          end5: int,
                          end3: int,
                          title: str | None = None):
        return super()._subregion_kwargs(end5, end3) | dict(
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
            yield self.subregion(end5, end3)

    @property
    def ct_title(self):
        """ Header line for the CT file."""
        return f"{self.region.length}\t{self.title}"

    @cached_property
    def ct_data(self):
        """ Convert the connectivity table to a DataFrame. """
        # Make an index the same length as the region and starting
        # from 1 (CT files must start at index 1).
        index = pd.Index(self.region.range_one, name=IDX_FIELD)
        # Adjust the numbers of the paired bases (i.e. pairs > 0) such
        # that they also index from 1.
        pairs = self.table.values.copy()
        pairs[pairs > UNPAIRED] -= self.region.end5 - 1
        # Generate the data for the connectivity table.
        data = {BASE_FIELD: self.seq.array,
                PREV_FIELD: index.values - 1,
                NEXT_FIELD: index.values + 1,
                PAIR_FIELD: pairs,
                POS_NAME: self.region.range_int}
        # The last value of the next field must be 0.
        if index.size > 0:
            data[NEXT_FIELD][-1] = 0
        # Assemble the data into a DataFrame.
        return pd.DataFrame.from_dict({field: pd.Series(values, index=index)
                                       for field, values in data.items()})

    @property
    def ct_text(self):
        """ Connectivity table as text. """
        data = self.ct_data.reset_index()
        return f"{self.ct_title}\n{data.to_string(index=False, header=False)}\n"

    @property
    def db_title(self):
        """ Header line for the DB file. """
        return f"{DB_NAME_MARK}{self.title}"

    @cached_property
    def db_structure(self):
        """ Dot-bracket string (structure only). """
        return format_db_structure(self.pairs,
                                   self.region.length,
                                   self.region.end5)

    def get_db_text(self, sequence: bool):
        """ Dot-bracket record. """
        lines = [self.db_title]
        if sequence:
            lines.append(str(self.seq))
        lines.append(self.db_structure)
        return "".join([f"{line}\n" for line in lines])
