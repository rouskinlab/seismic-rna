from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .base import RNARegion
from .db import DB_NAME_MARK, format_db_string, parse_db_string
from .pair import (UNPAIRED,
                   find_root_pairs,
                   pairs_to_dict,
                   pairs_to_table,
                   renumber_pairs,
                   table_to_pairs)
from ..seq import POS_NAME, DNA, RNA, Region, intersect
from ..validate import require_isinstance, require_equal

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

    @classmethod
    def from_db_string(cls,
                       db_string: str,
                       seq: DNA | RNA | str, *,
                       seq5: int = 1,
                       ref: str,
                       reg: str,
                       **kwargs):
        """ Create an RNAStructure from a dot-bracket string. """
        return cls(pairs=parse_db_string(db_string, seq5),
                   region=Region(ref, seq, seq5=seq5, name=reg),
                   **kwargs)

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
        return set(table_to_pairs(self.table))

    @cached_property
    def dict(self):
        """ Map from each paired base to its partner. """
        return pairs_to_dict(self.pairs)

    @cached_property
    def roots(self):
        """ All pairs that are not contained by any other pair. """
        return find_root_pairs(self.pairs)

    @cached_property
    def is_paired(self):
        """ Whether each base is paired. """
        return self.table != UNPAIRED

    @cached_property
    def is_unpaired(self):
        """ Whether each base is unpaired. """
        return self.table == UNPAIRED

    @cached_property
    def is_paired_internally(self):
        """ Whether each base is paired and between two other base pairs
        (no bulges or other unpaired bases next to it). """
        is_internal = self.is_paired.copy()
        # A base can be between two other base pairs only if both bases
        # to its 5' and its 3' are paired.
        paired5 = np.concatenate([[False], is_internal.values[:-1]])
        paired3 = np.concatenate([is_internal.values[1:], [False]])
        is_internal &= (paired5 & paired3)
        # A base can be between two other base pairs only if its partner
        # is also between two base pairs.
        is_internal_index = is_internal[is_internal].index
        is_internal.loc[is_internal_index] = (
                is_internal.loc[is_internal_index]
                &
                is_internal.loc[self.table[is_internal]].values
        )
        # Pairs between adjacent bases are not considered internal
        # because they are not flanked by a base pair on both sides.
        table_is_internal = self.table[is_internal]
        distance = np.abs(
            table_is_internal
            -
            table_is_internal.index.get_level_values(POS_NAME)
        )
        is_internal.loc[distance.index] = distance > 1
        # To handle pseudoknots, a base can be between two other base
        # pairs only if its 5' base and partner's 3' base are partners,
        # and its 3' base and partner's 5' base are partners.
        is_internal_pos = is_internal[is_internal].index.get_level_values(
            POS_NAME
        )
        internal_partners = self.table.loc[is_internal_pos]
        pseudoknot = (np.logical_or(
            internal_partners - 1 != self.table.loc[is_internal_pos + 1].values,
            internal_partners + 1 != self.table.loc[is_internal_pos - 1].values
        ))
        is_internal.loc[pseudoknot.loc[pseudoknot].index] = False
        return is_internal

    @cached_property
    def is_paired_terminally(self):
        """ Whether each base is paired and terminates a consecutive
        stretch of base pairs (i.e. is not internally paired). """
        return self.is_paired & ~self.is_paired_internally

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
    def db_string(self):
        """ Dot-bracket string (structure only). """
        return format_db_string(self.pairs,
                                self.region.length,
                                self.region.end5)

    def get_db_text(self, sequence: bool):
        """ Dot-bracket record. """
        lines = [self.db_title]
        if sequence:
            lines.append(str(self.seq))
        lines.append(self.db_string)
        return "".join([f"{line}\n" for line in lines])


def calc_wfmi(struct1: RNAStructure,
              struct2: RNAStructure,
              terminal_pairs: bool = True):
    """ Weighted Fowlkes-Mallows index between two structures. """
    require_isinstance("struct1", struct1, RNAStructure)
    require_isinstance("struct2", struct2, RNAStructure)
    region = struct1.region
    require_equal("struct1.region",
                  region,
                  struct2.region,
                  "struct2.region")
    # Whether each base is paired in each structure.
    if terminal_pairs:
        is_paired1 = struct1.is_paired
        is_paired2 = struct2.is_paired
    else:
        is_paired1 = struct1.is_paired_internally
        is_paired2 = struct2.is_paired_internally
    # Number of positions that are unpaired in both structures.
    n_unpaired = np.count_nonzero(struct1.is_unpaired & struct2.is_unpaired)
    # Number of positions with structural information.
    n_total = n_unpaired + np.count_nonzero(is_paired1 | is_paired2)
    if n_total == 0:
        # There are no positions with structural information.
        return np.nan
    # Calculate the Fowlkes-Mallows index (FMI).
    n_paired1 = np.count_nonzero(is_paired1)
    n_paired2 = np.count_nonzero(is_paired2)
    if n_paired1 > 0 and n_paired2 > 0:
        # FMI = TP / √(TP+FP)*(TP+FN)
        both_paired = is_paired1 & is_paired2
        n_paired12 = np.count_nonzero(
            struct1.table[both_paired] == struct2.table[both_paired]
        )
        fmi = n_paired12 / np.sqrt(n_paired1 * n_paired2)
    else:
        # At least one structure has 0 base pairs.
        fmi = 0.
    # Weight the FMI by the fraction of positions that are paired in at
    # least one structure.
    f_unpaired = n_unpaired / n_total
    return float(f_unpaired + (1. - f_unpaired) * fmi)
