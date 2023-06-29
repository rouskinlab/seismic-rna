from functools import cache, cached_property
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from . import path
from .sect import Section
from .seq import write_fasta

logger = getLogger(__name__)


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
        """ Sequence in RNA format. """
        return self.section.seq.tr()

    @property
    def seq_record(self):
        return self.section.ref_sect, self.seq


class RnaProfile(RnaSection):
    """ Mutational profile of an RNA from a specific sample. """

    def __init__(self, title: str, section: Section,
                 sample: str, data_sect: str, reacts: pd.Series):
        super().__init__(title, section)
        self.sample = sample
        self.dms_sect = data_sect
        if np.any(reacts < 0.) or np.any(reacts > 1.):
            raise ValueError(f"Got reactivities outside [0, 1]: {reacts}")
        self.reacts = reacts.reindex(self.section.range)

    def get_ceiling(self, quantile: float) -> float:
        """ Compute the maximum reactivity given a quantile. """
        if not 0. < quantile <= 1.:
            raise ValueError("Quantile for reactivity ceiling must be in "
                             f"(0, 1], but got {quantile}")
        return np.nanquantile(self.reacts.values, quantile)

    def normalize(self, quantile: float):
        """ Normalize the reactivities to a given quantile. """
        return self.reacts / self.get_ceiling(quantile)

    def winsorize(self, quantile: float):
        """ Normalize and winsorize the reactivities. """
        return pd.Series(np.clip(self.normalize(quantile), 0., 1.),
                         index=self.reacts.index)

    @cache
    def get_dir(self, out_dir: Path):
        """ Get the directory in which this RNA's files will go. """
        return path.builddir(path.ModSeg, path.SampSeg, path.RefSeg,
                             path.SectSeg, path.FoldSectSeg,
                             top=out_dir, module=path.MOD_STRUCT,
                             sample=self.sample, ref=self.ref,
                             sect=self.dms_sect, fold_sect=self.sect)

    @cache
    def get_file(self, out_dir: Path, segment: path.Segment, **kwargs):
        """ Get the path to a file of the RNA sequence. """
        return self.get_dir(out_dir).joinpath(segment.build(**kwargs))

    def get_fasta(self, out_dir: Path):
        """ Get the path to the FASTA file of the RNA sequence. """
        return self.get_file(out_dir, path.FastaSeg,
                             ref=self.title, ext=path.FASTA_EXTS[0])

    def ct_file(self, out_dir: Path):
        """ Get the path to the CT file of the RNA. """
        return self.get_file(out_dir, path.ConnectTableSeg,
                             struct=self.title, ext=path.CT_EXT)

    def dot_file(self, out_dir: Path):
        """ Get the path to the DOT file of the RNA. """
        return self.get_file(out_dir, path.DotBracketSeg,
                             struct=self.title, ext=path.DOT_EXTS[0])

    def dms_file(self, out_dir: Path):
        """ Get the path to the DMS data file of the RNA. """
        return self.get_file(out_dir, path.DmsReactsSeg,
                             reacts=self.title, ext=path.DMS_EXT)

    def varnac_file(self, out_dir: Path):
        """ Get the path to the VARNA color file of the RNA. """
        return self.get_file(out_dir, path.VarnaColorSeg,
                             reacts=self.title, ext=path.TXT_EXT)

    def to_fasta(self, out_dir: Path):
        """ Write the RNA sequence to a FASTA file. """
        fasta = self.get_fasta(out_dir)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_dms(self, out_dir: Path, quantile: float):
        """ Write the DMS reactivities to a DMS file. """
        # The DMS reactivities must be numbered starting from 1 at the
        # beginning of the section, even if the section does not start
        # at 1. Renumber the section from 1.
        dms = pd.Series(self.winsorize(quantile).values,
                        index=self.section.range_int_one)
        # Drop bases with missing data to make RNAstructure ignore them.
        dms.dropna(inplace=True)
        # Write the DMS reactivities to the DMS file.
        dms_file = self.dms_file(out_dir)
        dms.to_csv(dms_file, sep="\t", header=False)
        return dms_file

    def to_varnac(self, out_dir: Path, quantile: float):
        """ Write the VARNA colors to a file. """
        # Normalize and winsorize the DMS reactivities.
        varnac = self.winsorize(quantile)
        # Fill missing values with -1, to signify no data.
        varnac.fillna(-1., inplace=True)
        # Write the values to the VARNA color file.
        varnac_file = self.varnac_file(out_dir)
        varnac.to_csv(varnac_file, header=False, index=False)
        return varnac_file


class Rna2dStructure(RnaSection):
    """ RNA secondary structure. """

    IDX_FIELD = "n"
    BASE_FIELD = "Base"
    PRIOR_FIELD = "n-1"
    NEXT_FIELD = "n+1"
    PAIR_FIELD = "Pair"
    POS_FIELD = "Position"

    def __init__(self,
                 title: str,
                 section: Section,
                 pairs: Iterable[tuple[int, int]]):
        super().__init__(title, section)
        self.title = title
        self.pairs = set(pairs)

    @property
    def header(self):
        return f"{self.section.length}\t{self.title}"

    @cached_property
    def partners(self):
        """ Return a Series of every position in the section and the
        position to which it pairs, or 0 if it does not pair. """
        # Initialize the series of pairs to 0 for every position.
        partners = pd.Series(0, index=pd.Index(self.section.range,
                                               name=self.POS_FIELD))

        def add_pair(at: int, to: int):
            """ Add a base pair at position `at` to position `to`. """
            # Find the current pairing partner at this position.
            try:
                to2 = partners.loc[at]
            except KeyError:
                raise ValueError(f"Position {at} is not in {self.section}")
            if to2 == 0:
                # There is no pairing partner at this position: add it.
                partners.loc[at] = to
            elif to2 == to:
                # A previous partner matches the current partner.
                logger.warning(f"Pair {at}-{to} was given twice")
            else:
                # A previous partner conflicts with the current partner.
                raise ValueError(f"Position {at} was given pairs with both "
                                 f"{to2} and {to}")

        # Add all base pairs (in both directions) to the table.
        for pos1, pos2 in self.pairs:
            add_pair(pos1, pos2)
            add_pair(pos2, pos1)
        return partners

    @property
    def ct_data(self):
        """ Return the connectivity table as a DataFrame. """
        # Make an index the same length as the section and starting
        # from 1 (CT files must start at index 1).
        index = pd.RangeIndex(1, self.section.length + 1, name=self.IDX_FIELD)
        # Adjust the numbers of the paired bases (i.e. pairs > 0) such
        # that they also index from 1.
        pairs = self.partners.values.copy()
        pairs[pairs > 0] -= self.section.end5 - 1
        # Generate the data for the connectivity table.
        data = {
            self.BASE_FIELD: self.seq.to_str_array(),
            self.PRIOR_FIELD: index.values - 1,
            self.NEXT_FIELD: index.values + 1,
            self.PAIR_FIELD: pairs,
            self.POS_FIELD: self.section.range,
        }
        # Assemble the data into a DataFrame.
        return pd.DataFrame.from_dict({field: pd.Series(values, index=index)
                                       for field, values in data.items()})

    @property
    def ct_text(self):
        """ Return the connectivity table as text. """
        data = self.ct_data.reset_index()
        return f"{self.header}\n{data.to_string(index=False, header=False)}\n"


class RnaState(Rna2dStructure):
    """ RNA sequence and mutation rates. """

    def __init__(self,
                 title: str,
                 section: Section,
                 pairs: Iterable[tuple[int, int]],
                 sample: str,
                 mus: pd.Series):
        """
        Parameters
        ----------
        mus: pd.DataFrame
            Mutation rates of the RNA molecule.
        """
        super().__init__(title, section, pairs)
        self.sample = sample
        self.mus = mus


class RnaEnsemble(RnaSection):
    # FIXME: write this
    pass


def parse_ct_file(ct_file: Path):
    """ Yield RNA secondary structures from a CT file. """
    # FIXME: write this
    return
