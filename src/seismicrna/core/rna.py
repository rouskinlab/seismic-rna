from functools import cache, cached_property
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from . import path
from .sect import Section, POS_NAME
from .seq import write_fasta
from .pair import pairs_to_partners, parse_ct_pairs

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


class RnaProfile(RnaSection):
    """ Mutational profile of an RNA from a specific sample. """

    def __init__(self, *args, sample: str, data_sect: str, reacts: pd.Series,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.sample = sample
        self.data_sect = data_sect
        if np.any(reacts < 0.) or np.any(reacts > 1.):
            raise ValueError(f"Got reactivities outside [0, 1]: {reacts}")
        self.reacts = reacts.reindex(self.section.range)

    @cache
    def get_dir(self, out_dir: Path):
        """ Get the directory in which this RNA's files will go. """
        return path.builddir(path.ModSeg, path.SampSeg, path.RefSeg,
                             path.SectSeg, path.FoldSectSeg,
                             top=out_dir, module=path.MOD_STRUCT,
                             sample=self.sample, ref=self.ref,
                             sect=self.data_sect, fold_sect=self.sect)

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

    def varna_color_file(self, out_dir: Path):
        """ Get the path to the VARNA color file of the RNA. """
        return self.get_file(out_dir, path.VarnaColorSeg,
                             reacts=self.title, ext=path.TXT_EXT)

    def to_fasta(self, out_dir: Path):
        """ Write the RNA sequence to a FASTA file. """
        fasta = self.get_fasta(out_dir)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_dms(self, out_dir: Path):
        """ Write the DMS reactivities to a DMS file. """
        # The DMS reactivities must be numbered starting from 1 at the
        # beginning of the section, even if the section does not start
        # at 1. Renumber the section from 1.
        dms = self.reacts.copy()
        dms.index = self.section.range_one
        # Drop bases with missing data to make RNAstructure ignore them.
        dms.dropna(inplace=True)
        # Write the DMS reactivities to the DMS file.
        dms_file = self.dms_file(out_dir)
        dms.to_csv(dms_file, sep="\t", header=False)
        return dms_file

    def to_varna_color_file(self, out_dir: Path):
        """ Write the VARNA colors to a file. """
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.reacts.fillna(-1.)
        # Write the values to the VARNA color file.
        varna_color_file = self.varna_color_file(out_dir)
        varna_color.to_csv(varna_color_file, header=False, index=False)
        return varna_color_file


class Rna2dStructure(RnaSection):
    """ RNA secondary structure. """

    def __init__(self, *args, pairs: Iterable[tuple[int, int]], **kwargs):
        super().__init__(*args, **kwargs)
        self.partners = pairs_to_partners(pairs, self.section)

    @cached_property
    def pairs(self) -> list[tuple[int, int]]:
        """ List of tuples of the 5' and 3' position in each pair. """
        order_5to3 = np.logical_and(self.partners != 0,
                                    self.partners > self.partners.index)
        return list(self.partners.loc[order_5to3].items())

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
        pairs = self.partners.values.copy()
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


class RnaState(Rna2dStructure, RnaProfile):
    """ RNA secondary structure with mutation rates. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @cached_property
    def roc(self):
        return compute_roc_curve(self.partners != 0, self.reacts)

    @cached_property
    def auc(self):
        fpr, tpr = self.roc
        return compute_auc_roc(fpr, tpr)


# Helper functions #####################################################


def parse_ct_structures(ct_file: Path, section: Section):
    section_rna_seq = section.seq.tr()
    for title, seq, pairs in parse_ct_pairs(ct_file, section.end5):
        if seq != section_rna_seq:
            raise ValueError(f"Expected {section_rna_seq}, but got {seq}")
        yield Rna2dStructure(title=title, section=section, pairs=pairs)


def compute_roc_curve(paired: pd.Series, reacts: pd.Series):
    """ Compute the receiver operating characteristic (ROC) curve to
    compute how well chemical reactivities agree with a structure. """
    # Use only positions with non-missing reactivities.
    reacts_use = reacts.loc[np.logical_not(np.isnan(reacts.values))]
    # Count the number of positions with non-missing reactivities.
    n_use = reacts_use.size
    # Sort the positions in ascending order by reactivity.
    order = reacts_use.sort_values().index.get_level_values(POS_NAME)
    paired_use = paired.loc[order].astype(bool, copy=False)
    # Compute the cumulative number of paired bases at each threshold.
    paired_cumsum = np.hstack([[0], np.cumsum(paired_use)])
    # Count the paired positions.
    n_paired = paired_cumsum[-1]
    # Get the false positive rate: (paired and reactive) / paired.
    fpr = 1. - paired_cumsum / n_paired
    # Get the true positive rate: (unpaired and reactive) / unpaired.
    tpr = 1. - (np.arange(n_use + 1) - paired_cumsum) / (n_use - n_paired)
    # Traditionally, false positive rate is plotted on the x-axis and
    # true positive rate on the y-axis of an ROC curve.
    return fpr, tpr


def compute_auc_roc(fpr: np.ndarray, tpr: np.ndarray):
    """ Compute the area under the curve (AUC) of the receiver operating
    characteristic (ROC). """
    return -np.vdot(np.diff(fpr), tpr[1:])
