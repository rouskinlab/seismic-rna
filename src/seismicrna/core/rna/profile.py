from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd

from .section import RnaSection
from .. import path
from ..seq import write_fasta


class RnaProfile(RnaSection):
    """ Mutational profile of an RNA from a specific sample. """

    def __init__(self,
                 *args,
                 sample: str,
                 data_sect: str,
                 reacts: pd.Series,
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
        return path.builddir(*path.FOLD_SECT_DIR_SEGS,
                             top=out_dir,
                             cmd=path.CMD_FLD_DiR,
                             sample=self.sample,
                             ref=self.ref,
                             sect=self.data_sect,
                             fold_sect=self.sect)

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
