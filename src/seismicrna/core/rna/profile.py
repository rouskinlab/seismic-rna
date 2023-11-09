from pathlib import Path

import numpy as np
import pandas as pd

from .section import RnaSection
from .. import path
from ..seq import write_fasta


class RnaProfile(RnaSection):
    """ Mutational profile of an RNA from a specific sample. """

    @classmethod
    def dir_segs(cls):
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.SectSeg,
                path.FoldSectSeg)

    def __init__(self, *,
                 sample: str,
                 data_sect: str,
                 data: pd.Series,
                 **kwargs):
        super().__init__(**kwargs)
        self.sample = sample
        self.data_sect = data_sect
        if np.any(data < 0.) or np.any(data > 1.):
            raise ValueError(f"Got reactivities outside [0, 1]: {data}")
        self.data = data.reindex(self.section.range)

    def dir_fields(self, top: Path):
        return dict(top=top,
                    cmd=path.CMD_FOLD_DIR,
                    sample=self.sample,
                    ref=self.ref,
                    sect=self.data_sect,
                    fold_sect=self.sect)

    def get_dir(self, top: Path):
        """ Get the directory in which this RNA's files will go. """
        return path.builddir(*self.dir_segs(), **self.dir_fields(top))

    def get_file(self, top: Path, segment: path.Segment, **kwargs):
        """ Get the path to a file of the RNA sequence. """
        return self.get_dir(top).joinpath(segment.build(**kwargs))

    def get_fasta(self, top: Path):
        """ Get the path to the FASTA file of the RNA sequence. """
        return self.get_file(top,
                             path.FastaSeg,
                             ref=self.title,
                             ext=path.FASTA_EXTS[0])

    def ct_file(self, top: Path):
        """ Get the path to the CT file of the RNA. """
        return self.get_file(top,
                             path.ConnectTableSeg,
                             struct=self.title,
                             ext=path.CT_EXT)

    def dot_file(self, top: Path):
        """ Get the path to the DOT file of the RNA. """
        return self.get_file(top,
                             path.DotBracketSeg,
                             struct=self.title,
                             ext=path.DOT_EXTS[0])

    def dms_file(self, top: Path):
        """ Get the path to the DMS data file of the RNA. """
        return self.get_file(top,
                             path.DmsReactsSeg,
                             reacts=self.title,
                             ext=path.DMS_EXT)

    def varna_color_file(self, top: Path):
        """ Get the path to the VARNA color file of the RNA. """
        return self.get_file(top,
                             path.VarnaColorSeg,
                             reacts=self.title,
                             ext=path.TXT_EXT)

    def to_fasta(self, top: Path):
        """ Write the RNA sequence to a FASTA file. """
        fasta = self.get_fasta(top)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_dms(self, top: Path):
        """ Write the DMS reactivities to a DMS file. """
        # The DMS reactivities must be numbered starting from 1 at the
        # beginning of the section, even if the section does not start
        # at 1. Renumber the section from 1.
        dms = self.data.copy()
        dms.index = self.section.range_one
        # Drop bases with missing data to make RNAstructure ignore them.
        dms.dropna(inplace=True)
        # Write the DMS reactivities to the DMS file.
        dms_file = self.dms_file(top)
        dms.to_csv(dms_file, sep="\t", header=False)
        return dms_file

    def to_varna_color_file(self, top: Path):
        """ Write the VARNA colors to a file. """
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.data.fillna(-1.)
        # Write the values to the VARNA color file.
        varna_color_file = self.varna_color_file(top)
        varna_color.to_csv(varna_color_file, header=False, index=False)
        return varna_color_file

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
