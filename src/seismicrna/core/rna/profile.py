from functools import cached_property
from pathlib import Path

import pandas as pd

from .base import RNASection
from .. import path
from ..seq import write_fasta


class RNAProfile(RNASection):
    """ Mutational profile of an RNA. """

    def __init__(self, *,
                 sample: str,
                 data_sect: str,
                 data_name: str,
                 data: pd.Series,
                 **kwargs):
        """
        Parameters
        ----------
        sample: str
            Name of the sample from which the mutational profile comes.
        data_sect: str
            Name of the section from which the mutational profile comes.
        data_name: str
            Name of the mutational profile (e.g. "cluster_2-1").
        data: pandas.Series
            Data for the mutational profile (i.e. mutation rates).
        """
        super().__init__(**kwargs)
        self.sample = sample
        self.data_sect = data_sect
        self.data_name = data_name
        if not isinstance(data, pd.Series):
            raise TypeError(
                f"Expected data to be a Series, but got {type(data).__name__}"
            )
        if data.min() < 0. or data.max() > 1.:
            raise ValueError(f"Got mutation rates outside [0, 1]:\n{data}")
        self.data = data.reindex(self.section.range)

    @cached_property
    def init_args(self):
        return super().init_args | dict(sample=self.sample,
                                        data_sect=self.data_sect,
                                        data_name=self.data_name,
                                        data=self.data)

    def _renumber_from_args(self, seq5: int):
        return super()._renumber_from_args(seq5) | dict(
            data=pd.Series(self.data.values,
                           index=self.section.renumber_from(seq5).range)
        )

    @property
    def profile(self):
        """ Name of the mutational profile. """
        return f"{self.data_sect}__{self.data_name}"

    def _get_dir_fields(self, top: Path):
        """ Get the path fields for the directory of this RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        dict[str, str | pathlib.Path]
            Path fields.
        """
        return {path.TOP: top,
                path.CMD: path.CMD_FOLD_DIR,
                path.SAMP: self.sample,
                path.REF: self.ref,
                path.SECT: self.sect}

    def _get_dir(self, top: Path):
        """ Get the directory in which to write files of this RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Parent directory of files for this RNA.
        """
        return path.builddir(*path.SECT_DIR_SEGS, **self._get_dir_fields(top))

    def _get_file(self, top: Path, file_seg: path.Segment, **file_fields):
        """ Get the path to a file of the RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.
        file_seg: path.Segment
            Segment of the file component of the path.
        **file_fields
            Fields for the file segment.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_dir(top).joinpath(file_seg.build(**file_fields))

    def get_fasta(self, top: Path):
        """ Get the path to the FASTA file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.FastaSeg,
                              ref=self.profile,
                              ext=path.FASTA_EXTS[0])

    def get_ct_file(self, top: Path):
        """ Get the path to the connectivity table (CT) file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.ConnectTableSeg,
                              profile=self.profile,
                              ext=path.CT_EXT)

    def get_db_file(self, top: Path):
        """ Get the path to the dot-bracket (DB) file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.DotBracketSeg,
                              profile=self.profile,
                              ext=path.DOT_EXTS[0])

    def get_dms_file(self, top: Path):
        """ Get the path to the DMS data file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            DMS data file.
        """
        return self._get_file(top,
                              path.DmsReactsSeg,
                              profile=self.profile,
                              ext=path.DMS_EXT)

    def get_varna_color_file(self, top: Path):
        """ Get the path to the VARNA color file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            Path of the file.
        """
        return self._get_file(top,
                              path.VarnaColorSeg,
                              profile=self.profile,
                              ext=path.TXT_EXT)

    def to_fasta(self, top: Path):
        """ Write the RNA sequence to a FASTA file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            File into which the RNA sequence was written.
        """
        fasta = self.get_fasta(top)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_dms(self, top: Path):
        """ Write the DMS reactivities to a DMS file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            File into which the DMS reactivities were written.
        """
        # The DMS reactivities must be numbered starting from 1 at the
        # beginning of the section, even if the section does not start
        # at 1. Renumber the section from 1.
        dms = self.data.copy()
        dms.index = self.section.range_one
        # Drop bases with missing data to make RNAstructure ignore them.
        dms.dropna(inplace=True)
        # Write the DMS reactivities to the DMS file.
        dms_file = self.get_dms_file(top)
        dms.to_csv(dms_file, sep="\t", header=False)
        return dms_file

    def to_varna_color_file(self, top: Path):
        """ Write the VARNA colors to a file.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.

        Returns
        -------
        pathlib.Path
            File into which the VARNA colors were written.
        """
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.data.fillna(-1.)
        # Write the values to the VARNA color file.
        varna_color_file = self.get_varna_color_file(top)
        varna_color.to_csv(varna_color_file,
                           float_format="%f",
                           header=False,
                           index=False)
        return varna_color_file

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
