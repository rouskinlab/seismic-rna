from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from ..core import path
from ..core.mu import calc_pseudoenergies
from ..core.rna import RNAProfile
from ..core.seq import write_fasta
from ..core.validate import require_atleast, require_between


class RNAFoldProfile(RNAProfile):

    @classmethod
    def from_profile(cls, profile: RNAProfile, **kwargs):
        """ Make an RNAFoldProfile from an RNAProfile. """
        return cls(**profile.init_args, **kwargs)

    def __init__(self, *,
                 fold_temp: float | int,
                 fold_fpaired: float | int,
                 mu_eps: float | int,
                 **kwargs):
        """
        Parameters
        ----------
        fold_temp: float | int
            Predict structures at this temperature (Kelvin).
        fold_fpaired: float | int
            Assume this is the fraction of paired bases.
        mu_eps: float | int
            Clip mutation rates to the interval [mu_eps, 1 - mu_eps]
            during folding and structure prediction.
        """
        super().__init__(**kwargs)
        require_atleast(
            "fold_temp", fold_temp, 0., classes=(float, int)
        )
        self.fold_temp = fold_temp
        require_between(
            "fold_fpaired", fold_fpaired, 0., 1., classes=(float, int)
        )
        self.fold_fpaired = fold_fpaired
        require_between(
            "mu_eps", mu_eps, 0.0, 0.5, classes=(float, int)
        )
        self.mu_eps = mu_eps

    @cached_property
    def mus_clipped(self):
        """ Mutation rates after clipping to [mu_eps, 1 - mu_eps]. """
        return pd.Series(np.clip(self.mus, self.mu_eps, 1. - self.mu_eps),
                         self.mus.index)

    @cached_property
    def pseudoenergies(self):
        """ Pseudoenergies (kcal/mol) for structure prediction. """
        return calc_pseudoenergies(
            self.mus_clipped.dropna(),
            self.fold_temp,
            self.fold_fpaired
        ).reindex(self.mus_clipped.index)

    @cached_property
    def intercept(self):
        """ Intercept parameter (kcal/mol) for structure prediction. """
        pseudoenergies = self.pseudoenergies[np.isfinite(self.pseudoenergies)]
        if pseudoenergies.size == 0:
            return 0.
        return float(pseudoenergies.min())

    @cached_property
    def slope(self):
        """ Slope parameter (kcal/mol) for structure prediction. """
        pseudoenergies = self.pseudoenergies[np.isfinite(self.pseudoenergies)]
        if pseudoenergies.size == 0:
            return 0.
        return float((pseudoenergies.max() - self.intercept) / np.log(2.))

    @cached_property
    def pseudomus(self):
        """ Pseudo-mutation rates for structure prediction. """
        if self.slope == 0.:
            # This happens if all pseudoenergies equal the intercept.
            return pd.Series(0., self.pseudoenergies.index)
        # During folding, the pseudoenergies will be calculated as:
        # pseudoenergies = slope * log(pseudomus + 1) + intercept
        # Therefore:
        # pseudomus = exp((pseudoenergies - intercept) / slope) - 1
        return np.expm1((self.pseudoenergies - self.intercept) / self.slope)

    def _get_dir_fields(self, top: Path, branch: str):
        """ Get the path fields for the directory of this RNA.

        Parameters
        ----------
        top: pathlib.Path
            Top-level directory.
        branch: str
            Branch to add (optional, for folding).

        Returns
        -------
        dict[str, str | pathlib.Path]
            Path fields.
        """
        return {path.TOP: top,
                path.SAMPLE: self.sample,
                path.STEP: path.FOLD_STEP,
                path.BRANCHES: path.add_branch(path.FOLD_STEP,
                                               branch,
                                               self.branches),
                path.REF: self.ref,
                path.REG: self.reg}

    def _get_dir(self, top: Path, branch: str):
        """ Get the directory in which to write files of this RNA. """
        return path.builddir(path.REG_DIR_SEGS,
                             self._get_dir_fields(top, branch))

    def _get_file(self,
                  top: Path,
                  branch: str,
                  path_seg: path.PathSegment,
                  path_fields):
        """ Get the path to a file of the RNA. """
        return self._get_dir(top, branch).joinpath(path_seg.build(path_fields))

    def get_fasta(self, top: Path, branch: str):
        """ Get the path to the FASTA file. """
        return self._get_file(top,
                              branch,
                              path.FastaSeg,
                              {path.REF: self.profile,
                               path.EXT: path.FASTA_EXTS[0]})

    def get_ct_file(self, top: Path, branch: str):
        """ Get the path to the connectivity table (CT) file. """
        return self._get_file(top,
                              branch,
                              path.ConnectTableSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.CT_EXT})

    def get_db_file(self, top: Path, branch: str):
        """ Get the path to the dot-bracket (DB) file. """
        return self._get_file(top,
                              branch,
                              path.DotBracketSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.DOT_EXTS[0]})

    def get_pseudodata_file(self, top: Path, branch: str):
        """ Get the path to the pseudo-mutation rates file. """
        return self._get_file(top,
                              branch,
                              path.PseudoMusSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.PSEUDOMUS_EXT})

    def get_vienna_file(self, top: Path, branch: str):
        """ Get the path to the vienna file. """
        return self._get_file(top,
                              branch,
                              path.ViennaSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.VIENNA_EXT})

    def get_command_file(self, top: Path, branch: str):
        """ Get the path to the vienna command file. """
        return self._get_file(top,
                              branch,
                              path.CommandSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.COMMAND_EXT})

    def get_varna_color_file(self, top: Path, branch: str):
        """ Get the path to the VARNA color file. """
        return self._get_file(top,
                              branch,
                              path.VarnaColorSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.TXT_EXT})

    def write_fasta(self, top: Path, branch: str):
        """ Write the RNA sequence to a FASTA file. """
        fasta = self.get_fasta(top, branch)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def to_pseudomus(self, top: Path, branch: str):
        """ Write the pseudo-mutation rates to a file. """
        # The pseudo-mutation rates must be numbered starting from 1 at
        # the beginning of the region, even if the region does not start
        # at 1. Renumber the region from 1.
        pseudomus = self.pseudomus.copy()
        pseudomus.index = self.region.range_one
        # Drop bases with missing data to make RNAstructure ignore them.
        pseudomus.dropna(inplace=True)
        # Write the pseudo-mutation rates to the pseudo-mutation rates file.
        pseudodata_file = self.get_pseudodata_file(top, branch)
        pseudomus.to_csv(pseudodata_file, sep="\t", header=False)
        return pseudodata_file

    def to_varna_color_file(self, top: Path, branch: str):
        """ Write the VARNA colors to a file. """
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.mus.fillna(-1.)
        # Write the values to the VARNA color file.
        varna_color_file = self.get_varna_color_file(top, branch)
        varna_color.to_csv(varna_color_file,
                           float_format="%f",
                           header=False,
                           index=False)
        return varna_color_file
