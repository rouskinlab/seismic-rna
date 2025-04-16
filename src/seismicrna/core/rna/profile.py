from functools import cached_property
from pathlib import Path

import pandas as pd

from .base import RNARegion
from .. import path
from ..arg import opt_fold_fpaired, opt_fold_temp
from ..mu import calc_pseudoenergies
from ..seq import write_fasta
from ..validate import require_isinstance, require_between


class RNAProfile(RNARegion):
    """ Mutational profile of an RNA. """

    def __init__(self, *,
                 sample: str,
                 branches: dict[str, str],
                 mus_reg: str,
                 mus_name: str,
                 mus: pd.Series,
                 fold_temp: float | int = opt_fold_temp.default,
                 fold_fpaired: float | int = opt_fold_fpaired.default,
                 **kwargs):
        """
        Parameters
        ----------
        sample: str
            Name of the sample from which the mutational profile comes.
        branches: dict[str, str]
            Branches of the workflow.
        mus_reg: str
            Name of the region from which the mutational profile comes.
        mus_name: str
            Name of the mutational profile (e.g. "cluster_2-1").
        mus: pd.Series
            Data for the mutational profile (i.e. mutation rates).
        fold_temp: float | int
            Predict structures at this temperature (Kelvin).
        fold_fpaired: float | int
            Assume this is the fraction of paired bases.
        """
        super().__init__(**kwargs)
        self.sample = sample
        self.branches = branches
        self.mus_reg = mus_reg
        self.mus_name = mus_name
        self.fold_temp = fold_temp
        self.fold_fpaired = fold_fpaired
        require_isinstance("data", mus, pd.Series)
        if mus.size > 0:
            require_between("data.min()", mus.min(), 0., 1.)
            require_between("data.max()", mus.max(), 0., 1.)
        self.mus = mus.reindex(self.region.range)

    @cached_property
    def init_args(self):
        return super().init_args | dict(sample=self.sample,
                                        branches=self.branches,
                                        mus_reg=self.mus_reg,
                                        mus_name=self.mus_name,
                                        mus=self.mus)

    def _renumber_from_args(self, seq5: int):
        return super()._renumber_from_args(seq5) | dict(
            data=pd.Series(self.mus.values,
                           index=self.region.renumber_from(seq5).range)
        )

    @property
    def profile(self):
        """ Name of the mutational profile. """
        return f"{self.mus_reg}__{self.mus_name}"
    
    @cached_property
    def pseudoenergies(self):
        """ Pseudoenergies (kcal/mol) for structure prediction. """
        return calc_pseudoenergies(self.mus.dropna(),
                                   self.fold_temp,
                                   self.fold_fpaired).reindex(self.mus.index)
    
    @cached_property
    def intercept_param(self):
        """ Intercept parameter (kcal/mol) for structure prediction. """
        if self.pseudoenergies.size == 0:
            return 0.
        return float(self.pseudoenergies.min())
    
    @cached_property
    def slope_param(self):
        """ Slope parameter (kcal/mol) for structure prediction. """
        if self.pseudoenergies.size == 0:
            return 1.
        return float(self.pseudoenergies.max()) - self.intercept_param

    @cached_property
    def pseudomus(self):
        """ Pseudo-mutation rates for structure prediction. """
        if self.slope_param == 0.:
            # This happens if all pseudoenergies equal the intercept.
            return self.pseudoenergies - self.intercept_param
        # During folding, the pseudoenergies will be calculated as:
        # pseudoenergies = slope * pseudodata + intercept
        # Thus: pseudodata = (pseudoenergies - intercept) / slope
        return (self.pseudoenergies - self.intercept_param) / self.slope_param

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

    def get_varna_color_file(self, top: Path, branch: str):
        """ Get the path to the VARNA color file. """
        return self._get_file(top,
                              branch,
                              path.VarnaColorSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.TXT_EXT})

    def to_fasta(self, top: Path, branch: str):
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
