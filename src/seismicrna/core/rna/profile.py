from functools import cached_property
from pathlib import Path

import pandas as pd

from .base import RNARegion
from .. import path
from ..logs import logger
from ..mu import (PAIRED_ALPHA,
                  PAIRED_BETA,
                  UNPAIRED_ALPHA,
                  UNPAIRED_BETA,
                  DEFAULT_MEAN_A,
                  DEFAULT_COV_A,
                  fit_beta_mixture_model)
from ..seq import BASE_NAME, RNA, write_fasta


class RNAProfile(RNARegion):
    """ Mutational profile of an RNA. """

    def __init__(self, *,
                 sample: str,
                 branches: dict[str, str],
                 data_reg: str,
                 data_name: str,
                 data: pd.Series,
                 **kwargs):
        """
        Parameters
        ----------
        sample: str
            Name of the sample from which the mutational profile comes.
        branches: dict[str, str]
            Branches of the workflow.
        data_reg: str
            Name of the region from which the mutational profile comes.
        data_name: str
            Name of the mutational profile (e.g. "cluster_2-1").
        data: pandas.Series
            Data for the mutational profile (i.e. mutation rates).
        """
        super().__init__(**kwargs)
        self.sample = sample
        self.branches = branches
        self.data_reg = data_reg
        self.data_name = data_name
        if not isinstance(data, pd.Series):
            raise TypeError(
                f"Expected data to be a Series, but got {type(data).__name__}"
            )
        if data.min() < 0. or data.max() > 1.:
            raise ValueError(f"Got mutation rates outside [0, 1]:\n{data}")
        self.data = data.reindex(self.region.range)

    @cached_property
    def init_args(self):
        return super().init_args | dict(sample=self.sample,
                                        branches=self.branches,
                                        data_reg=self.data_reg,
                                        data_name=self.data_name,
                                        data=self.data)

    def _renumber_from_args(self, seq5: int):
        return super()._renumber_from_args(seq5) | dict(
            data=pd.Series(self.data.values,
                           index=self.region.renumber_from(seq5).range)
        )

    @property
    def profile(self):
        """ Name of the mutational profile. """
        return f"{self.data_reg}__{self.data_name}"

    @cached_property
    def beta_params(self):
        """ Beta parameters for the mutational profile. """
        beta_params = dict()
        for base in RNA.four():
            logger.detail(f"Fitting parameters for base {repr(base)}")
            is_base = self.data.index.get_level_values(BASE_NAME) == base
            beta_params[base] = fit_beta_mixture_model(self.data.loc[is_base],
                                                       DEFAULT_MEAN_A,
                                                       DEFAULT_COV_A)
        return pd.DataFrame.from_dict(beta_params, orient="index")

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

    def get_mus_file(self, top: Path, branch: str):
        """ Get the path to the mutation rate data file. """
        return self._get_file(top,
                              branch,
                              path.FoldMusSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.FOLD_MUS_EXT})

    def get_beta_file(self, top: Path, branch: str):
        """ Get the path to the beta parameter data file. """
        return self._get_file(top,
                              branch,
                              path.FoldBetaSeg,
                              {path.PROFILE: self.profile,
                               path.EXT: path.FOLD_BETA_EXT})

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

    def write_mus(self, top: Path, branch: str):
        """ Write the mutation rates to a file. """
        # The mutation rates must be numbered starting from 1 at the
        # beginning of the region, even if the region does not start
        # at 1. Renumber the region from 1.
        mus = self.data.copy()
        mus.index = self.region.range_one
        # Drop bases with missing data to make RNAstructure ignore them.
        mus.dropna(inplace=True)
        # Write the mutation rates to the file.
        mus_file = self.get_mus_file(top, branch)
        mus.to_csv(mus_file, sep="\t", header=False)
        return mus_file

    def write_beta_params(self, top: Path, branch: str):
        """ Write the beta parameters to a file. """
        lines = ["# Parameters: "
                 "A paired alpha, beta; A unpaired alpha, beta; "
                 "C paired alpha, beta; C unpaired alpha, beta; "
                 "G paired alpha, beta; G unpaired alpha, beta; "
                 "U paired alpha, beta; U unpaired alpha, beta"]
        for base, params in self.beta_params.iterrows():
            lines.extend([f"{params[PAIRED_ALPHA]} {params[PAIRED_BETA]}",
                          f"{params[UNPAIRED_ALPHA]} {params[UNPAIRED_BETA]}"])
        beta_params_file = self.get_beta_file(top, branch)
        with open(beta_params_file, "x") as f:
            f.write("\n".join(lines))
        return beta_params_file

    def write_varna_color_file(self, top: Path, branch: str):
        """ Write the VARNA colors to a file. """
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.data.fillna(-1.)
        # Write the values to the VARNA color file.
        varna_color_file = self.get_varna_color_file(top, branch)
        varna_color.to_csv(varna_color_file,
                           float_format="%f",
                           header=False,
                           index=False)
        return varna_color_file
