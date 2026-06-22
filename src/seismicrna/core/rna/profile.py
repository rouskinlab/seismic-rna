from __future__ import annotations
from functools import cached_property
from pathlib import Path


from .base import RNARegion
from .. import path
from ..logs import format_sample_reference_region

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


class RNAProfile(RNARegion):
    """Mutational profile of an RNA."""

    def __init__(
        self,
        *,
        sample: str,
        branches: dict[str, str],
        mus_reg: str,
        mus_name: str,
        mus: pd.Series,
        **kwargs,
    ):
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
        mus: pandas.Series
            Data for the mutational profile (i.e. mutation rates).
        """
        import pandas as pd

        super().__init__(**kwargs)
        self.sample = sample
        self.branches = branches
        self.mus_reg = mus_reg
        self.mus_name = mus_name
        if not isinstance(mus, pd.Series):
            raise TypeError(
                f"Expected data to be a Series, but got {type(mus).__name__}"
            )
        if mus.min() < 0.0 or mus.max() > 1.0:
            raise ValueError(f"Got mutation rates outside [0, 1]:\n{mus}")
        self.mus = mus.reindex(self.region.range)

    def __str__(self):
        srr = format_sample_reference_region(self.sample, self.ref, self.reg)
        return f"{type(self).__name__} of {srr}"

    @cached_property
    def init_args(self):
        return super().init_args | dict(
            sample=self.sample,
            branches=self.branches,
            mus_reg=self.mus_reg,
            mus_name=self.mus_name,
            mus=self.mus,
        )

    def _renumber_from_args(self, seq5: int):
        import pandas as pd

        return super()._renumber_from_args(seq5) | dict(
            mus=pd.Series(self.mus.values, index=self.region.renumber_from(seq5).range)
        )

    @property
    def profile(self):
        """Name of the mutational profile."""
        return f"{self.mus_reg}__{self.mus_name}"

    def _get_dir_fields(self, top: Path, branch: str):
        """Get the path fields for the directory of this RNA.

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
        return {
            path.TOP: top,
            path.SAMPLE: self.sample,
            path.STEP: path.FOLD_STEP,
            path.BRANCHES: path.add_branch(path.FOLD_STEP, branch, self.branches),
            path.REF: self.ref,
            path.REG: self.reg,
        }

    def _get_dir(self, top: Path, branch: str):
        """Get the directory in which to write files of this RNA."""
        return path.builddir(path.REG_DIR_SEGS, self._get_dir_fields(top, branch))

    def _get_file(
        self, top: Path, branch: str, path_seg: path.PathSegment, path_fields
    ):
        """Get the path to a file of the RNA."""
        return self._get_dir(top, branch).joinpath(path_seg.build(path_fields))

    def get_ct_file(self, top: Path, branch: str):
        """Get the path to the connectivity table (CT) file."""
        return self._get_file(
            top,
            branch,
            path.ConnectTableSeg,
            {path.PROFILE: self.profile, path.EXT: path.CT_EXT},
        )

    def get_db_file(self, top: Path, branch: str):
        """Get the path to the dot-bracket (DB) file."""
        return self._get_file(
            top,
            branch,
            path.DotBracketSeg,
            {path.PROFILE: self.profile, path.EXT: path.DOT_EXTS[0]},
        )
