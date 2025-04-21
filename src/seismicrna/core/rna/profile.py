from functools import cached_property

import numpy as np
import pandas as pd

from .base import RNARegion
from ..validate import require_isinstance, require_between


class RNAProfile(RNARegion):
    """ Mutational profile of an RNA. """

    def __init__(self, *,
                 sample: str,
                 branches: dict[str, str],
                 mus_reg: str,
                 mus_name: str,
                 mus: pd.Series,
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
        require_isinstance("mus", mus, pd.Series)
        if np.count_nonzero(np.isnan(mus)) < mus.size:
            require_between("mus.min()", mus.min(), 0., 1.)
            require_between("mus.max()", mus.max(), 0., 1.)
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
