from functools import cached_property

from .profile import RNAProfile
from .roc import compute_auc, compute_roc_curve, compute_rolling_auc
from .struct import RNAStructure
from .. import path
from ..error import InconsistentValueError


class RNAState(RNAStructure, RNAProfile):
    """ RNA secondary structure with mutation rates. """

    @classmethod
    def from_struct_profile(cls, struct: RNAStructure, profile: RNAProfile):
        """ Make an RNAState from an RNAStructure and an RNAProfile. """
        if struct.region.ref != profile.region.ref:
            raise InconsistentValueError(
                "Reference names differ between "
                f"structure ({repr(struct.region.ref)}) "
                f"and profile ({repr(profile.region.ref)})"
            )
        return cls(region=struct.region,
                   title=struct.title,
                   pairs=struct.pairs,
                   sample=profile.sample,
                   branches=path.add_branch(path.FOLD_STEP,
                                            struct.branch,
                                            profile.branches),
                   data_reg=profile.data_reg,
                   data_name=profile.data_name,
                   data=profile.data)

    @cached_property
    def roc(self):
        return compute_roc_curve(self.is_paired, self.data)

    @cached_property
    def auc(self):
        return compute_auc(*self.roc)

    def rolling_auc(self, size: int, min_data: int = 2):
        return compute_rolling_auc(self.is_paired, self.data, size, min_data)
