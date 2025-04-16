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
                   mus_reg=profile.mus_reg,
                   mus_name=profile.mus_name,
                   mus=profile.mus,
                   fold_temp=profile.fold_temp,
                   fold_fpaired=(struct.is_paired.mean()
                                 if struct.is_paired.size > 0
                                 else profile.fold_fpaired))

    def _get_structs_args(self, terminal_pairs: bool):
        if terminal_pairs:
            return self.is_paired,
        return self.is_paired_internally, self.is_unpaired

    def calc_roc(self, terminal_pairs: bool = True):
        """ Calculate the receiver operating characteristic (ROC) curve.

        Parameters
        ----------
        terminal_pairs: bool
            Whether to count terminal base pairs as paired (if True)
            or as neither paired nor unpaired (if False).
        """
        return compute_roc_curve(self.data,
                                 *self._get_structs_args(terminal_pairs))

    def calc_auc(self, terminal_pairs: bool = True):
        """ Calculate the area under the ROC curve (AUC-ROC).

        Parameters
        ----------
        terminal_pairs: bool
            Whether to count terminal base pairs as paired (if True)
            or as neither paired nor unpaired (if False).
        """
        return compute_auc(*self.calc_roc(terminal_pairs))

    def calc_auc_rolling(self,
                         size: int,
                         min_data: int = 2,
                         terminal_pairs: bool = True):
        """ Calculate the area under the ROC curve (AUC-ROC).
        
        Parameters
        ----------
        size: int
            Size of the window.
        min_data: int
            Minimum number of data in a window to use it, otherwise NaN.
        terminal_pairs: bool
            Whether to count terminal base pairs as paired (if True)
            or as neither paired nor unpaired (if False).
        """
        return compute_rolling_auc(self.data,
                                   *self._get_structs_args(terminal_pairs),
                                   size=size,
                                   min_data=min_data)
