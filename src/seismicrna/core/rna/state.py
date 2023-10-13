from functools import cached_property

from .profile import RnaProfile
from .roc import compute_auc_roc, compute_roc_curve
from .struct import Rna2dStructure


class RnaState(Rna2dStructure, RnaProfile):
    """ RNA secondary structure with mutation rates. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @cached_property
    def roc(self):
        return compute_roc_curve(self.table != 0, self.reacts)

    @cached_property
    def auc(self):
        return compute_auc_roc(*self.roc)
