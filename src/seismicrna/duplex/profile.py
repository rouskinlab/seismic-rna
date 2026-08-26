from __future__ import annotations
from functools import cached_property
from pathlib import Path

from ..core.arg.cli import FOLD_ENERGY_METHOD_CORDERO, FOLD_ENERGY_METHOD_DEIGAN
from ..core.error import IncompatibleOptionsError
from ..core.mu.dist import calc_pseudoenergies
from ..core.mu.scale import winsorize
from ..core.rna.profile import RNAProfile
from ..core.seq.region import BASE_NAME
from ..core.seq.xna import BASEA, BASEC
from ..fold.profile import RNAFoldProfile
from .dataset import LINKER_LENGTH

# Character separating the two strands of a duplex in the input to
# RNAcofold and in its output.
STRAND_SEP = "&"


class RNACofoldProfile(RNAFoldProfile):
    """Mutational profile of a duplex of two RNAs to cofold together.

    The two strands are joined (5' strand first) by a short unpaired
    linker into a single sequence numbered from 1, so that the duplex is
    represented exactly like any other folded region and its outputs are
    compatible with the rest of SEISMIC-RNA (e.g. drawing and graphing).
    The duplex is folded by RNAcofold on the two strands alone; the
    linker only exists in the output, where it turns each inter-strand
    base pair into a drawable loop.
    """

    @classmethod
    def from_duplex(cls, profile: RNAProfile, cut: int, **kwargs):
        """Make a duplex profile from an already-fused (duplex) profile and
        its strand-break position."""
        return cls(
            region=profile.region,
            sample=profile.sample,
            branches=profile.branches,
            mus_reg=profile.mus_reg,
            mus_name=profile.mus_name,
            mus=profile.mus,
            cut=cut,
            **kwargs,
        )

    def __init__(self, *, cut: int, fold_fpaired: float | int = 0.5, **kwargs):
        """
        Parameters
        ----------
        cut: int
            Number of positions in the 5' strand, i.e. the position of
            the strand break (before the linker) within the fused
            sequence.
        fold_fpaired: float | int
            Assumed fraction of paired bases, used when computing Cordero
            pseudoenergies.
        """
        super().__init__(**kwargs)
        self.cut = cut
        self.fold_fpaired = fold_fpaired

    @cached_property
    def mus_normalized(self):
        """Mutation rates after normalizing and winsorizing each strand
        independently, since the two strands can come from different
        samples with different overall reactivities."""
        import pandas as pd

        strand5 = winsorize(self.mus.iloc[: self.cut], self.fold_quantile)
        # The linker (all missing) stays with the 3' strand and is left
        # untouched by winsorize.
        strand3 = winsorize(self.mus.iloc[self.cut :], self.fold_quantile)
        return pd.concat([strand5, strand3])

    def _strand_pseudoenergies(self, strand_mus):
        """Cordero pseudoenergies (kcal/mol) for one strand's A and C
        reactivities (DMS informs A/C)."""
        import numpy as np

        mus = strand_mus.dropna()
        bases = mus.index.get_level_values(BASE_NAME)
        mus_ac = mus.loc[np.isin(bases, [BASEA, BASEC])]
        if mus_ac.empty:
            # A data-less strand (e.g. a bare partner sequence) has no
            # reactivities to fit.
            return mus_ac
        return calc_pseudoenergies(mus_ac, self.fold_temp_k, self.fold_fpaired)

    @cached_property
    def pseudoenergies(self):
        """Cordero pseudoenergies (kcal/mol) for structure prediction.
        Each strand is fit independently (its own scale factor), because
        the two strands can come from different samples whose overall
        reactivities differ; pooling them would bias the fit. The
        resulting energies share one kcal/mol scale, so a single
        pseudomus/shapeMethod still reproduces them for RNAcofold."""
        import pandas as pd

        strand5 = self.mus.iloc[: self.cut]
        # The linker (all missing) rides with the 3' strand and drops out.
        strand3 = self.mus.iloc[self.cut :]
        return pd.concat(
            [self._strand_pseudoenergies(strand5), self._strand_pseudoenergies(strand3)]
        ).reindex(self.mus.index)

    @cached_property
    def intercept_param(self):
        """Intercept parameter (kcal/mol) for structure prediction."""
        if self.pseudoenergies.size == 0:
            return 0.0
        return float(self.pseudoenergies.min())

    @cached_property
    def slope_param(self):
        """Slope parameter (kcal/mol) for structure prediction."""
        if self.pseudoenergies.size == 0:
            return 1.0
        return float(self.pseudoenergies.max()) - self.intercept_param

    @cached_property
    def pseudomus(self):
        """Pseudo-mutation rates for structure prediction."""
        import numpy as np

        if self.slope_param == 0.0:
            # This happens if all pseudoenergies equal the intercept.
            return self.pseudoenergies - self.intercept_param
        # RNAcofold applies the Deigan pseudoenergy
        #   pseudoenergies = slope * ln(pseudomus + 1) + intercept
        # so invert it to reproduce the target pseudoenergies exactly:
        #   pseudomus = exp((pseudoenergies - intercept) / slope) - 1
        return np.expm1((self.pseudoenergies - self.intercept_param) / self.slope_param)

    def _unsupported_method_error(self):
        return IncompatibleOptionsError(
            f"RNAcofold cannot fold a duplex with the {self.fold_energy_method} "
            f"energy method; use {FOLD_ENERGY_METHOD_DEIGAN} (SHAPE reactivities) "
            f"or {FOLD_ENERGY_METHOD_CORDERO} (DMS pseudo-energies)"
        )

    @property
    def _fold_data(self):
        """Per-position data to feed RNAcofold (over the fused region):
        the winsorized reactivities for Deigan (applied directly), or the
        Cordero pseudo-mutation rates for DMS."""
        if self.fold_energy_method == FOLD_ENERGY_METHOD_DEIGAN:
            return self.mus_normalized
        if self.fold_energy_method == FOLD_ENERGY_METHOD_CORDERO:
            return self.pseudomus
        raise self._unsupported_method_error()

    @cached_property
    def cofold_shape_method(self):
        """SHAPE incorporation method for RNAcofold's ``--shapeMethod``
        (RNAcofold has no modern ``--sp-strategy`` interface, so only the
        Deigan ``Dm..b..`` form is available). Deigan passes the user's
        slope/intercept straight through; Cordero passes the slope and
        intercept that back-transform the DMS pseudo-energies."""
        if self.fold_energy_method == FOLD_ENERGY_METHOD_DEIGAN:
            return f"Dm{round(self.deigan_slope, 3)}b{round(self.deigan_intercept, 3)}"
        if self.fold_energy_method == FOLD_ENERGY_METHOD_CORDERO:
            # RNAcofold applies slope * ln(pseudomus + 1) + intercept.
            return f"Dm{round(self.slope_param, 3)}b{round(self.intercept_param, 3)}"
        raise self._unsupported_method_error()

    def write_fasta(self, top: Path, branch: str):
        """Write the two strands to a FASTA file, separated by ``&`` and
        without the linker, as the input for RNAcofold."""
        seq = str(self.seq)
        strand5 = seq[: self.cut]
        strand3 = seq[self.cut + LINKER_LENGTH :]
        fasta = self.get_fasta(top, branch)
        with open(fasta, "w") as f:
            f.write(f">{self.region.ref_reg}\n{strand5}{STRAND_SEP}{strand3}\n")
        return fasta

    def write_mus_file(self, top: Path, branch: str):
        """Write the folding data (reactivities for Deigan, or Cordero
        pseudo-mutation rates) for RNAcofold, numbered along the two
        strands alone (i.e. without the linker positions)."""
        import numpy as np
        import pandas as pd

        # Drop the linker from the folding data and renumber the two
        # strands continuously starting from 1, matching the sequence
        # given to RNAcofold.
        values = self._fold_data.values
        strands = np.concatenate(
            [values[: self.cut], values[self.cut + LINKER_LENGTH :]]
        )
        mus = pd.Series(strands, index=range(1, strands.size + 1))
        mus.dropna(inplace=True)
        mus_file = self.get_mus_file(top, branch)
        mus.to_csv(mus_file, sep="\t", header=False)
        return mus_file
