from functools import cached_property
from pathlib import Path

import pandas as pd

from ..core import path
from ..core.arg import (
    FOLD_ENERGY_METHODS,
    FOLD_ENERGY_METHOD_CORDERO,
    FOLD_ENERGY_METHOD_DEIGAN,
    FOLD_ENERGY_METHOD_EDDY,
)
from ..core.logs import logger
from ..core.mu import winsorize
from ..core.rna import RNAProfile
from ..core.seq import write_fasta
from ..core.validate import (
    require_atleast,
    require_between,
    require_isin,
    require_isinstance,
)

ZERO_CELSIUS = 273.15  # Kelvin


def celsius_to_kelvin(temp_c: float | int):
    """Convert a temperature from Celsius to Kelvin."""
    require_atleast(
        "temp_c", temp_c, -ZERO_CELSIUS, f"-{ZERO_CELSIUS} °C", classes=(float, int)
    )
    return temp_c + ZERO_CELSIUS


def kelvin_to_celsius(temp_k: float | int):
    """Convert a temperature from Kelvin to Celsius."""
    require_atleast("temp_k", temp_k, 0.0, "0 K", classes=(float, int))
    return temp_k - ZERO_CELSIUS


def guess_temperature_to_celsius(temp: float | int):
    """Guess whether a temperature is in Celsius or Kelvin and return
    as Celsius."""
    require_isinstance("temp", temp, (float, int))
    if temp < -ZERO_CELSIUS:
        raise ValueError(f"Temperature is too low to be in Celsius or Kelvin: {temp}")
    if temp < (ZERO_CELSIUS + 100.0) / 2.0:
        # Assume Celsius if the temperature is less than the midpoint of
        # 100 (boiling point of water in Celsius) and 273.15 (freezing
        # point of water in Kelvin).
        return temp
    # Assume Kelvin if the temperature is greater than the midpoint.
    logger.warning(
        "Assuming temperature {} is in Kelvin: if you meant it to be in "
        "Celsius, please enter it as {}",
        temp,
        round(celsius_to_kelvin(temp), 2),
    )
    return kelvin_to_celsius(temp)


class RNAFoldProfile(RNAProfile):
    @classmethod
    def from_profile(cls, profile: RNAProfile, **kwargs):
        """Make an RNAFoldProfile from an RNAProfile."""
        return cls(**profile.init_args, **kwargs)

    def __init__(
        self,
        *,
        fold_temp: float | int,
        fold_energy_method: str,
        fold_quantile: float | int,
        deigan_slope: float | int,
        deigan_intercept: float | int,
        **kwargs,
    ):
        """
        Parameters
        ----------
        fold_temp: float | int
            Predict structures at this temperature (Kelvin).
        fold_energy_method: str
            Method for normalizing the reactivities for folding.
        fold_quantile: float
            Normalize mutation rates to this quantile, then winsorize.
        """
        super().__init__(**kwargs)
        # Folding temperature (Celsius)
        self.fold_temp_c = guess_temperature_to_celsius(fold_temp)
        require_isin(
            "fold_energy_method",
            fold_energy_method,
            FOLD_ENERGY_METHODS,
            "FOLD_ENERGY_METHODS",
        )
        self.fold_energy_method = fold_energy_method
        require_between("fold_quantile", fold_quantile, 0.0, 1.0, classes=(float, int))
        self.fold_quantile = float(fold_quantile)
        require_isinstance("deigan_slope", deigan_slope, (float, int))
        self.deigan_slope = float(deigan_slope)
        require_isinstance("deigan_intercept", deigan_intercept, (float, int))
        self.deigan_intercept = float(deigan_intercept)

    @property
    def fold_temp_k(self):
        """Folding temperature (Kelvin)."""
        return celsius_to_kelvin(self.fold_temp_c)

    def get_rnastructure_shape_args(self, top: Path, branch: str):
        """Get the SHAPE/DMS arguments for Fold/ShapeKnots."""
        mus_file = self.get_mus_file(top, branch)
        if self.fold_energy_method == FOLD_ENERGY_METHOD_DEIGAN:
            # For Deigan, use the SHAPE file, even if data aren't SHAPE.
            return dict(
                dms_file=None,
                shape_file=mus_file,
                deigan_slope=self.deigan_slope,
                deigan_intercept=self.deigan_intercept,
            )
        if self.fold_energy_method == FOLD_ENERGY_METHOD_CORDERO:
            # For Cordero, use the DMS file with no slope/intercept.
            return dict(
                dms_file=mus_file,
                shape_file=None,
                deigan_slope=None,
                deigan_intercept=None,
            )
        raise ValueError(
            "Invalid fold_energy_method for Fold/ShapeKnots: "
            f"{repr(self.fold_energy_method)}"
        )

    @cached_property
    def rnafold_sp_strategy(self):
        """--sp-strategy string for RNAFold."""
        if self.fold_energy_method == FOLD_ENERGY_METHOD_DEIGAN:
            # For Deigan (D), specify slope and intercept.
            return f"Dm{round(self.deigan_slope, 3)}b{round(self.deigan_intercept, 3)}"
        if self.fold_energy_method == FOLD_ENERGY_METHOD_EDDY:
            # For Eddy (E), specify temperature in Celsius.
            return f"Et{round(self.fold_temp_c, 3)}"
        raise ValueError(
            f"Invalid fold_energy_method for RNAFold: {repr(self.fold_energy_method)}"
        )

    @cached_property
    def mus_normalized(self):
        """Mutation rates after normalizing and winsorizing."""
        return winsorize(self.mus, self.fold_quantile)

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

    def get_fasta(self, top: Path, branch: str):
        """Get the path to the FASTA file."""
        return self._get_file(
            top,
            branch,
            path.FastaSeg,
            {path.REF: self.profile, path.EXT: path.FASTA_EXTS[0]},
        )

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

    def get_mus_file(self, top: Path, branch: str):
        """Get the path to the mutation rates file."""
        return self._get_file(
            top,
            branch,
            path.MusSeg,
            {path.PROFILE: self.profile, path.EXT: path.MUS_EXT},
        )

    def get_vienna_file(self, top: Path, branch: str):
        """Get the path to the vienna file."""
        return self._get_file(
            top,
            branch,
            path.ViennaSeg,
            {path.PROFILE: self.profile, path.EXT: path.VIENNA_EXT},
        )

    def get_varna_color_file(self, top: Path, branch: str):
        """Get the path to the VARNA color file."""
        return self._get_file(
            top,
            branch,
            path.VarnaColorSeg,
            {path.PROFILE: self.profile, path.EXT: path.TXT_EXT},
        )

    def write_fasta(self, top: Path, branch: str):
        """Write the RNA sequence to a FASTA file."""
        fasta = self.get_fasta(top, branch)
        write_fasta(fasta, [self.seq_record])
        return fasta

    def write_mus_file(self, top: Path, branch: str):
        """Write the mutation rates to a file."""
        # The mutation rates must be numbered starting from 1 at the
        # beginning of the region, even if the region does not start
        # at 1. Renumber the region from 1.
        mus = pd.Series(self.mus_normalized.values, self.region.range_one)
        # Drop bases with missing data to make RNAstructure ignore them.
        mus.dropna(inplace=True)
        # Write the pseudo-mutation rates to the pseudo-mutation rates file.
        mus_file = self.get_mus_file(top, branch)
        mus.to_csv(mus_file, sep="\t", header=False)
        return mus_file

    def write_varna_color_file(self, top: Path, branch: str):
        """Write the VARNA colors to a file."""
        # Fill missing reactivities with -1, to signify no data.
        varna_color = self.mus_normalized.fillna(-1.0)
        # Write the values to the VARNA color file.
        varna_color_file = self.get_varna_color_file(top, branch)
        varna_color.to_csv(
            varna_color_file, float_format="%f", header=False, index=False
        )
        return varna_color_file
