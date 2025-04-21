from __future__ import annotations

from abc import ABC

from ..core import path
from ..core.io import RegFileIO
from ..core.report import (RegReport,
                           ProfileF,
                           FoldViennaF,
                           FoldTempF,
                           FoldFPairedF,
                           FoldMuEpsF,
                           FoldMaxDistF,
                           FoldMinFreeEnergyF,
                           FoldMaxStructsF,
                           FoldPercent,
                           CommandsChecksumF,
                           ConstraintChecksumF)


class FoldIO(RegFileIO, ABC):

    @classmethod
    def get_step(cls):
        return path.FOLD_STEP


class FoldReport(RegReport, FoldIO, ABC):

    @classmethod
    def get_file_seg_type(cls):
        return path.FoldRepSeg

    @classmethod
    def get_param_report_fields(cls):
        return [ProfileF,
                FoldTempF,
                FoldViennaF,
                FoldFPairedF,
                FoldMuEpsF,
                FoldMaxDistF,
                FoldMinFreeEnergyF,
                FoldMaxStructsF,
                FoldPercent,
                *super().get_param_report_fields()]

    @classmethod
    def get_checksum_report_fields(cls):
        return [CommandsChecksumF,
                ConstraintChecksumF]
