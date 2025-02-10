from __future__ import annotations

from abc import ABC

from ..core import path
from ..core.io import RegFileIO
from ..core.report import (RegReport,
                           ProfileF,
                           Quantile,
                           FoldTempF,
                           FoldMaxDistF,
                           FoldMinFreeEnergyF,
                           FoldMaxStructsF,
                           FoldPercent)


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
                Quantile,
                FoldTempF,
                FoldMaxDistF,
                FoldMinFreeEnergyF,
                FoldMaxStructsF,
                FoldPercent,
                *super().get_param_report_fields()]
