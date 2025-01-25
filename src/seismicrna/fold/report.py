from __future__ import annotations

from ..core import path
from ..core.report import (Report,
                           SampleF,
                           RefF,
                           RegF,
                           ProfileF,
                           Quantile,
                           FoldTempF,
                           FoldMaxDistF,
                           FoldMinFreeEnergyF,
                           FoldMaxStructsF,
                           FoldPercent)


class FoldReport(Report):

    @classmethod
    def fields(cls):
        return [SampleF,
                RefF,
                RegF,
                ProfileF,
                Quantile,
                FoldTempF,
                FoldMaxDistF,
                FoldMinFreeEnergyF,
                FoldMaxStructsF,
                FoldPercent] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.FoldRepSeg

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg, path.RefSeg, path.RegSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.FOLD_STEP}
