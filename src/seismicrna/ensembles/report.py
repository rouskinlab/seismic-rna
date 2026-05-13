from ..core import path
from ..core.report import (RegReport,
                           TileLengthF,
                           TileMinOverlapF,
                           EraseTilesF,
                           PairFdrF,
                           MinPairsF,
                           ThresholdMultiplierF,
                           MinClusterLengthF,
                           MaxClusterLengthF,
                           GapModeF)
from .io import EnsemblesIO


class EnsemblesReport(RegReport, EnsemblesIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.EnsemblesRepSeg

    @classmethod
    def get_param_report_fields(cls):
        return [
            TileLengthF,
            TileMinOverlapF,
            EraseTilesF,
            PairFdrF,
            MinPairsF,
            ThresholdMultiplierF,
            MinClusterLengthF,
            MaxClusterLengthF,
            GapModeF,
            *super().get_param_report_fields()
        ]
