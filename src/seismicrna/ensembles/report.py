from ..core import path
from ..core.report import (
    RegReport,
    TileLengthF,
    TileMinOverlapF,
    EraseTilesF,
    PairFdrF,
    MinPairsF,
    ThresholdMultiplierF,
    MinClusterLengthF,
    MaxClusterLengthF,
    GapModeF,
    TileCoordsF,
    NumSignifPairsF,
    NumModulesF,
    ModuleCoordsF,
    ClusterDirsF,
    BestKsF,
)
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
            *super().get_param_report_fields(),
        ]

    @classmethod
    def get_result_report_fields(cls):
        return [
            TileCoordsF,
            NumSignifPairsF,
            NumModulesF,
            ModuleCoordsF,
            BestKsF,
            ClusterDirsF,
            *super().get_result_report_fields(),
        ]
