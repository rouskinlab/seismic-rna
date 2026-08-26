from ..core import path
from ..core.report import (
    RegReport,
    # Region ends.
    End5F,
    End3F,
    # Domain-detection parameters.
    TileLengthF,
    TileMinOverlapF,
    EraseTilesF,
    BandWidthF,
    MinPairCoverageF,
    MinExpectBothF,
    PairFdrF,
    MinFoldChangeF,
    DetectFdrF,
    MergeFdrF,
    WidenF,
    FillF,
    MaxDomainLengthF,
    MinDomainLengthF,
    MinPairsF,
    # Results.
    TileCoordsF,
    NumPositivePairsF,
    NullBridgeRateF,
    NumDomainsF,
    DomainCoordsF,
)
from .io import FilterScanIO


class FilterScanReport(RegReport, FilterScanIO):
    @classmethod
    def get_file_seg_type(cls):
        return path.FilterScanRepSeg

    @classmethod
    def get_param_report_fields(cls):
        return [
            # Region 5' and 3' ends.
            End5F,
            End3F,
            TileLengthF,
            TileMinOverlapF,
            EraseTilesF,
            BandWidthF,
            MinPairCoverageF,
            MinExpectBothF,
            PairFdrF,
            MinFoldChangeF,
            DetectFdrF,
            MergeFdrF,
            WidenF,
            FillF,
            MaxDomainLengthF,
            MinDomainLengthF,
            MinPairsF,
            *super().get_param_report_fields(),
        ]

    @classmethod
    def get_result_report_fields(cls):
        return [
            TileCoordsF,
            NumPositivePairsF,
            NullBridgeRateF,
            NumDomainsF,
            DomainCoordsF,
            *super().get_result_report_fields(),
        ]
