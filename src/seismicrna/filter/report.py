from abc import ABC

from .io import FilterIO, FilterBatchIO
from ..core import path
from ..core.join import JoinReport
from ..core.report import (
    RegReport,
    BatchedReport,
    ProbeF,
    End5F,
    End3F,
    CountMutsF,
    CountRefsF,
    MaskAF,
    MaskCF,
    MaskGF,
    MaskUF,
    MaskPolyAF,
    MaskPosF,
    DiscontigF,
    NumDiscontigF,
    MinNInfoPosF,
    MaxFMutPosF,
    NumPosInitF,
    NumPosCutAF,
    NumPosCutCF,
    NumPosCutGF,
    NumPosCutUF,
    NumPosCutNF,
    NumPosCutPolyAF,
    NumPosCutListF,
    NumPosCutLoInfoF,
    NumPosCutHiMutF,
    NumPosKeptF,
    PosCutAF,
    PosCutCF,
    PosCutGF,
    PosCutUF,
    PosCutNF,
    PosCutPolyAF,
    PosCutListF,
    PosCutLoInfoF,
    PosCutHiMutF,
    PosKeptF,
    MinNCovReadF,
    MinFCovReadF,
    MinFInfoReadF,
    MaxFMutReadF,
    MinMutGapF,
    MutCollisionsF,
    NumReadsInitF,
    NumReadCutListF,
    NumReadsLoNCovF,
    NumReadsLoFCovF,
    NumReadsLoInfoF,
    NumReadsHiMutF,
    NumReadsCloseMutF,
    NumReadsKeptF,
    QuickUnbiasF,
    QuickUnbiasThreshF,
    MaxFilterIterF,
    NumFilterIterF,
    JoinedRegionsF,
)


class BaseFilterReport(RegReport, FilterIO, ABC):
    @classmethod
    def get_file_seg_type(cls):
        return path.FilterRepSeg


class FilterReport(BatchedReport, BaseFilterReport):
    @classmethod
    def _get_batch_types(cls):
        return [FilterBatchIO]

    @classmethod
    def get_param_report_fields(cls):
        return [
            # Region 5' and 3' ends.
            End5F,
            End3F,
            # Probe type.
            ProbeF,
            # Types of mutations and matches to count.
            CountMutsF,
            CountRefsF,
            # Position filtering parameters.
            MaskAF,
            MaskCF,
            MaskGF,
            MaskUF,
            MaskPolyAF,
            MaskPosF,
            MinNInfoPosF,
            MaxFMutPosF,
            # Read filtering parameters.
            MinNCovReadF,
            MinFCovReadF,
            DiscontigF,
            MinFInfoReadF,
            MaxFMutReadF,
            MinMutGapF,
            MutCollisionsF,
            # Iterations.
            MaxFilterIterF,
            # Observer bias correction.
            QuickUnbiasF,
            QuickUnbiasThreshF,
            *super().get_param_report_fields(),
        ]

    @classmethod
    def get_result_report_fields(cls):
        return [
            # Position filtering results.
            NumPosInitF,
            NumPosCutAF,
            NumPosCutCF,
            NumPosCutGF,
            NumPosCutUF,
            NumPosCutNF,
            NumPosCutPolyAF,
            NumPosCutListF,
            NumPosCutLoInfoF,
            NumPosCutHiMutF,
            NumPosKeptF,
            PosCutAF,
            PosCutCF,
            PosCutGF,
            PosCutUF,
            PosCutNF,
            PosCutPolyAF,
            PosCutListF,
            PosCutLoInfoF,
            PosCutHiMutF,
            PosKeptF,
            # Read filtering results.
            NumReadsInitF,
            NumReadCutListF,
            NumReadsLoNCovF,
            NumReadsLoFCovF,
            NumDiscontigF,
            NumReadsLoInfoF,
            NumReadsHiMutF,
            NumReadsCloseMutF,
            NumReadsKeptF,
            # Iterations.
            NumFilterIterF,
            *super().get_result_report_fields(),
        ]


class JoinFilterReport(JoinReport, BaseFilterReport):
    @classmethod
    def get_param_report_fields(cls):
        return [JoinedRegionsF, *super().get_param_report_fields()]
