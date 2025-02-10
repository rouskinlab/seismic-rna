from abc import ABC

from .io import MaskIO, MaskBatchIO
from ..core import path
from ..core.join import JoinReport
from ..core.report import (RegReport,
                           BatchedReport,
                           End5F,
                           End3F,
                           CountMutsF,
                           CountRefsF,
                           ExclGUF,
                           ExclPolyAF,
                           ExclListPosF,
                           DiscontigF,
                           NumDiscontigF,
                           MinNInfoPosF,
                           MaxFMutPosF,
                           NumPosInitF,
                           NumPosCutGUF,
                           NumPosCutPolyAF,
                           NumPosCutListF,
                           NumPosCutLoInfoF,
                           NumPosCutHiMutF,
                           NumPosKeptF,
                           PosCutGUF,
                           PosCutPolyAF,
                           PosCutListF,
                           PosCutLoInfoF,
                           PosCutHiMutF,
                           PosKeptF,
                           MinNCovReadF,
                           MinFInfoReadF,
                           MaxFMutReadF,
                           MinMutGapF,
                           NumReadsInitF,
                           NumReadCutListF,
                           NumReadsLoNCovF,
                           NumReadsLoInfoF,
                           NumReadsHiMutF,
                           NumReadsCloseMutF,
                           NumReadsKeptF,
                           QuickUnbiasF,
                           QuickUnbiasThreshF,
                           MaxMaskIterF,
                           NumMaskIterF,
                           JoinedRegionsF)


class BaseMaskReport(RegReport, MaskIO, ABC):

    @classmethod
    def get_file_seg_type(cls):
        return path.MaskRepSeg


class MaskReport(BatchedReport, BaseMaskReport):

    @classmethod
    def _get_batch_types(cls):
        return [MaskBatchIO]

    @classmethod
    def get_param_report_fields(cls):
        return [
            # Region 5' and 3' ends.
            End5F,
            End3F,
            # Types of mutations and matches to count.
            CountMutsF,
            CountRefsF,
            # Position filtering parameters.
            ExclGUF,
            ExclPolyAF,
            ExclListPosF,
            MinNInfoPosF,
            MaxFMutPosF,
            # Read filtering parameters.
            MinNCovReadF,
            DiscontigF,
            MinFInfoReadF,
            MaxFMutReadF,
            MinMutGapF,
            # Iterations.
            MaxMaskIterF,
            # Observer bias correction.
            QuickUnbiasF,
            QuickUnbiasThreshF,
            *super().get_param_report_fields()
        ]

    @classmethod
    def get_result_report_fields(cls):
        return [
            # Position filtering results.
            NumPosInitF,
            NumPosCutGUF,
            NumPosCutPolyAF,
            NumPosCutListF,
            NumPosCutLoInfoF,
            NumPosCutHiMutF,
            NumPosKeptF,
            PosCutGUF,
            PosCutPolyAF,
            PosCutListF,
            PosCutLoInfoF,
            PosCutHiMutF,
            PosKeptF,
            # Read filtering results.
            NumReadsInitF,
            NumReadCutListF,
            NumReadsLoNCovF,
            NumDiscontigF,
            NumReadsLoInfoF,
            NumReadsHiMutF,
            NumReadsCloseMutF,
            NumReadsKeptF,
            # Iterations.
            NumMaskIterF,
            *super().get_result_report_fields()
        ]


class JoinMaskReport(JoinReport, BaseMaskReport):

    @classmethod
    def get_param_report_fields(cls):
        return [JoinedRegionsF,
                *super().get_param_report_fields()]
