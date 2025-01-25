from .io import MaskIO, MaskBatchIO
from ..core import path
from ..core.join import JoinReport
from ..core.report import (BatchedReport,
                           SampleF,
                           RefF,
                           RegF,
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


class MaskReport(BatchedReport, MaskIO):

    @classmethod
    def file_seg_type(cls):
        return path.MaskRepSeg

    @classmethod
    def _batch_types(cls):
        return [MaskBatchIO]

    @classmethod
    def fields(cls):
        return [
            # Sample, reference, and region information.
            SampleF,
            RefF,
            RegF,
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
            # Read filtering parameters.
            MinNCovReadF,
            DiscontigF,
            MinFInfoReadF,
            MaxFMutReadF,
            MinMutGapF,
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
            MaxMaskIterF,
            NumMaskIterF,
            # Observer bias correction.
            QuickUnbiasF,
            QuickUnbiasThreshF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.MASK_STEP}


class JoinMaskReport(JoinReport):

    @classmethod
    def file_seg_type(cls):
        return path.MaskRepSeg

    @classmethod
    def fields(cls):
        return [
            # Sample and reference.
            SampleF,
            RefF,
            RegF,
            # Joined data.
            JoinedRegionsF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.MASK_STEP}
