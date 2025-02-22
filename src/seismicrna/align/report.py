from __future__ import annotations

from abc import ABC

from ..core import path
from ..core.report import (Report,
                           RefF,
                           RefFastaChecksumF,
                           FastqChecksumsF,
                           XamChecksumF,
                           IsDemultF,
                           IsPairedEndF,
                           PhredEncF,
                           UseFastpF,
                           Fastp5F,
                           Fastp3F,
                           FastpWF,
                           FastpMF,
                           FastpPolyGF,
                           FastpPolyGMinLenF,
                           FastpPolyXF,
                           FastpPolyXMinLenF,
                           FastpAdapterTrimmingF,
                           FastpAdapter1F,
                           FastpAdapter2F,
                           FastpAdapterFastaF,
                           FastpDetectAdapterForPEF,
                           FastpMinLengthF,
                           Bowtie2Local,
                           Bowtie2Discord,
                           Bowtie2Mixed,
                           Bowtie2Dovetail,
                           Bowtie2Contain,
                           Bowtie2Un,
                           Bowtie2ScoreMin,
                           Bowtie2MinLengthF,
                           Bowtie2MaxLengthF,
                           Bowtie2GBarF,
                           Bowtie2SeedLength,
                           Bowtie2SeedInterval,
                           Bowtie2ExtTries,
                           Bowtie2Reseed,
                           Bowtie2Dpad,
                           Bowtie2Orient,
                           SepStrandsF,
                           F1R2FwdF,
                           RevLabelF,
                           MinMapQualF,
                           MinReadsF,
                           AlignReadsInitF,
                           ReadsTrimF,
                           ReadsAlignF,
                           ReadsDedupF,
                           ReadsRefsF)


class BaseAlignReport(Report, ABC):

    @classmethod
    def get_step(cls):
        return path.ALIGN_STEP

    @classmethod
    def get_param_report_fields(cls):
        return [IsDemultF,
                IsPairedEndF,
                PhredEncF,
                UseFastpF,
                Fastp5F,
                Fastp3F,
                FastpWF,
                FastpMF,
                FastpPolyGF,
                FastpPolyGMinLenF,
                FastpPolyXF,
                FastpPolyXMinLenF,
                FastpAdapterTrimmingF,
                FastpAdapter1F,
                FastpAdapter2F,
                FastpAdapterFastaF,
                FastpDetectAdapterForPEF,
                FastpMinLengthF,
                Bowtie2Local,
                Bowtie2Discord,
                Bowtie2Mixed,
                Bowtie2Dovetail,
                Bowtie2Contain,
                Bowtie2ScoreMin,
                Bowtie2MinLengthF,
                Bowtie2MaxLengthF,
                Bowtie2GBarF,
                Bowtie2SeedLength,
                Bowtie2SeedInterval,
                Bowtie2ExtTries,
                Bowtie2Reseed,
                Bowtie2Dpad,
                Bowtie2Orient,
                Bowtie2Un,
                MinMapQualF,
                SepStrandsF,
                F1R2FwdF,
                RevLabelF,
                MinReadsF,
                *super().get_param_report_fields()]

    @classmethod
    def get_result_report_fields(cls):
        return [AlignReadsInitF,
                ReadsTrimF,
                ReadsAlignF,
                ReadsDedupF,
                ReadsRefsF,
                *super().get_result_report_fields()]

    @classmethod
    def get_checksum_report_fields(cls):
        return [RefFastaChecksumF,
                FastqChecksumsF]


class AlignSampleReport(BaseAlignReport):

    @classmethod
    def get_file_seg_type(cls):
        return path.AlignSampleRepSeg

    def __init__(self, *,
                 ref: str | None = None,
                 demultiplexed: bool,
                 **kwargs):
        if ref is not None:
            raise TypeError(f"Got an unexpected reference name: {repr(ref)}")
        if demultiplexed:
            raise ValueError(f"{type(self).__name__} cannot be demultiplexed")
        super().__init__(demultiplexed=demultiplexed, **kwargs)


class AlignRefReport(BaseAlignReport):

    @classmethod
    def get_file_seg_type(cls):
        return path.AlignRefRepSeg

    @classmethod
    def get_ident_report_fields(cls):
        return [*super().get_ident_report_fields(),
                RefF]

    def __init__(self, *,
                 ref: str,
                 demultiplexed: bool,
                 **kwargs):
        if not isinstance(ref, str):
            raise TypeError(f"Expected a reference name, but got {repr(ref)}")
        if not demultiplexed:
            raise ValueError(f"{type(self).__name__} must be demultiplexed")
        super().__init__(ref=ref, demultiplexed=demultiplexed, **kwargs)


class SplitReport(Report, ABC):

    @classmethod
    def get_file_seg_type(cls):
        return path.SplitRepSeg

    @classmethod
    def get_step(cls):
        return path.ALIGN_STEP

    @classmethod
    def get_param_report_fields(cls):
        return [Bowtie2Local,
                Bowtie2Discord,
                Bowtie2Mixed,
                Bowtie2Dovetail,
                Bowtie2Contain,
                Bowtie2ScoreMin,
                Bowtie2MinLengthF,
                Bowtie2MaxLengthF,
                Bowtie2GBarF,
                Bowtie2SeedLength,
                Bowtie2SeedInterval,
                Bowtie2ExtTries,
                Bowtie2Reseed,
                Bowtie2Dpad,
                Bowtie2Orient,
                Bowtie2Un,
                MinMapQualF,
                SepStrandsF,
                F1R2FwdF,
                RevLabelF,
                MinReadsF,
                *super().get_param_report_fields()]

    @classmethod
    def get_result_report_fields(cls):
        return [ReadsRefsF,
                *super().get_result_report_fields()]

    @classmethod
    def get_checksum_report_fields(cls):
        return [RefFastaChecksumF,
                XamChecksumF]
