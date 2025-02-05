from __future__ import annotations

from abc import ABC, abstractmethod

from ..core import path
from ..core.report import (Report,
                           SampleF,
                           RefF,
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


class AlignReport(Report, ABC):

    @classmethod
    @abstractmethod
    def fields(cls):
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
                AlignReadsInitF,
                ReadsTrimF,
                ReadsAlignF,
                ReadsDedupF,
                ReadsRefsF] + super().fields()

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.ALIGN_STEP}


class AlignSampleReport(AlignReport):

    @classmethod
    def fields(cls):
        return [SampleF] + super().fields()

    @classmethod
    def file_seg_type(cls):
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


class AlignRefReport(AlignReport):

    @classmethod
    def fields(cls):
        return [SampleF, RefF] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.AlignRefRepSeg

    def __init__(self, *,
                 ref: str,
                 demultiplexed: bool,
                 **kwargs):
        if not isinstance(ref, str):
            raise TypeError(f"Expected a reference name, but got {repr(ref)}")
        if not demultiplexed:
            raise ValueError(f"{type(self).__name__} must be demultiplexed")
        super().__init__(ref=ref, demultiplexed=demultiplexed, **kwargs)
