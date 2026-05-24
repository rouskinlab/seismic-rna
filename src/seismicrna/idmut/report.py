from __future__ import annotations

from abc import ABC
from pathlib import Path

from .io import IDmutIO, ReadNamesBatchIO, IDmutBatchIO, RefseqIO
from ..core import path
from ..core.report import (
    RefReport,
    BatchedReport,
    RefF,
    SampleF,
    BranchesF,
    PooledSamplesF,
    MinMapQualF,
    MinReadsF,
    NumReadsXamF,
    NumReadsRelF,
    Insert3F,
    AmbindelF,
    OverhangsF,
    MinPhredF,
    PhredEncF,
    ClipEnd5F,
    ClipEnd3F,
    RefseqChecksumF,
    MinPearsonPoolF,
    MaxMarcdPoolF,
)

BATCH_INDEX_COL = "Read Name"


class BaseIDmutReport(RefReport, IDmutIO, ABC):
    @classmethod
    def get_file_seg_type(cls):
        return path.IDmutRepSeg


class IDmutReport(BatchedReport, BaseIDmutReport):
    @classmethod
    def get_param_report_fields(cls):
        return [
            MinMapQualF,
            PhredEncF,
            MinPhredF,
            Insert3F,
            AmbindelF,
            OverhangsF,
            ClipEnd5F,
            ClipEnd3F,
            MinReadsF,
            *super().get_param_report_fields(),
        ]

    @classmethod
    def get_result_report_fields(cls):
        return [NumReadsXamF, NumReadsRelF, *super().get_result_report_fields()]

    @classmethod
    def get_checksum_report_fields(cls):
        return [*super().get_checksum_report_fields(), RefseqChecksumF]

    @classmethod
    def _get_batch_types(cls):
        return [ReadNamesBatchIO, IDmutBatchIO]

    def get_refseq_file(self, top: Path):
        return RefseqIO.build_path(
            {
                path.TOP: top,
                path.SAMPLE: self.get_field(SampleF),
                path.BRANCHES: self.get_path(BranchesF),
                path.REF: self.get_field(RefF),
            }
        )


class PoolReport(BaseIDmutReport):
    @classmethod
    def get_param_report_fields(cls):
        return [
            PooledSamplesF,
            MinPearsonPoolF,
            MaxMarcdPoolF,
            *super().get_param_report_fields(),
        ]
