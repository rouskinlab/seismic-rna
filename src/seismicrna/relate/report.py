from __future__ import annotations

from functools import cache
from pathlib import Path

from .io import RelateIO, ReadNamesBatchIO, RelateBatchIO
from ..core import path
from ..core.io import RefIO
from ..core.report import (Report,
                           BatchedRefseqReport,
                           RefF,
                           SampleF,
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
                           ClipEnd3F)

BATCH_INDEX_COL = "Read Name"


class RelateReport(BatchedRefseqReport, RelateIO):

    @classmethod
    def fields(cls):
        return [SampleF,
                RefF,
                MinMapQualF,
                PhredEncF,
                MinPhredF,
                Insert3F,
                AmbindelF,
                OverhangsF,
                ClipEnd5F,
                ClipEnd3F,
                MinReadsF,
                NumReadsXamF,
                NumReadsRelF] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def _batch_types(cls):
        return [ReadNamesBatchIO, RelateBatchIO]

    def refseq_file(self, top: Path):
        return refseq_file_path(top,
                                self.get_field(SampleF),
                                self.get_field(RefF))


class PoolReport(Report, RefIO):

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def fields(cls):
        return [
            # Sample and reference.
            SampleF,
            RefF,
            # Pooled samples.
            PooledSamplesF,
        ] + super().fields()

    @classmethod
    def path_segs(cls):
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.RelateRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.RELATE_STEP}


@cache
def refseq_file_seg_types():
    return RelateReport.seg_types()[:-1] + (path.RefseqFileSeg,)


@cache
def refseq_file_auto_fields():
    return {**RelateReport.auto_fields(), path.EXT: path.BROTLI_PICKLE_EXT}


def refseq_file_path(top: Path, sample: str, ref: str):
    return path.build(*refseq_file_seg_types(),
                      **refseq_file_auto_fields(),
                      top=top,
                      sample=sample,
                      ref=ref)
