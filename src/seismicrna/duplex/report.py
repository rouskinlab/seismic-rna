from __future__ import annotations

from .io import DuplexIO
from ..core import path
from ..core.report import (
    RegReport,
    End5F,
    End3F,
    RefseqChecksumF,
    DuplexCutF,
    DuplexKsF,
    DuplexSourcesF,
)


class DuplexReport(RegReport, DuplexIO):
    """Report for a duplex of two source datasets into one duplex.

    A duplex report records only the fused identity, the two strands'
    combined sequence, the strand-break position, and pointers to the
    source tables; the per-position data live in the duplex's own CSV, and
    the read-level data stay in the source datasets (no duplication).
    """

    @classmethod
    def get_file_seg_type(cls):
        return path.DuplexRepSeg

    @classmethod
    def get_param_report_fields(cls):
        return [
            End5F,
            End3F,
            DuplexCutF,
            DuplexKsF,
            DuplexSourcesF,
            *super().get_param_report_fields(),
        ]

    @classmethod
    def get_checksum_report_fields(cls):
        return [RefseqChecksumF, *super().get_checksum_report_fields()]
