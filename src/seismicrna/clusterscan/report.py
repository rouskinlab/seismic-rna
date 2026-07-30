from ..core import path
from ..core.report import RegReport, DomainCoordsF, BestKsF, MergedDomainsF
from .io import ClusterScanIO


class ClusterScanReport(RegReport, ClusterScanIO):
    @classmethod
    def get_file_seg_type(cls):
        return path.ClusterScanRepSeg

    @classmethod
    def get_result_report_fields(cls):
        return [
            DomainCoordsF,
            BestKsF,
            MergedDomainsF,
            *super().get_result_report_fields(),
        ]
