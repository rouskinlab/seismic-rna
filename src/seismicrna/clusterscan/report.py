from ..core import path
from ..core.report import RegReport, BestKsF, MergedDomainsF
from .io import ClusterScanIO


class ClusterScanReport(RegReport, ClusterScanIO):
    @classmethod
    def get_file_seg_type(cls):
        return path.ClusterScanRepSeg

    @classmethod
    def get_result_report_fields(cls):
        return [BestKsF, MergedDomainsF, *super().get_result_report_fields()]
