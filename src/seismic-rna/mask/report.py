from ..core import path
from ..core.report import BatchReport


class MaskReport(BatchReport):
    __slots__ = (
        # Sample, reference, and section information.
        "sample", "ref", "sect", "end5", "end3",
        # Batch information.
        "checksums", "n_batches",
        # Types of mutations and matches to count.
        "count_muts", "count_refs",
        # Position filtering parameters.
        "exclude_gu", "exclude_polya", "exclude_pos",
        "min_ninfo_pos", "max_fmut_pos",
        # Position filtering results.
        "n_pos_init",
        "n_pos_gu", "n_pos_polya", "n_pos_user",
        "n_pos_min_ninfo", "n_pos_max_fmut", "n_pos_kept",
        "pos_gu", "pos_polya", "pos_user",
        "pos_min_ninfo", "pos_max_fmut", "pos_kept",
        # Read filtering parameters.
        "min_finfo_read", "max_fmut_read", "min_mut_gap",
        # Read filtering results.
        "n_reads_init",
        "n_reads_min_finfo", "n_reads_max_fmut", "n_reads_min_gap",
        "n_reads_kept",
    )

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.SectSeg, path.MaskRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.MOD: path.MOD_MASK}

    @classmethod
    def get_batch_seg(cls):
        return path.MaskBatSeg
