from __future__ import annotations
from abc import ABC
from functools import cached_property


from ..core.array import calc_inverse, get_length
from ..core.batch.read import ReadBatch
from ..core.batch.muts import RegionMutsBatch
from ..core.seq.region import Region
from ..core.validate import require_isinstance


class PartialReadBatch(ReadBatch, ABC):
    @cached_property
    def max_read(self):
        return self.read_nums.max(initial=0)

    @cached_property
    def read_indexes(self):
        return calc_inverse(self.read_nums, what="read_nums", verify=False)


class PartialRegionMutsBatch(PartialReadBatch, RegionMutsBatch, ABC):
    pass


class FilterReadBatch(PartialReadBatch):
    def __init__(self, *, read_nums: np.ndarray, **kwargs):
        import numpy as np

        require_isinstance("read_nums", read_nums, np.ndarray)
        self._read_nums = read_nums
        super().__init__(**kwargs)

    @property
    def read_nums(self):
        return self._read_nums

    @cached_property
    def num_reads(self):
        return get_length(self.read_nums, "read_nums")


class FilterMutsBatch(FilterReadBatch, PartialRegionMutsBatch):
    @property
    def read_weights(self):
        # Masked reads are excluded by removing them upstream, not by
        # weighting. The weighted counting path requires a per-cluster
        # Series num_reads, which this batch (scalar num_reads) does not
        # provide, so per-read weighting of masked reads is unsupported.
        if self.masked_reads_bool.any():
            raise NotImplementedError(
                "Per-read weighting of masked reads is not supported for "
                f"{type(self).__name__}; remove masked reads instead"
            )
        return None


def apply_filters(
    batch: RegionMutsBatch,
    read_nums: np.ndarray | None = None,
    region: Region | None = None,
    sanitize: bool = False,
):
    """Apply read/position filters to a batch, returning a FilterMutsBatch.

    Parameters
    ----------
    batch: RegionMutsBatch
        Source batch to filter.
    read_nums: np.ndarray or None, optional
        Array of read numbers to keep; if None, all reads are kept.
    region: Region or None, optional
        Region to clip reads to; if None, the batch's existing region
        is used.
    sanitize: bool, optional
        Whether to run extra validation checks on the new batch
        (default False).

    Returns
    -------
    FilterMutsBatch
        A new batch containing only the filtered reads and positions.
    """
    import numpy as np

    require_isinstance("batch", batch, RegionMutsBatch)
    # Determine which reads to use.
    if read_nums is not None:
        # Drop the reads that do not appear in read_nums.
        read_indexes = batch.read_indexes[read_nums]
        seg_end5s = batch.seg_end5s[read_indexes]
        seg_end3s = batch.seg_end3s[read_indexes]
        # Determine which reads were dropped.
        dropped_reads = np.setdiff1d(batch.read_nums, read_nums)
    else:
        # Use all reads.
        read_nums = batch.read_nums
        seg_end5s = batch.seg_end5s
        seg_end3s = batch.seg_end3s
        dropped_reads = None
    # Determine the region of the new batch.
    if region is not None:
        require_isinstance("region", region, Region)
        # Clip the read coordinates to the region bounds.
        seg_end5s = seg_end5s.clip(region.end5, region.end3 + 1)
        seg_end3s = seg_end3s.clip(region.end5 - 1, region.end3)
    else:
        # Use the same region as the given batch.
        region = batch.region
    # Put only the remaining positions and reads into muts.
    muts = {
        pos: (
            {
                mut: np.setdiff1d(pos_mut_reads, dropped_reads, assume_unique=True)
                for mut, pos_mut_reads in batch.muts[pos].items()
            }
            if dropped_reads is not None
            else batch.muts[pos]
        )
        for pos in region.unmasked_int
    }
    return FilterMutsBatch(
        batch=batch.batch,
        read_nums=read_nums,
        region=region,
        seg_end5s=seg_end5s,
        seg_end3s=seg_end3s,
        muts=muts,
        sanitize=sanitize,
    )
