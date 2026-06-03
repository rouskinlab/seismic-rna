import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .batch import FilterMutsBatch, apply_filters
from .dataset import FilterMutsDataset
from .io import FilterBatchIO
from .report import FilterReport
from .table import FilterBatchTabulator
from ..core import path
from ..core.arg import MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE, PROBE_NONE
from ..core.batch import RegionMutsBatch
from ..core.dataset import MissingBatchTypeError
from ..core.error import IncompatibleValuesError
from ..core.lists import PositionList
from ..core.logs import logger
from ..core.rel import RelPattern
from ..core.report import filter_iter_no_convergence
from ..core.seq import (
    BASEA,
    BASEC,
    BASEG,
    BASET,
    BASEN,
    FIELD_REF,
    POS_NAME,
    Region,
    index_to_pos,
)
from ..core.table import MUTAT_REL, INFOR_REL
from ..core.tmp import release_to_out, with_tmp_dir
from ..core.write import need_write
from ..idmut.dataset import IDmutMutsDataset, PoolMutsDataset, load_read_names_dataset


class Filterer(object):
    """Filter batches of relationships."""

    PATTERN_KEY = "pattern"
    DROP_READ_INIT = "read-init"
    DROP_READ_LIST = "read-exclude"
    DROP_READ_NCOV = "read-ncov"
    DROP_READ_FCOV = "read-fcov"
    DROP_READ_DISCONTIG = "read-discontig"
    DROP_READ_FINFO = "read-finfo"
    DROP_READ_FMUT = "read-fmut"
    DROP_READ_GAP = "read-gap"
    DROP_READ_KEPT = "read-kept"
    MASK_POS_NINFO = "pos-ninfo"
    MASK_POS_FMUT = "pos-fmut"
    CHECKSUM_KEY = FilterReport.get_batch_type().btype()

    def __init__(
        self,
        dataset: IDmutMutsDataset | PoolMutsDataset,
        region: Region,
        pattern: RelPattern,
        *,
        max_filter_iter: int,
        mask_polya: int,
        mask_a: bool,
        mask_c: bool,
        mask_g: bool,
        mask_u: bool,
        mask_pos: list[tuple[str, int]],
        mask_pos_file: list[Path],
        drop_read: list[str],
        drop_read_file: list[Path],
        drop_discontig: bool,
        min_ncov_read: int,
        min_fcov_read: float,
        min_finfo_read: float,
        max_fmut_read: float,
        min_mut_gap: int,
        mut_collisions: str,
        probe: str,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        quick_unbias: bool,
        quick_unbias_thresh: float,
        count_read: bool,
        brotli_level: int,
        top: Path,
        branch: str,
        self_contained: bool,
        num_cpus: int = 1,
    ):
        """
        Parameters
        ----------
        dataset: IDmutMutsDataset or PoolMutsDataset
            Dataset of read-level mutation data to be filtered.
        region: Region
            Region to filter over.
        pattern: RelPattern
            Relationship pattern defining which bases count as
            mutations.
        max_filter_iter: int
            Maximum number of filtering iterations; 0 means unlimited.
        mask_polya: int
            Mask positions in poly(A) runs of at least this length.
        mask_a: bool
            Whether to mask adenine positions.
        mask_c: bool
            Whether to mask cytosine positions.
        mask_g: bool
            Whether to mask guanine positions.
        mask_u: bool
            Whether to mask uracil/thymine positions.
        mask_pos: list[tuple[str, int]]
            ``(ref, position)`` pairs to pre-exclude from the region.
        mask_pos_file: list[Path]
            Files of positions to pre-exclude.
        drop_read: list[str]
            Read names to pre-exclude.
        drop_read_file: list[Path]
            Files of read names to pre-exclude.
        drop_discontig: bool
            Whether to drop reads with discontiguous mate pairs.
        min_ncov_read: int
            Minimum number of covered positions required per read.
        min_fcov_read: float
            Minimum fraction of region positions that must be covered per read.
        min_finfo_read: float
            Minimum fraction of informative positions required per read.
        max_fmut_read: float
            Maximum fraction of mutated positions allowed per read.
        min_mut_gap: int
            Minimum gap in nucleotides between adjacent mutations in a
            read; reads with closer mutations are handled per
            ``mut_collisions``.
        mut_collisions: str
            How to handle reads with mutations closer than
            ``min_mut_gap`` (e.g. drop or merge).
        probe: str
            Probe type used for auto-selecting defaults.
        min_ninfo_pos: int
            Minimum number of informative reads required per position.
        max_fmut_pos: float
            Maximum mutation fraction allowed per position.
        quick_unbias: bool
            Whether to use the fast approximate bias-correction.
        quick_unbias_thresh: float
            Convergence threshold for the quick unbias algorithm.
        count_read: bool
            Whether to count reads for the read-level table.
        brotli_level: int
            Brotli compression level for batch files.
        top: Path
            Top-level output directory.
        branch: str
            Branch label appended to the filter step in output paths.
        num_cpus: int, optional
            Number of CPUs for parallel batch processing (default 1).
        """
        # Set the general parameters.
        self._began = datetime.now()
        self.dataset = dataset
        self.region = Region(
            dataset.ref,
            dataset.refseq,
            end5=region.end5,
            end3=region.end3,
            name=region.name,
        )
        self.pattern = pattern
        self.max_iter = max_filter_iter
        self._iter = 0
        self._converged = False
        # Set the parameters for excluding positions from the region.
        if probe != PROBE_NONE and not 0 < mask_polya <= 5:
            logger.warning(
                "It is not recommended to keep sequences of 5 or "
                "more consecutive As because of an artifact during "
                "RT that causes low reactivity. See Kladwang et al. "
                "(https://doi.org/10.1021/acs.biochem.0c00020)."
            )
        self.mask_polya = mask_polya
        self.mask_a = mask_a
        self.mask_c = mask_c
        self.mask_g = mask_g
        self.mask_u = mask_u
        self.mask_pos = self._get_mask_pos(mask_pos, mask_pos_file)
        self.drop_read = self._get_drop_read(drop_read, drop_read_file)
        # Set the parameters for filtering reads.
        if min_mut_gap > 0 and not drop_discontig:
            raise ValueError(
                "The observer bias correction does not work with "
                "discontiguous reads. If you need discontiguous "
                "reads, disable bias correction with the option "
                "--min-mut-gap=0 (but be warned that disabling "
                "bias correction can produce misleading results, "
                "especially with clustering)."
            )
        self.drop_discontig = drop_discontig
        self.min_ncov_read = min_ncov_read
        self.min_fcov_read = min_fcov_read
        self.min_mut_gap = min_mut_gap
        self.mut_collisions = mut_collisions
        self.probe = probe
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        self._n_reads = defaultdict(int)
        # Set the parameters for filtering positions.
        self.min_ninfo_pos = min_ninfo_pos
        self.max_fmut_pos = max_fmut_pos
        # Set the parameters for observer bias correction.
        self.quick_unbias = quick_unbias
        self.quick_unbias_thresh = quick_unbias_thresh
        # Set the parameters for saving files.
        self.top = top
        self.count_read = count_read
        self.brotli_level = brotli_level
        self.self_contained = self_contained
        self.checksums = [""] * dataset.num_batches
        # After the first iteration, self.dataset will become the new,
        # filtered dataset, which will also have a branch for the filter
        # step, so calculate the branches using the original dataset.
        self.branches = path.add_branch(path.FILTER_STEP, branch, dataset.branches)
        # Parallelization
        self.num_cpus = num_cpus

    # This property can change: do not cache it.
    @property
    def n_reads_init(self):
        return self._n_reads[self.DROP_READ_INIT]

    # This property can change: do not cache it.
    @property
    def n_reads_list(self):
        return self._n_reads[self.DROP_READ_LIST]

    # This property can change: do not cache it.
    @property
    def n_reads_min_ncov(self):
        return self._n_reads[self.DROP_READ_NCOV]

    # This property can change: do not cache it.
    @property
    def n_reads_min_fcov(self):
        return self._n_reads[self.DROP_READ_FCOV]

    # This property can change: do not cache it.
    @property
    def n_reads_discontig(self):
        return self._n_reads[self.DROP_READ_DISCONTIG]

    # This property can change: do not cache it.
    @property
    def n_reads_min_finfo(self):
        return self._n_reads[self.DROP_READ_FINFO]

    # This property can change: do not cache it.
    @property
    def n_reads_max_fmut(self):
        return self._n_reads[self.DROP_READ_FMUT]

    # This property can change: do not cache it.
    @property
    def n_reads_min_gap(self):
        return self._n_reads[self.DROP_READ_GAP]

    # This property can change: do not cache it.
    @property
    def n_reads_kept(self):
        """Number of reads kept."""
        return self._n_reads[self.DROP_READ_KEPT]

    # This property can change: do not cache it.
    @property
    def pos_a(self):
        """Positions masked for having base A."""
        return self.region.get_mask(BASEA, missing_ok=True)

    # This property can change: do not cache it.
    @property
    def pos_c(self):
        """Positions masked for having base C."""
        return self.region.get_mask(BASEC, missing_ok=True)

    # This property can change: do not cache it.
    @property
    def pos_g(self):
        """Positions masked for having base G."""
        return self.region.get_mask(BASEG, missing_ok=True)

    # This property can change: do not cache it.
    @property
    def pos_u(self):
        """Positions masked for having base T or U."""
        return self.region.get_mask(BASET, missing_ok=True)

    # This property can change: do not cache it.
    @property
    def pos_n(self):
        """Positions masked for having base N."""
        return self.region.get_mask(BASEN, missing_ok=True)

    # This property can change: do not cache it.
    @property
    def pos_polya(self):
        """Positions masked for lying in a poly(A) sequence."""
        return self.region.get_mask(self.region.MASK_POLYA)

    # This property can change: do not cache it.
    @property
    def pos_list(self):
        """Positions masked arbitrarily from a list."""
        return self.region.get_mask(self.region.MASK_LIST)

    # This property can change: do not cache it.
    @property
    def pos_min_ninfo(self):
        """Positions masked for having too few informative reads."""
        return self.region.get_mask(self.MASK_POS_NINFO)

    # This property can change: do not cache it.
    @property
    def pos_max_fmut(self):
        """Positions masked for having too many mutations."""
        return self.region.get_mask(self.MASK_POS_FMUT)

    # This property can change: do not cache it.
    @property
    def pos_kept(self):
        """Positions kept."""
        return self.region.unmasked_int

    # This property can change: do not cache it.
    @property
    def _force_write(self):
        """Whether to force-write each file."""
        return self._iter > 1

    # This property can change: do not cache it.
    @property
    def read_names_dataset(self):
        """Dataset of the read names."""
        return load_read_names_dataset(self.dataset.report_file)

    def _get_mask_pos(
        self, mask_pos: Iterable[tuple[str, int]], mask_pos_file: Iterable[str | Path]
    ):
        """List all positions to mask."""
        # Collect the positions to mask from the list.
        dataset_ref = self.dataset.ref
        mask_pos = np.array(
            [pos for ref, pos in mask_pos if ref == dataset_ref], dtype=int
        )
        logger.trace(
            f"Got {mask_pos.size} positions listed individually "
            f"to pre-exclude for reference {repr(dataset_ref)}"
        )
        # List positions to exclude from file(s).
        for file in path.find_files_chain(mask_pos_file, [path.PositionListSeg]):
            file_data = PositionList.load_data(file)
            ref_rows = file_data[FIELD_REF] == dataset_ref
            file_pos = file_data.loc[ref_rows, POS_NAME].values
            logger.trace(
                f"Got {file_pos.size} positions in {file} "
                f"to pre-exclude for reference {repr(dataset_ref)}"
            )
            if file_pos.size > 0:
                mask_pos = np.concatenate([mask_pos, file_pos])
        # Drop redundant positions and sort the remaining ones.
        mask_pos = np.unique(np.asarray(mask_pos, dtype=int))
        # Keep only the positions in the region.
        mask_pos = mask_pos[
            np.logical_and(mask_pos >= self.region.end5, mask_pos <= self.region.end3)
        ]
        logger.trace(
            f"Got {mask_pos.size} positions to pre-exclude "
            f"for reference {repr(dataset_ref)}"
        )
        return mask_pos

    @staticmethod
    def _get_drop_read(drop_read: Iterable[str], drop_read_file: Iterable[str | Path]):
        """List all reads to drop."""
        # Ensure that the given read names are all strings.
        drop_read = set(map(str, drop_read))
        # List reads to exclude from file(s).
        for file in path.find_files_chain(drop_read_file, [path.ReadListSeg]):
            with open(file) as f:
                # Ensure every read name is unique.
                drop_read.update(map(str.rstrip, f))
        drop_read = np.asarray(list(drop_read), dtype=str)
        logger.trace(f"Got {drop_read.size} reads to pre-exclude")
        return drop_read

    def _drop_predefined_reads(self, batch: RegionMutsBatch):
        """Drop reads from a predefined list."""
        if self.drop_read.size == 0:
            # Pre-exclude no reads.
            logger.trace(f"{self} skipped pre-excluding reads in {batch}")
            return batch
        # Load the names of the reads in this batch.
        try:
            names_batch = self.read_names_dataset.get_batch(batch.batch)
        except MissingBatchTypeError:
            raise IncompatibleValuesError(
                "Reads can be excluded with --drop-read or --drop-read-file "
                "only if IDmut was run using the option --write-read-names"
            ) from None
        if names_batch.num_reads != batch.num_reads:
            raise ValueError(
                f"Expected {batch.num_reads} read names,but got {names_batch.num_reads}"
            )
        # Find the numbers of the reads to keep.
        reads = batch.read_nums[
            np.isin(names_batch.names, self.drop_read, assume_unique=True, invert=True)
        ]
        logger.trace(f"{self} kept {reads.size} reads after pre-excluding")
        # Return a new batch of only those reads.
        return apply_filters(batch, reads)

    def _drop_min_ncov_read(self, batch: RegionMutsBatch):
        """Drop reads with too few covered positions."""
        if self.min_ncov_read < 1:
            raise ValueError(f"min_ncov_read must be ≥ 1, but got {self.min_ncov_read}")
        # Find the reads with sufficiently many covered positions.
        reads = batch.read_nums[
            batch.cover_per_read.values.sum(axis=1) >= self.min_ncov_read
        ]
        logger.trace(
            f"{self} kept {reads.size} reads with coverage "
            f"≥ {self.min_ncov_read} in {batch}"
        )
        # Return a new batch of only those reads.
        return apply_filters(batch, reads)

    def _drop_min_fcov_read(self, batch: RegionMutsBatch):
        """Drop reads covering too small a fraction of the region."""
        if not 0.0 <= self.min_fcov_read <= 1.0:
            raise ValueError(
                f"min_fcov_read must be in [0, 1], but got {self.min_fcov_read}"
            )
        if self.min_fcov_read == 0.0:
            logger.trace(
                f"{self} skipped filtering reads with insufficient "
                f"coverage fractions in {batch}"
            )
            return batch
        ncov = batch.cover_per_read.values.sum(axis=1)
        n_pos = self.region.size
        fcov = ncov / n_pos if n_pos > 0 else np.zeros(ncov.size)
        reads = batch.read_nums[fcov >= self.min_fcov_read]
        logger.trace(
            f"{self} kept {reads.size} reads with coverage "
            f"fractions ≥ {self.min_fcov_read} in {batch}"
        )
        return apply_filters(batch, reads)

    def _drop_discontig(self, batch: RegionMutsBatch):
        """Drop reads with discontiguous mates."""
        if not self.drop_discontig:
            # Keep discontiguous reads.
            logger.trace(
                f"{self} skipped filtering reads with discontiguous mates in {batch}"
            )
            return batch
        # Find the reads with contiguous mates.
        reads = batch.read_nums[batch.contiguous]
        logger.trace(f"{self} kept {reads.size} reads with contiguous mates in {batch}")
        # Return a new batch of only those reads.
        return apply_filters(batch, reads)

    def _drop_min_finfo_read(self, batch: RegionMutsBatch):
        """Drop reads with too few informative positions."""
        if not 0.0 <= self.min_finfo_read <= 1.0:
            raise ValueError(
                f"min_finfo_read must be ≥ 0, ≤ 1, but got {self.min_finfo_read}"
            )
        if self.min_finfo_read == 0.0:
            # All reads have sufficiently many informative positions.
            logger.trace(
                f"{self} skipped filtering reads with insufficient "
                f"informative fractions in {batch}"
            )
            return batch
        # Find the reads with sufficiently many informative positions.
        info, muts = batch.count_per_read(self.pattern)
        finfo_read = info.values / batch.cover_per_read.values.sum(axis=1)
        reads = info.index.values[finfo_read >= self.min_finfo_read]
        logger.trace(
            f"{self} kept {reads.size} reads with informative "
            f"fractions ≥ {self.min_finfo_read} in {batch}"
        )
        # Return a new batch of only those reads.
        return apply_filters(batch, reads)

    def _drop_max_fmut_read(self, batch: RegionMutsBatch):
        """Drop reads with too many mutations."""
        if not 0.0 <= self.max_fmut_read <= 1.0:
            raise ValueError(
                f"max_fmut_read must be ≥ 0, ≤ 1, but got {self.max_fmut_read}"
            )
        if self.max_fmut_read == 1.0:
            # All reads have sufficiently few mutations.
            logger.trace(
                f"{self} skipped filtering reads with excessive "
                f"mutation fractions in {batch}"
            )
            return batch
        # Find the reads with sufficiently few mutations.
        info, muts = batch.count_per_read(self.pattern)
        with np.errstate(invalid="ignore"):
            fmut_read = muts.values / info.values
        reads = info.index.values[fmut_read <= self.max_fmut_read]
        logger.trace(
            f"{self} kept {reads.size} reads with mutated "
            f"fractions ≤ {self.max_fmut_read} in {batch}"
        )
        # Return a new batch of only those reads.
        return apply_filters(batch, reads)

    def _filter_min_mut_gap(self, batch: RegionMutsBatch):
        """Filter out mutations that are too close by either dropping
        reads or merging mutations."""
        if not self.min_mut_gap >= 0:
            raise ValueError(f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_mut_gap == 0:
            # No read can have a pair of mutations that are too close.
            logger.trace(
                f"{self} skipped filtering pairs of mutations too close in {batch}"
            )
            return batch
        if self.mut_collisions == MUT_COLLISIONS_DROP:
            # Drop reads with pairs of mutations that are too close.
            reads = batch.reads_noclose_muts(self.pattern, self.min_mut_gap)
            logger.trace(
                f"{self} kept {reads.size} reads with no two mutations "
                f"separated by < {self.min_mut_gap} nt in {batch}"
            )
            return apply_filters(batch, reads)
        if self.mut_collisions == MUT_COLLISIONS_MERGE:
            # Merge nearby mutations into a single mutation.
            logger.trace(
                f"{self} merged mutations closer than {self.min_mut_gap} nt in {batch}"
            )
            return FilterMutsBatch(
                batch=batch.batch,
                read_nums=batch.read_nums,
                region=batch.region,
                seg_end5s=batch.seg_end5s,
                seg_end3s=batch.seg_end3s,
                muts=batch.merge_close_muts(self.pattern, self.min_mut_gap),
            )
        raise ValueError(f"Invalid mut_collisions: {repr(self.mut_collisions)}")

    def _mask_predefined_positions(self):
        """Mask predefined positions."""
        self.region.mask_n()
        if self.mask_a:
            self.region.mask_a()
        if self.mask_c:
            self.region.mask_c()
        if self.mask_g:
            self.region.mask_g()
        if self.mask_u:
            self.region.mask_t()
        self.region.mask_polya(self.mask_polya)
        self.region.mask_list(self.mask_pos)

    def _get_batch_num_path(self, batch_num: int):
        return FilterBatchIO.build_path(
            {
                path.TOP: self.top,
                path.SAMPLE: self.dataset.sample,
                path.BRANCHES: self.branches,
                path.REF: self.dataset.ref,
                path.REG: self.region.name,
                path.BATCH: batch_num,
            }
        )

    def _get_n_reads_path(self, batch_num: int):
        """Get the path to the file of the number of reads dropped
        for a batch."""
        return self._get_batch_num_path(batch_num).with_suffix(path.JSON_EXT)

    def _get_checksum_path(self, batch_num: int):
        """Get the path to the file of the checksum for a batch."""
        return self._get_batch_num_path(batch_num).with_suffix(path.TXT_EXT)

    def _filter_batch_reads(self, batch_num: int, **kwargs):
        """Drop the reads in the batch that do not pass the filters
        and return a new batch without those reads."""
        n_reads = dict()
        # Load the batch.
        batch = self.dataset.get_batch(batch_num)
        # Determine the initial number of reads in the batch.
        n_reads[self.DROP_READ_INIT] = (n := batch.num_reads)
        if self._iter == 1:
            # Mask positions in the dataset that are masked in the region.
            batch = apply_filters(batch, region=self.region)
            # Drop reads from a predefined list.
            batch = self._drop_predefined_reads(batch)
            n_reads[self.DROP_READ_LIST] = n - (n := batch.num_reads)
        # Drop reads with too few covered positions.
        batch = self._drop_min_ncov_read(batch)
        n_reads[self.DROP_READ_NCOV] = n - (n := batch.num_reads)
        # Drop reads covering too small a fraction of the region.
        batch = self._drop_min_fcov_read(batch)
        n_reads[self.DROP_READ_FCOV] = n - (n := batch.num_reads)
        # Drop reads with discontiguous mates.
        batch = self._drop_discontig(batch)
        n_reads[self.DROP_READ_DISCONTIG] = n - (n := batch.num_reads)
        # Drop reads with too few informative positions.
        batch = self._drop_min_finfo_read(batch)
        n_reads[self.DROP_READ_FINFO] = n - (n := batch.num_reads)
        # Drop reads with too many mutations.
        batch = self._drop_max_fmut_read(batch)
        n_reads[self.DROP_READ_FMUT] = n - (n := batch.num_reads)
        # Filter out mutations that are too close together.
        batch = self._filter_min_mut_gap(batch)
        n_reads[self.DROP_READ_GAP] = n - (n := batch.num_reads)
        # Record the number of reads remaining after filtering.
        n_reads[self.DROP_READ_KEPT] = n
        # Save the batch.
        sc_kwargs = (
            dict(
                region=batch.region,
                seg_end5s=batch.seg_end5s,
                seg_end3s=batch.seg_end3s,
                muts=batch.muts,
            )
            if self.self_contained
            else {}
        )
        batch_file = FilterBatchIO(
            sample=self.dataset.sample,
            branches=self.branches,
            ref=self.dataset.ref,
            reg=self.region.name,
            batch=batch.batch,
            read_nums=batch.read_nums,
            **sc_kwargs,
        )
        _, checksum = batch_file.save(
            self.top, brotli_level=self.brotli_level, force=self._force_write
        )
        # Save the checksum.
        checksum_file = self._get_checksum_path(batch_num)
        with open(checksum_file, "x") as f:
            f.write(checksum)
            logger.trace(f"Wrote checksum {repr(checksum)} to {checksum_file}")
        # Save the read counts.
        n_reads_file = self._get_n_reads_path(batch_num)
        with open(n_reads_file, "x") as f:
            json.dump(n_reads, f)
            logger.trace(f"Wrote number of reads {n_reads} to {n_reads_file}")
        return batch.count_all(**kwargs)

    def _filter_positions(self, info: pd.Series, muts: pd.Series):
        """Mask the positions that do not pass the filters."""
        # Mask the positions with insufficient informative reads.
        if not 1 <= self.min_ninfo_pos:
            raise ValueError(f"min_ninfo_pos must be ≥ 1, but got {self.min_ninfo_pos}")
        self.region.add_mask(
            self.MASK_POS_NINFO, index_to_pos(info.index[info < self.min_ninfo_pos])
        )
        # Mask the positions with excessive mutation fractions.
        if not 0.0 <= self.max_fmut_pos <= 1.0:
            raise ValueError(
                f"max_fmut_pos must be ≥ 0 and ≤ 1, but got {self.max_fmut_pos}"
            )
        self.region.add_mask(
            self.MASK_POS_FMUT,
            index_to_pos(info.index[(muts / info) > self.max_fmut_pos]),
        )

    def _filter_iteration(self):
        """Run an iteration of filtering."""
        # Drop reads that fail to pass the filters.
        tabulator = FilterBatchTabulator(
            get_batch_count_all=self._filter_batch_reads,
            num_batches=self.dataset.num_batches,
            top=self.top,
            sample=self.dataset.sample,
            branches=self.branches,
            refseq=self.dataset.refseq,
            region=self.region,
            pattern=self.pattern,
            min_mut_gap=self.min_mut_gap,
            mut_collisions=self.mut_collisions,
            quick_unbias=self.quick_unbias,
            quick_unbias_thresh=self.quick_unbias_thresh,
            count_pos=True,
            count_read=self.count_read,
            validate=False,
            num_cpus=self.num_cpus,
        )
        # Count the informative and mutated bases.
        info = tabulator.data_per_pos[INFOR_REL]
        muts = tabulator.data_per_pos[MUTAT_REL]
        # Load the checksums and the number of reads.
        self._n_reads[self.DROP_READ_KEPT] = 0
        for batch_num in range(self.dataset.num_batches):
            # Checksums
            checksum_file = self._get_checksum_path(batch_num)
            with open(checksum_file) as f:
                self.checksums[batch_num] = f.read()
            checksum_file.unlink()
            # Number of reads
            n_reads_file = self._get_n_reads_path(batch_num)
            with open(n_reads_file) as f:
                for category, n_reads in json.load(f).items():
                    if self._iter == 1 or category != self.DROP_READ_INIT:
                        self._n_reads[category] += n_reads
                        logger.trace(
                            f"{self} batch {batch_num} had {n_reads} "
                            f"{repr(category)} on iteration {self._iter}"
                        )
            n_reads_file.unlink()
        logger.trace(
            f"{self} on iteration {self._iter} counted "
            + "\n".join(
                f"{category}: {n_reads}" for category, n_reads in self._n_reads.items()
            )
        )
        if self.n_reads_kept == 0:
            logger.warning(f"No reads remained for {self}")
        # Filter out positions based on the new reads.
        self._filter_positions(info, muts)
        if self.pos_kept.size == 0:
            logger.warning(f"No positions remained for {self}")
        return tabulator

    def run_filtering(self):
        if self._iter > 0:
            raise ValueError(f"{self} already filtered the data")
        # Mask positions based on the parameters.
        self._mask_predefined_positions()
        unmasked_curr = self.pos_kept
        self._iter = 1
        while True:
            with logger.debug.begin(f"{self} iteration {self._iter}"):
                unmasked_prev = unmasked_curr
                tabulator = self._filter_iteration()
                unmasked_curr = self.pos_kept
                logger.trace(f"{self} kept {unmasked_curr.size} position(s)")
            # Filtering has converged if the same positions were
            # masked before and after this iteration.
            if np.array_equal(unmasked_prev, unmasked_curr):
                self._converged = True
                logger.debug(f"{self} converged on iteration {self._iter}")
            # Create and save the report after the opportunity to set
            # self._converged to True (so that the report will have the
            # correct value of self._converged) and before returning.
            report = self.create_report()
            report_saved = report.save(self.top, force=self._force_write)
            if self._converged or self._iter >= self.max_iter > 0:
                if not self._converged:
                    logger.warning(
                        f"{self} did not converge within {self.max_iter} iteration(s)"
                    )
                return tabulator, report_saved
            # The first iteration uses the dataset from the IDmut step.
            # Each subsequent iteration uses the nascent Filter dataset
            # in order to resume where the previous iteration ended.
            if self._iter == 1:
                # To be able to load, the nascent Filter dataset must have
                # access to the original IDmut dataset.
                self.dataset.link_data_dirs_to_tmp(self.top)
            self.dataset = FilterMutsDataset(report_saved)
            self._iter += 1

    def create_report(self):
        return FilterReport(
            sample=self.dataset.sample,
            branches=self.branches,
            ref=self.dataset.ref,
            reg=self.region.name,
            end5=self.region.end5,
            end3=self.region.end3,
            checksums={self.CHECKSUM_KEY: self.checksums},
            n_batches=self.dataset.num_batches,
            count_refs=self.pattern.nos,
            count_muts=self.pattern.yes,
            probe=self.probe,
            mask_a=self.mask_a,
            mask_c=self.mask_c,
            mask_g=self.mask_g,
            mask_u=self.mask_u,
            mask_polya=self.mask_polya,
            mask_pos=self.mask_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            max_filter_iter=self.max_iter,
            n_pos_init=self.region.length,
            n_pos_a=self.pos_a.size,
            n_pos_c=self.pos_c.size,
            n_pos_g=self.pos_g.size,
            n_pos_u=self.pos_u.size,
            n_pos_n=self.pos_n.size,
            n_pos_polya=self.pos_polya.size,
            n_pos_list=self.pos_list.size,
            n_pos_min_ninfo=self.pos_min_ninfo.size,
            n_pos_max_fmut=self.pos_max_fmut.size,
            n_pos_kept=self.pos_kept.size,
            pos_a=self.pos_a,
            pos_c=self.pos_c,
            pos_g=self.pos_g,
            pos_u=self.pos_u,
            pos_n=self.pos_n,
            pos_polya=self.pos_polya,
            pos_list=self.pos_list,
            pos_min_ninfo=self.pos_min_ninfo,
            pos_max_fmut=self.pos_max_fmut,
            pos_kept=self.pos_kept,
            min_ncov_read=self.min_ncov_read,
            min_fcov_read=self.min_fcov_read,
            min_finfo_read=self.min_finfo_read,
            max_fmut_read=self.max_fmut_read,
            min_mut_gap=self.min_mut_gap,
            mut_collisions=self.mut_collisions,
            drop_discontig=self.drop_discontig,
            n_reads_init=self.n_reads_init,
            n_reads_list=self.n_reads_list,
            n_reads_min_ncov=self.n_reads_min_ncov,
            n_reads_min_fcov=self.n_reads_min_fcov,
            n_reads_discontig=self.n_reads_discontig,
            n_reads_min_finfo=self.n_reads_min_finfo,
            n_reads_max_fmut=self.n_reads_max_fmut,
            n_reads_min_gap=self.n_reads_min_gap,
            n_reads_kept=self.n_reads_kept,
            quick_unbias=self.quick_unbias,
            quick_unbias_thresh=self.quick_unbias_thresh,
            n_filter_iter=(
                self._iter if self._converged else filter_iter_no_convergence
            ),
            began=self._began,
            ended=datetime.now(),
        )

    def __str__(self):
        return f"{type(self).__name__}: {self.dataset} over {self.region}"


@with_tmp_dir(pass_keep_tmp=False)
def filter_region(
    dataset: IDmutMutsDataset | PoolMutsDataset,
    region: Region,
    *,
    branch: str,
    tmp_dir: Path,
    count_del: bool,
    count_ins: bool,
    no_mut: Iterable[str],
    only_mut: Iterable[str],
    filter_pos_table: bool,
    filter_read_table: bool,
    self_contained: bool,
    force: bool,
    num_cpus: int,
    **kwargs,
):
    """Filter reads, positions, and relationships in a dataset."""
    # Check if the report file already exists.
    branches = path.add_branch(path.FILTER_STEP, branch, dataset.branches)
    report_file = FilterReport.build_path(
        {
            path.TOP: dataset.top,
            path.SAMPLE: dataset.sample,
            path.BRANCHES: branches,
            path.REF: dataset.ref,
            path.REG: region.name,
        }
    )
    if need_write(report_file, force):
        pattern = RelPattern.from_counts(count_del, count_ins, no_mut, only_mut)
        filterer = Filterer(
            dataset,
            region,
            pattern,
            top=tmp_dir,
            branch=branch,
            count_read=filter_read_table,
            self_contained=self_contained,
            num_cpus=num_cpus,
            **kwargs,
        )
        tabulator, report_saved = filterer.run_filtering()
        tabulator.write_tables(pos=filter_pos_table, read=filter_read_table)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
    return report_file.parent
