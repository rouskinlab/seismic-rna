import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .batch import apply_mask
from .dataset import MaskMutsDataset
from .io import MaskBatchIO
from .report import MaskReport
from .table import MaskBatchTabulator
from ..core import path
from ..core.arg import docdef
from ..core.batch import RegionMutsBatch
from ..core.dataset import MissingBatchTypeError
from ..core.error import IncompatibleValuesError
from ..core.lists import PositionList
from ..core.logs import logger
from ..core.rel import RelPattern, HalfRelPattern
from ..core.report import mask_iter_no_convergence
from ..core.seq import FIELD_REF, POS_NAME, Region, index_to_pos
from ..core.table import MUTAT_REL, INFOR_REL
from ..core.tmp import release_to_out, with_tmp_dir
from ..core.write import need_write
from ..relate.dataset import (RelateMutsDataset,
                              PoolDataset,
                              load_read_names_dataset)


class Masker(object):
    """ Mask batches of relationships. """

    PATTERN_KEY = "pattern"
    MASK_READ_INIT = "read-init"
    MASK_READ_LIST = "read-exclude"
    MASK_READ_NCOV = "read-ncov"
    MASK_READ_DISCONTIG = "read-discontig"
    MASK_READ_FINFO = "read-finfo"
    MASK_READ_FMUT = "read-fmut"
    MASK_READ_GAP = "read-gap"
    MASK_READ_KEPT = "read-kept"
    MASK_POS_NINFO = "pos-ninfo"
    MASK_POS_FMUT = "pos-fmut"
    CHECKSUM_KEY = MaskReport.get_batch_type().btype()

    # This decorator is used to make Masker instances easier to create
    # using the Python API, since the arguments with default values do
    # not need to be specified; but it can also obscure bugs caused by
    # failing to pass an argument to __init__, so it's best to comment
    # out @docdef.auto() when developing the source code.
    @docdef.auto()
    def __init__(self,
                 dataset: RelateMutsDataset | PoolDataset,
                 region: Region,
                 pattern: RelPattern, *,
                 max_mask_iter: int,
                 mask_polya: int,
                 mask_gu: bool,
                 mask_pos: list[tuple[str, int]],
                 mask_pos_file: list[Path],
                 mask_read: list[str],
                 mask_read_file: list[Path],
                 mask_discontig: bool,
                 min_ncov_read: int,
                 min_finfo_read: float,
                 max_fmut_read: float,
                 min_mut_gap: int,
                 min_ninfo_pos: int,
                 max_fmut_pos: float,
                 quick_unbias: bool,
                 quick_unbias_thresh: float,
                 count_read: bool,
                 brotli_level: int,
                 top: Path,
                 branch: str,
                 max_procs: int = 1):
        # Set the general parameters.
        self._began = datetime.now()
        self.dataset = dataset
        self.region = Region(dataset.ref,
                             dataset.refseq,
                             end5=region.end5,
                             end3=region.end3,
                             name=region.name)
        self.pattern = pattern
        self.max_iter = max_mask_iter
        self._iter = 0
        self._converged = False
        # Set the parameters for excluding positions from the region.
        if not 0 < mask_polya <= 5:
            logger.warning("It is not recommended to keep sequences of 5 or "
                           "more consecutive As because of an artifact during "
                           "RT that causes low reactivity. See Kladwang et al. "
                           "(https://doi.org/10.1021/acs.biochem.0c00020).")
        self.mask_polya = mask_polya
        self.mask_gu = mask_gu
        self.mask_pos = self._get_mask_pos(mask_pos, mask_pos_file)
        self.mask_read = self._get_mask_read(mask_read, mask_read_file)
        # Set the parameters for filtering reads.
        if min_mut_gap > 0 and not mask_discontig:
            raise ValueError("The observer bias correction does not work with "
                             "discontiguous reads. If you need discontiguous "
                             "reads, disable bias correction with the option "
                             "--min-mut-gap=0 (but be warned that disabling "
                             "bias correction can produce misleading results, "
                             "especially with clustering).")
        self.mask_discontig = mask_discontig
        self.min_ncov_read = min_ncov_read
        self.min_mut_gap = min_mut_gap
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
        self.checksums = [""] * dataset.num_batches
        # After the first iteration, self.dataset will become the new,
        # masked dataset, which will also have a branch for the mask
        # step, so calculate the branches using the original dataset.
        self.branches = path.add_branch(path.MASK_STEP,
                                        branch,
                                        dataset.branches)
        # Parallelization
        self.max_procs = max_procs

    # This property can change: do not cache it.
    @property
    def n_reads_init(self):
        return self._n_reads[self.MASK_READ_INIT]

    # This property can change: do not cache it.
    @property
    def n_reads_list(self):
        return self._n_reads[self.MASK_READ_LIST]

    # This property can change: do not cache it.
    @property
    def n_reads_min_ncov(self):
        return self._n_reads[self.MASK_READ_NCOV]

    # This property can change: do not cache it.
    @property
    def n_reads_discontig(self):
        return self._n_reads[self.MASK_READ_DISCONTIG]

    # This property can change: do not cache it.
    @property
    def n_reads_min_finfo(self):
        return self._n_reads[self.MASK_READ_FINFO]

    # This property can change: do not cache it.
    @property
    def n_reads_max_fmut(self):
        return self._n_reads[self.MASK_READ_FMUT]

    # This property can change: do not cache it.
    @property
    def n_reads_min_gap(self):
        return self._n_reads[self.MASK_READ_GAP]

    # This property can change: do not cache it.
    @property
    def n_reads_kept(self):
        """ Number of reads kept. """
        return self._n_reads[self.MASK_READ_KEPT]

    # This property can change: do not cache it.
    @property
    def pos_gu(self):
        """ Positions masked for having a G or U base. """
        return (self.region.get_mask(self.region.MASK_GU)
                if self.mask_gu
                else np.array([], dtype=int))

    # This property can change: do not cache it.
    @property
    def pos_polya(self):
        """ Positions masked for lying in a poly(A) sequence. """
        return self.region.get_mask(self.region.MASK_POLYA)

    # This property can change: do not cache it.
    @property
    def pos_list(self):
        """ Positions masked arbitrarily from a list. """
        return self.region.get_mask(self.region.MASK_LIST)

    # This property can change: do not cache it.
    @property
    def pos_min_ninfo(self):
        """ Positions masked for having too few informative reads. """
        return self.region.get_mask(self.MASK_POS_NINFO)

    # This property can change: do not cache it.
    @property
    def pos_max_fmut(self):
        """ Positions masked for having too many mutations. """
        return self.region.get_mask(self.MASK_POS_FMUT)

    # This property can change: do not cache it.
    @property
    def pos_kept(self):
        """ Positions kept. """
        return self.region.unmasked_int

    # This property can change: do not cache it.
    @property
    def _force_write(self):
        """ Whether to force-write each file. """
        return self._iter > 1

    # This property can change: do not cache it.
    @property
    def read_names_dataset(self):
        """ Dataset of the read names. """
        return load_read_names_dataset(self.dataset.report_file)

    def _get_mask_pos(self,
                      mask_pos: Iterable[tuple[str, int]],
                      mask_pos_file: Iterable[str | Path]):
        """ List all positions to mask. """
        # Collect the positions to mask from the list.
        dataset_ref = self.dataset.ref
        mask_pos = np.array([pos for ref, pos in mask_pos
                             if ref == dataset_ref],
                            dtype=int)
        logger.detail(f"Got {mask_pos.size} positions listed individually "
                      f"to pre-exclude for reference {repr(dataset_ref)}")
        # List positions to exclude from file(s).
        for file in path.find_files_chain(mask_pos_file,
                                          [path.PositionListSeg]):
            file_data = PositionList.load_data(file)
            ref_rows = file_data[FIELD_REF] == dataset_ref
            file_pos = file_data.loc[ref_rows, POS_NAME].values
            logger.detail(f"Got {file_pos.size} positions in {file} "
                          f"to pre-exclude for reference {repr(dataset_ref)}")
            if file_pos.size > 0:
                mask_pos = np.concatenate([mask_pos, file_pos])
        # Drop redundant positions and sort the remaining ones.
        mask_pos = np.unique(np.asarray(mask_pos, dtype=int))
        # Keep only the positions in the region.
        mask_pos = mask_pos[np.logical_and(mask_pos >= self.region.end5,
                                           mask_pos <= self.region.end3)]
        logger.detail(f"Got {mask_pos.size} positions to pre-exclude "
                      f"for reference {repr(dataset_ref)}")
        return mask_pos

    @staticmethod
    def _get_mask_read(mask_read: Iterable[str],
                       mask_read_file: Iterable[str | Path]):
        """ List all reads to pre-exclude. """
        # Ensure that the given read names are all strings.
        mask_read = set(map(str, mask_read))
        # List reads to exclude from file(s).
        for file in path.find_files_chain(mask_read_file,
                                          [path.ReadListSeg]):
            with open(file) as f:
                # Ensure every read name is unique.
                mask_read.update(map(str.rstrip, f))
        mask_read = np.asarray(list(mask_read), dtype=str)
        logger.detail(f"Got {mask_read.size} reads to pre-exclude")
        return mask_read

    def _filter_exclude_read(self, batch: RegionMutsBatch):
        """ Filter out reads in the list of pre-excluded read. """
        if self.mask_read.size == 0:
            # Pre-exclude no reads.
            logger.detail(f"{self} skipped pre-excluding reads in {batch}")
            return batch
        # Load the names of the reads in this batch.
        try:
            names_batch = self.read_names_dataset.get_batch(batch.batch)
        except MissingBatchTypeError:
            raise IncompatibleValuesError(
                "Reads can be excluded with --mask-read or --mask-read-file "
                "only if relate was run using the option --write-read-names"
            ) from None
        if names_batch.num_reads != batch.num_reads:
            raise ValueError(f"Expected {batch.num_reads} read names,"
                             f"but got {names_batch.num_reads}")
        # Find the numbers of the reads to keep.
        reads = batch.read_nums[np.isin(names_batch.names,
                                        self.mask_read,
                                        assume_unique=True,
                                        invert=True)]
        logger.detail(f"{self} kept {reads.size} reads after pre-excluding")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_min_ncov_read(self, batch: RegionMutsBatch):
        """ Filter out reads with too few covered positions. """
        if self.min_ncov_read < 1:
            raise ValueError(f"min_ncov_read must be ≥ 1, but got "
                             f"{self.min_ncov_read}")
        # Find the reads with sufficiently many covered positions.
        reads = batch.read_nums[batch.cover_per_read.values.sum(axis=1)
                                >= self.min_ncov_read]
        logger.detail(f"{self} kept {reads.size} reads with coverage "
                      f"≥ {self.min_ncov_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_discontig(self, batch: RegionMutsBatch):
        """ Filter out reads with discontiguous mates. """
        if not self.mask_discontig:
            # Keep discontiguous reads.
            logger.detail(f"{self} skipped filtering reads with "
                          f"discontiguous mates in {batch}")
            return batch
        # Find the reads with contiguous mates.
        reads = batch.read_nums[batch.contiguous]
        logger.detail(f"{self} kept {reads.size} reads with "
                      f"contiguous mates in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_min_finfo_read(self, batch: RegionMutsBatch):
        """ Filter out reads with too few informative positions. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be ≥ 0, ≤ 1, but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # All reads have sufficiently many informative positions.
            logger.detail(f"{self} skipped filtering reads with insufficient "
                          f"informative fractions in {batch}")
            return batch
        # Find the reads with sufficiently many informative positions.
        info, muts = batch.count_per_read(self.pattern)
        finfo_read = info.values / batch.cover_per_read.values.sum(axis=1)
        reads = info.index.values[finfo_read >= self.min_finfo_read]
        logger.detail(f"{self} kept {reads.size} reads with informative "
                      f"fractions ≥ {self.min_finfo_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_max_fmut_read(self, batch: RegionMutsBatch):
        """ Filter out reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be ≥ 0, ≤ 1, but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read == 1.:
            # All reads have sufficiently few mutations.
            logger.detail(f"{self} skipped filtering reads with excessive "
                          f"mutation fractions in {batch}")
            return batch
        # Find the reads with sufficiently few mutations.
        info, muts = batch.count_per_read(self.pattern)
        with np.errstate(invalid="ignore"):
            fmut_read = muts.values / info.values
        reads = info.index.values[fmut_read <= self.max_fmut_read]
        logger.detail(f"{self} kept {reads.size} reads with mutated "
                      f"fractions ≤ {self.max_fmut_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _mask_min_mut_gap(self, batch: RegionMutsBatch):
        """ Filter out reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}"
            )
        if self.min_mut_gap == 0:
            # No read can have a pair of mutations that are too close.
            logger.detail(f"{self} skipped filtering reads with pairs of "
                          f"mutations too close in {batch}")
            return batch
        reads = batch.reads_noclose_muts(self.pattern, self.min_mut_gap)
        logger.detail(f"{self} kept {reads.size} reads with no two mutations "
                      f"separated by < {self.min_mut_gap} nt in {batch}")
        return apply_mask(batch, reads)

    def _exclude_positions(self):
        """ Pre-exclude positions from the region. """
        self.region.mask_polya(self.mask_polya)
        if self.mask_gu:
            self.region.mask_gu()
        self.region.mask_list(self.mask_pos)

    def _get_batch_num_path(self, batch_num: int):
        return MaskBatchIO.build_path({path.TOP: self.top,
                                       path.SAMPLE: self.dataset.sample,
                                       path.BRANCHES: self.branches,
                                       path.REF: self.dataset.ref,
                                       path.REG: self.region.name,
                                       path.BATCH: batch_num})

    def _get_n_reads_path(self, batch_num: int):
        """ Get the path to the file of the number of reads masked out
        for a batch. """
        return self._get_batch_num_path(batch_num).with_suffix(path.JSON_EXT)

    def _get_checksum_path(self, batch_num: int):
        """ Get the path to the file of the checksum for a batch. """
        return self._get_batch_num_path(batch_num).with_suffix(path.TXT_EXT)

    def _filter_batch_reads(self, batch_num: int, **kwargs):
        """ Remove the reads in the batch that do not pass the filters
        and return a new batch without those reads. """
        n_reads = dict()
        # Load the batch.
        batch = self.dataset.get_batch(batch_num)
        # Determine the initial number of reads in the batch.
        n_reads[self.MASK_READ_INIT] = (n := batch.num_reads)
        if self._iter == 1:
            # Keep only the unmasked positions.
            batch = apply_mask(batch, region=self.region)
            # Remove pre-excluded reads.
            batch = self._filter_exclude_read(batch)
            n_reads[self.MASK_READ_LIST] = (n - (n := batch.num_reads))
        # Remove reads with too few covered positions.
        batch = self._filter_min_ncov_read(batch)
        n_reads[self.MASK_READ_NCOV] = (n - (n := batch.num_reads))
        # Remove reads with discontiguous mates.
        batch = self._filter_discontig(batch)
        n_reads[self.MASK_READ_DISCONTIG] = (n - (n := batch.num_reads))
        # Remove reads with too few informative positions.
        batch = self._filter_min_finfo_read(batch)
        n_reads[self.MASK_READ_FINFO] = (n - (n := batch.num_reads))
        # Remove reads with too many mutations.
        batch = self._filter_max_fmut_read(batch)
        n_reads[self.MASK_READ_FMUT] = (n - (n := batch.num_reads))
        # Remove reads with mutations too close together.
        batch = self._mask_min_mut_gap(batch)
        n_reads[self.MASK_READ_GAP] = (n - (n := batch.num_reads))
        # Record the number of reads remaining after filtering.
        n_reads[self.MASK_READ_KEPT] = n
        # Save the batch.
        batch_file = MaskBatchIO(sample=self.dataset.sample,
                                 branches=self.branches,
                                 ref=self.dataset.ref,
                                 reg=self.region.name,
                                 batch=batch.batch,
                                 read_nums=batch.read_nums)
        _, checksum = batch_file.save(self.top,
                                      brotli_level=self.brotli_level,
                                      force=self._force_write)
        # Save the checksum.
        checksum_file = self._get_checksum_path(batch_num)
        with open(checksum_file, "x") as f:
            f.write(checksum)
            logger.detail(f"Wrote checksum {repr(checksum)} to {checksum_file}")
        # Save the read counts.
        n_reads_file = self._get_n_reads_path(batch_num)
        with open(n_reads_file, "x") as f:
            json.dump(n_reads, f)
            logger.detail(f"Wrote number of reads {n_reads} to {n_reads_file}")
        return batch.count_all(**kwargs)

    def _filter_positions(self, info: pd.Series, muts: pd.Series):
        """ Remove the positions that do not pass the filters. """
        # Mask the positions with insufficient informative reads.
        if not 1 <= self.min_ninfo_pos:
            raise ValueError("min_ninfo_pos must be ≥ 1, "
                             f"but got {self.min_ninfo_pos}")
        self.region.add_mask(
            self.MASK_POS_NINFO,
            index_to_pos(info.index[info < self.min_ninfo_pos])
        )
        # Mask the positions with excessive mutation fractions.
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError("max_fmut_pos must be ≥ 0 and ≤ 1, "
                             f"but got {self.max_fmut_pos}")
        self.region.add_mask(
            self.MASK_POS_FMUT,
            index_to_pos(info.index[(muts / info) > self.max_fmut_pos])
        )

    def _mask_iteration(self):
        """ Run an iteration of masking. """
        # Mask out reads that fail to pass the filter parameters.
        tabulator = MaskBatchTabulator(
            get_batch_count_all=self._filter_batch_reads,
            num_batches=self.dataset.num_batches,
            top=self.top,
            sample=self.dataset.sample,
            branches=self.branches,
            refseq=self.dataset.refseq,
            region=self.region,
            pattern=self.pattern,
            min_mut_gap=self.min_mut_gap,
            quick_unbias=self.quick_unbias,
            quick_unbias_thresh=self.quick_unbias_thresh,
            count_pos=True,
            count_read=self.count_read,
            validate=False,
            max_procs=self.max_procs,
        )
        # Count the informative and mutated bases.
        info = tabulator.data_per_pos[INFOR_REL]
        muts = tabulator.data_per_pos[MUTAT_REL]
        # Load the checksums and the number of reads.
        self._n_reads[self.MASK_READ_KEPT] = 0
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
                    if self._iter == 1 or category != self.MASK_READ_INIT:
                        self._n_reads[category] += n_reads
                        logger.detail(
                            f"{self} batch {batch_num} had {n_reads} "
                            f"{repr(category)} on iteration {self._iter}"
                        )
            n_reads_file.unlink()
        logger.detail(f"{self} on iteration {self._iter} counted "
                      + "\n".join(f"{category}: {n_reads}"
                                  for category, n_reads
                                  in self._n_reads.items()))
        if self.n_reads_kept == 0:
            logger.warning(f"No reads remained for {self}")
        # Filter out positions based on the parameters.
        self._filter_positions(info, muts)
        if self.pos_kept.size == 0:
            logger.warning(f"No positions remained for {self}")
        return tabulator

    def mask(self):
        if self._iter > 0:
            raise ValueError(f"{self} already masked the data")
        # Exclude positions based on the parameters.
        self._exclude_positions()
        unmasked_curr = self.pos_kept
        self._iter = 1
        while True:
            logger.routine(f"Began {self} iteration {self._iter}")
            unmasked_prev = unmasked_curr
            tabulator = self._mask_iteration()
            unmasked_curr = self.pos_kept
            logger.detail(f"{self} kept {unmasked_curr.size} position(s)")
            logger.routine(f"Ended {self} iteration {self._iter}")
            # Masking has converged if either the same positions were
            # masked before and after this iteration.
            if np.array_equal(unmasked_prev, unmasked_curr):
                self._converged = True
                logger.routine(f"{self} converged on iteration {self._iter}")
            # Create and save the report after the opportunity to set
            # self._converged to True (so that the report will have the
            # correct value of self._converged) and before returning.
            report = self.create_report()
            report_saved = report.save(self.top, force=self._force_write)
            if self._converged or self._iter >= self.max_iter > 0:
                if not self._converged:
                    logger.warning(f"{self} did not converge "
                                   f"within {self.max_iter} iteration(s)")
                return tabulator, report_saved
            # The first iteration uses the dataset from the relate step.
            # Each subsequent iteration uses the nascent mask dataset
            # in order to resume where the previous iteration ended.
            if self._iter == 1:
                # To be able to load, the nascent mask dataset must have
                # access to the original relate dataset.
                self.dataset.link_data_dirs_to_tmp(self.top)
            self.dataset = MaskMutsDataset(report_saved)
            self._iter += 1

    def create_report(self):
        return MaskReport(
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
            mask_gu=self.mask_gu,
            mask_polya=self.mask_polya,
            mask_pos=self.mask_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            max_mask_iter=self.max_iter,
            n_pos_init=self.region.length,
            n_pos_gu=self.pos_gu.size,
            n_pos_polya=self.pos_polya.size,
            n_pos_list=self.pos_list.size,
            n_pos_min_ninfo=self.pos_min_ninfo.size,
            n_pos_max_fmut=self.pos_max_fmut.size,
            n_pos_kept=self.pos_kept.size,
            pos_gu=self.pos_gu,
            pos_polya=self.pos_polya,
            pos_list=self.pos_list,
            pos_min_ninfo=self.pos_min_ninfo,
            pos_max_fmut=self.pos_max_fmut,
            pos_kept=self.pos_kept,
            min_ncov_read=self.min_ncov_read,
            min_finfo_read=self.min_finfo_read,
            max_fmut_read=self.max_fmut_read,
            min_mut_gap=self.min_mut_gap,
            mask_discontig=self.mask_discontig,
            n_reads_init=self.n_reads_init,
            n_reads_list=self.n_reads_list,
            n_reads_min_ncov=self.n_reads_min_ncov,
            n_reads_discontig=self.n_reads_discontig,
            n_reads_min_finfo=self.n_reads_min_finfo,
            n_reads_max_fmut=self.n_reads_max_fmut,
            n_reads_min_gap=self.n_reads_min_gap,
            n_reads_kept=self.n_reads_kept,
            quick_unbias=self.quick_unbias,
            quick_unbias_thresh=self.quick_unbias_thresh,
            n_mask_iter=(self._iter
                         if self._converged
                         else mask_iter_no_convergence),
            began=self._began,
            ended=datetime.now(),
        )

    def __str__(self):
        return f"Mask {self.dataset} over {self.region}"


@with_tmp_dir(pass_keep_tmp=False)
def mask_region(dataset: RelateMutsDataset | PoolDataset,
                region: Region, *,
                branch: str,
                tmp_dir: Path,
                mask_del: bool,
                mask_ins: bool,
                mask_mut: Iterable[str],
                count_mut: Iterable[str],
                mask_pos_table: bool,
                mask_read_table: bool,
                force: bool,
                n_procs: int,
                **kwargs):
    """ Mask out certain reads, positions, and relationships. """
    # Check if the report file already exists.
    branches = path.add_branch(path.MASK_STEP, branch, dataset.branches)
    report_file = MaskReport.build_path({path.TOP: dataset.top,
                                         path.SAMPLE: dataset.sample,
                                         path.BRANCHES: branches,
                                         path.REF: dataset.ref,
                                         path.REG: region.name})
    if need_write(report_file, force):
        if count_mut:
            pattern = RelPattern(HalfRelPattern.from_counts(count_sub=False,
                                                            count_del=not mask_del,
                                                            count_ins=not mask_ins,
                                                            count=count_mut,
                                                            discount=mask_mut),
                                 HalfRelPattern.from_counts(count_ref=True,
                                                            discount=mask_mut))
        else:
            pattern = RelPattern.from_counts(not mask_del, not mask_ins, mask_mut)
        masker = Masker(dataset,
                        region,
                        pattern,
                        top=tmp_dir,
                        branch=branch,
                        count_read=mask_read_table,
                        max_procs=n_procs,
                        **kwargs)
        tabulator, report_saved = masker.mask()
        tabulator.write_tables(pos=mask_pos_table, read=mask_read_table)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
    return report_file.parent
