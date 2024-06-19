"""
Mask -- Write Module
"""

from collections import defaultdict
from datetime import datetime
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .batch import apply_mask
from .io import MaskBatchIO
from .report import MaskReport
from ..core.arg import docdef
from ..core.batch import SectionMutsBatch, accum_per_pos
from ..core.rel import RelPattern
from ..core.seq import FIELD_REF, POS_NAME, Section, index_to_pos
from ..core.tmp import release_to_out
from ..core.write import need_write
from ..pool.data import PoolDataset
from ..relate.data import RelateDataset

logger = getLogger(__name__)


class Masker(object):
    """ Mask batches of relation vectors. """

    PATTERN_KEY = "pattern"
    MASK_READ_INIT = "read-init"
    MASK_READ_NCOV = "read-ncov"
    MASK_READ_DISCONTIG = "read-discontig"
    MASK_READ_FINFO = "read-finfo"
    MASK_READ_FMUT = "read-fmut"
    MASK_READ_GAP = "read-gap"
    MASK_READ_KEPT = "read-kept"
    MASK_POS_NINFO = "pos-ninfo"
    MASK_POS_FMUT = "pos-fmut"
    CHECKSUM_KEY = MaskReport.get_batch_type().btype()

    @docdef.auto()
    def __init__(self,
                 dataset: RelateDataset | PoolDataset,
                 section: Section,
                 pattern: RelPattern, *,
                 mask_polya: int,
                 mask_gu: bool,
                 mask_pos: list[tuple[str, int]],
                 mask_pos_file: Path | None,
                 mask_discontig: bool,
                 min_ncov_read: int,
                 min_finfo_read: float,
                 max_fmut_read: float,
                 min_mut_gap: int,
                 min_ninfo_pos: int,
                 max_fmut_pos: float,
                 quick_unbias: bool,
                 quick_unbias_thresh: float,
                 brotli_level: int,
                 top: Path):
        # Set the general parameters.
        self.dataset = dataset
        self.section = Section(dataset.ref,
                               dataset.refseq,
                               end5=section.end5,
                               end3=section.end3,
                               name=section.name)
        self.pattern = pattern
        # Set the parameters for excluding positions from the section.
        if mask_polya <= 0 or mask_polya > 5:
            logger.warning("It is not recommended to keep sequences of 5 or "
                           "more consecutive As because of an artifact during "
                           "RT that causes low reactivity. See Kladwang et al. "
                           "(https://doi.org/10.1021/acs.biochem.0c00020).")
        self.mask_polya = mask_polya
        self.mask_gu = mask_gu
        self.mask_pos = self._get_mask_pos(mask_pos, mask_pos_file)
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
        self.brotli_level = brotli_level
        self.checksums = list()

    @property
    def n_reads_premask(self):
        return self._n_reads[self.MASK_READ_INIT]

    @property
    def n_reads_min_ncov(self):
        return self._n_reads[self.MASK_READ_NCOV]

    @property
    def n_reads_discontig(self):
        return self._n_reads[self.MASK_READ_DISCONTIG]

    @property
    def n_reads_min_finfo(self):
        return self._n_reads[self.MASK_READ_FINFO]

    @property
    def n_reads_max_fmut(self):
        return self._n_reads[self.MASK_READ_FMUT]

    @property
    def n_reads_min_gap(self):
        return self._n_reads[self.MASK_READ_GAP]

    @property
    def n_reads_kept(self):
        """ Number of reads kept. """
        return self._n_reads[self.MASK_READ_KEPT]

    @property
    def pos_gu(self):
        """ Positions masked for having a G or U base. """
        return (self.section.get_mask(self.section.MASK_GU)
                if self.mask_gu
                else np.array([], dtype=int))

    @property
    def pos_polya(self):
        """ Positions masked for lying in a poly(A) sequence. """
        return self.section.get_mask(self.section.MASK_POLYA)

    @property
    def pos_list(self):
        """ Positions masked arbitrarily from a list. """
        return self.section.get_mask(self.section.MASK_LIST)

    @property
    def pos_min_ninfo(self):
        """ Positions masked for having too few informative reads. """
        return self.section.get_mask(self.MASK_POS_NINFO)

    @property
    def pos_max_fmut(self):
        """ Positions masked for having too many mutations. """
        return self.section.get_mask(self.MASK_POS_FMUT)

    @property
    def pos_kept(self):
        """ Positions kept. """
        return self.section.unmasked_int

    @property
    def n_batches(self):
        """ Number of batches of reads. """
        return len(self.checksums)

    def _get_mask_pos(self,
                      mask_pos: list[tuple[str, int]],
                      mask_pos_file: Path | None):
        """ List all positions to mask. """
        # Collect the positions to mask from the list.
        mask_pos = np.array([pos for ref, pos in mask_pos
                             if ref == self.dataset.ref])
        # Read positions to exclude from a file.
        if mask_pos_file is not None:
            mask_pos = np.concatenate([mask_pos,
                                       pd.read_csv(
                                           mask_pos_file,
                                           index_col=[FIELD_REF, POS_NAME]
                                       ).loc[self.dataset.ref].index.values])
        # Drop redundant positions and sort the remaining ones.
        mask_pos = np.unique(np.asarray(mask_pos, dtype=int))
        # Keep only the positions in the section.
        return mask_pos[np.logical_and(mask_pos >= self.section.end5,
                                       mask_pos <= self.section.end3)]

    def _filter_min_ncov_read(self, batch: SectionMutsBatch):
        """ Filter out reads with too few covered positions. """
        if self.min_ncov_read < 1:
            raise ValueError(f"min_ncov_read must be ≥ 1, but got "
                             f"{self.min_ncov_read}")
        # Find the reads with sufficiently many covered positions.
        reads = batch.read_nums[batch.cover_per_read.values.sum(axis=1)
                                >= self.min_ncov_read]
        logger.debug(f"{self} kept {reads.size} reads with coverage "
                     f"≥ {self.min_ncov_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_discontig(self, batch: SectionMutsBatch):
        """ Filter out reads with discontiguous mates. """
        if not self.mask_discontig:
            # Keep discontiguous reads.
            logger.debug(f"{self} skipped filtering reads with "
                         f"discontiguous mates in {batch}")
            return batch
        # Find the reads with contiguous mates.
        reads = batch.read_nums[batch.contiguous]
        logger.debug(f"{self} kept {reads.size} reads with "
                     f"contiguous mates in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_min_finfo_read(self, batch: SectionMutsBatch):
        """ Filter out reads with too few informative positions. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be ≥ 0, ≤ 1, but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # All reads have sufficiently many informative positions.
            logger.debug(f"{self} skipped filtering reads with insufficient "
                         f"informative fractions in {batch}")
            return batch
        # Find the reads with sufficiently many informative positions.
        info, muts = batch.count_per_read(self.pattern)
        finfo_read = info.values / batch.cover_per_read.values.sum(axis=1)
        reads = info.index.values[finfo_read >= self.min_finfo_read]
        logger.debug(f"{self} kept {reads.size} reads with informative "
                     f"fractions ≥ {self.min_finfo_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_max_fmut_read(self, batch: SectionMutsBatch):
        """ Filter out reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be ≥ 0, ≤ 1, but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read == 1.:
            # All reads have sufficiently few mutations.
            logger.debug(f"{self} skipped filtering reads with excessive "
                         f"mutation fractions in {batch}")
            return batch
        # Find the reads with sufficiently few mutations.
        info, muts = batch.count_per_read(self.pattern)
        with np.errstate(invalid="ignore"):
            fmut_read = muts.values / info.values
        reads = info.index.values[fmut_read <= self.max_fmut_read]
        logger.debug(f"{self} kept {reads.size} reads with mutated "
                     f"fractions ≤ {self.max_fmut_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _mask_min_mut_gap(self, batch: SectionMutsBatch):
        """ Filter out reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_mut_gap == 0:
            # No read can have a pair of mutations that are too close.
            logger.debug(f"{self} skipped filtering reads with pairs of "
                         f"mutations too close in {batch}")
            return batch
        reads = batch.reads_noclose_muts(self.pattern, self.min_mut_gap)
        logger.debug(f"{self} kept {reads.size} reads with no two mutations "
                     f"separated by < {self.min_mut_gap} nt in {batch}")
        return apply_mask(batch, reads)

    def _exclude_positions(self):
        """ Exclude positions from the section. """
        self.section.mask_polya(self.mask_polya)
        if self.mask_gu:
            self.section.mask_gu()
        self.section.mask_list(self.mask_pos)

    def _filter_batch_reads(self, batch: SectionMutsBatch):
        """ Remove the reads in the batch that do not pass the filters
        and return a new batch without those reads. """
        # Determine the initial number of reads in the batch.
        self._n_reads[self.MASK_READ_INIT] += (n := batch.num_reads)
        # Keep only the unmasked positions.
        batch = apply_mask(batch, section=self.section)
        # Remove reads with too few covered positions.
        batch = self._filter_min_ncov_read(batch)
        self._n_reads[self.MASK_READ_NCOV] += (n - (n := batch.num_reads))
        # Remove reads with discontiguous mates.
        batch = self._filter_discontig(batch)
        self._n_reads[self.MASK_READ_DISCONTIG] += (n - (n := batch.num_reads))
        # Remove reads with too few informative positions.
        batch = self._filter_min_finfo_read(batch)
        self._n_reads[self.MASK_READ_FINFO] += (n - (n := batch.num_reads))
        # Remove reads with too many mutations.
        batch = self._filter_max_fmut_read(batch)
        self._n_reads[self.MASK_READ_FMUT] += (n - (n := batch.num_reads))
        # Remove reads with mutations too close together.
        batch = self._mask_min_mut_gap(batch)
        self._n_reads[self.MASK_READ_GAP] += (n - (n := batch.num_reads))
        # Record the number of reads remaining after filtering.
        self._n_reads[self.MASK_READ_KEPT] += n
        # Save the batch.
        batch_file = MaskBatchIO(sample=self.dataset.sample,
                                 ref=self.dataset.ref,
                                 sect=self.section.name,
                                 batch=batch.batch,
                                 read_nums=batch.read_nums)
        _, checksum = batch_file.save(self.top,
                                      brotli_level=self.brotli_level,
                                      force=True)
        self.checksums.append(checksum)
        return batch

    def _filter_positions(self, info: pd.Series, muts: pd.Series):
        """ Remove the positions that do not pass the filters. """
        # Mask the positions with insufficient informative reads.
        if not 1 <= self.min_ninfo_pos:
            raise ValueError("min_ninfo_pos must be ≥ 1, "
                             f"but got {self.min_ninfo_pos}")
        self.section.add_mask(
            self.MASK_POS_NINFO,
            index_to_pos(info.index[info < self.min_ninfo_pos])
        )
        # Mask the positions with excessive mutation fractions.
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError("max_fmut_pos must be ≥ 0 and ≤ 1, "
                             f"but got {self.max_fmut_pos}")
        self.section.add_mask(
            self.MASK_POS_FMUT,
            index_to_pos(info.index[(muts / info) > self.max_fmut_pos])
        )

    def mask(self):
        # Exclude positions based on the parameters.
        self._exclude_positions()
        if self.pos_kept.size == 0:
            logger.warning(f"No positions remained after excluding with {self}")
        # Filter out reads based on the parameters and count the number
        # of informative and mutated positions remaining.
        _, muts, info = accum_per_pos(map(self._filter_batch_reads,
                                          self.dataset.iter_batches()),
                                      self.dataset.refseq,
                                      self.pos_kept,
                                      {self.PATTERN_KEY: self.pattern})
        if self.n_reads_kept == 0:
            logger.warning(f"No reads remained after filtering with {self}")
        # Filter out positions based on the parameters.
        self._filter_positions(info[self.PATTERN_KEY], muts[self.PATTERN_KEY])
        if self.pos_kept.size == 0:
            logger.warning(f"No positions remained after filtering with {self}")

    def create_report(self, began: datetime, ended: datetime):
        return MaskReport(
            sample=self.dataset.sample,
            ref=self.dataset.ref,
            sect=self.section.name,
            end5=self.section.end5,
            end3=self.section.end3,
            checksums={self.CHECKSUM_KEY: self.checksums},
            n_batches=self.n_batches,
            count_refs=self.pattern.nos,
            count_muts=self.pattern.yes,
            mask_gu=self.mask_gu,
            mask_polya=self.mask_polya,
            mask_pos=self.mask_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            n_pos_init=self.section.length,
            n_pos_gu=self.pos_gu.size,
            n_pos_polya=self.pos_polya.size,
            n_pos_user=self.pos_list.size,
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
            n_reads_premask=self.n_reads_premask,
            n_reads_min_ncov=self.n_reads_min_ncov,
            n_reads_discontig=self.n_reads_discontig,
            n_reads_min_finfo=self.n_reads_min_finfo,
            n_reads_max_fmut=self.n_reads_max_fmut,
            n_reads_min_gap=self.n_reads_min_gap,
            n_reads_kept=self.n_reads_kept,
            quick_unbias=self.quick_unbias,
            quick_unbias_thresh=self.quick_unbias_thresh,
            began=began,
            ended=ended,
        )

    def __str__(self):
        return f"Mask {self.dataset} over {self.section} with {self.pattern}"


def mask_section(dataset: RelateDataset | PoolDataset,
                 section: Section,
                 mask_del: bool,
                 mask_ins: bool,
                 mask_mut: Iterable[str], *,
                 tmp_dir: Path,
                 force: bool,
                 **kwargs):
    """ Filter a section of a set of bit vectors. """
    # Check if the report file already exists.
    report_file = MaskReport.build_path(top=dataset.top,
                                        sample=dataset.sample,
                                        ref=dataset.ref,
                                        sect=section.name)
    if need_write(report_file, force):
        began = datetime.now()
        pattern = RelPattern.from_counts(not mask_del, not mask_ins, mask_mut)
        masker = Masker(dataset,
                        section,
                        pattern,
                        top=tmp_dir,
                        **kwargs)
        masker.mask()
        ended = datetime.now()
        report = masker.create_report(began, ended)
        report_saved = report.save(tmp_dir)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
    return report_file

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
