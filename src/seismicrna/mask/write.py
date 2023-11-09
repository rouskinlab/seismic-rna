"""
Mask -- Write Module
"""

from collections import defaultdict
from datetime import datetime
from logging import getLogger
from typing import Iterable

import numpy as np
import pandas as pd

from .batch import apply_mask
from .io import MaskBatchIO
from .report import MaskReport
from ..core.batch import RefseqMutsBatch, accum_per_pos
from ..core.io import DEFAULT_BROTLI_LEVEL
from ..core.rel import RelPattern
from ..core.seq import Section, index_to_pos
from ..core.write import need_write
from ..relate.data import RelateLoader

logger = getLogger(__name__)


class RelMasker(object):
    """ Mask batches of relation vectors. """

    PATTERN_KEY = "pattern"
    MASK_READ_INIT = "read-init"
    MASK_READ_FINFO = "read-finfo"
    MASK_READ_FMUT = "read-fmut"
    MASK_READ_GAP = "read-gap"
    MASK_READ_KEPT = "read-kept"
    MASK_POS_NINFO = "pos-ninfo"
    MASK_POS_FMUT = "pos-fmut"
    CHECKSUM_KEY = MaskReport.get_batch_type().btype()

    def __init__(self,
                 dataset: RelateLoader,
                 section: Section,
                 pattern: RelPattern, *,
                 exclude_polya: int = 0,
                 exclude_gu: bool = False,
                 exclude_pos: Iterable[tuple[str, int]] = (),
                 min_mut_gap: int = 0,
                 min_finfo_read: float = 0.,
                 max_fmut_read: float = 1.,
                 min_ninfo_pos: int = 0,
                 max_fmut_pos: float = 1.,
                 brotli_level: int = DEFAULT_BROTLI_LEVEL):
        """
        Parameters
        ----------
        dataset: RelateLoader
            Relation vector loader
        section: Section
            The section over which to mask
        pattern: RelPattern
            Relationship pattern
        exclude_polya: int = 0
            Exclude stretches of consecutive A bases at least this long.
            If 0, exclude no bases. Must be ≥ 0.
        exclude_gu: bool = False
            Whether to exclude G and U bases.
        exclude_pos: Iterable[tuple[str, int]] = ()
            Additional, arbitrary positions to exclude. Each position
            must be a tuple of (reference name, 1-indexed position).
        min_mut_gap: int = 0
            Filter out reads with any two mutations separated by fewer
            than `min_mut_gap` positions. Adjacent mutations have a
            gap of 0. If 0, keep all. Must be in [0, length_of_section).
        min_finfo_read: float = 0.0
            Filter out reads with less than this fraction of informative
            bases (i.e. match or mutation). If 0.0, keep all. Must be in
            [0, 1].
        max_fmut_read: float = 1.0
            Filter out reads with more than this fraction of mutated
            bases. If 1.0, keep all. Must be in [0, 1].
        min_ninfo_pos: int = 0
            Filter out positions with less than this number of informative
            bases. Must be ≥ 0.
        max_fmut_pos: float = 1.0
            Filter out positions with more than this fraction of mutated
            reads. Must be in [0, 1].
        """
        # Set the general parameters.
        self.dataset = dataset
        self.section = Section(dataset.ref,
                               dataset.refseq,
                               end5=section.end5,
                               end3=section.end3,
                               name=section.name)
        self.pattern = pattern
        # Set the parameters for excluding positions from the section.
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = np.array([pos for ref, pos in exclude_pos
                                     if ref == self.dataset.ref],
                                    dtype=int)
        # Set the parameters for filtering reads.
        self.min_mut_gap = min_mut_gap
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        self._n_reads = defaultdict(int)
        # Set the parameters for filtering positions.
        if not min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {min_ninfo_pos}")
        self.min_ninfo_pos = min_ninfo_pos
        if not 0. <= max_fmut_pos <= 1.:
            raise ValueError(
                f"max_fmut_pos must be in [0, 1], but got {max_fmut_pos}")
        self.max_fmut_pos = max_fmut_pos
        # Set the parameters for saving files.
        self.top = dataset.top
        self.brotli_level = brotli_level
        self.checksums = list()

    @property
    def n_reads_init(self):
        return self._n_reads[self.MASK_READ_INIT]

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
                if self.exclude_gu
                else np.ndarray([], dtype=int))

    @property
    def pos_polya(self):
        """ Positions masked for lying in a poly(A) sequence. """
        return self.section.get_mask(self.section.MASK_POLYA)

    @property
    def pos_user(self):
        """ Positions masked arbitrarily by the user. """
        return self.section.get_mask(self.section.MASK_POS)

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

    def _filter_min_finfo_read(self, batch: RefseqMutsBatch):
        """ Filter out reads with too few informative positions. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be in [0, 1], but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # All reads have sufficiently many informative positions.
            logger.debug(f"{self} skipped filtering reads with insufficient "
                         f"informative fractions in {batch}")
            return batch
        # Find the reads with sufficiently many informative positions.
        info, muts = batch.count_per_read(self.pattern)
        with np.errstate(invalid="ignore"):
            finfo_read = info.values / self.pos_kept.size
        reads = info.index[finfo_read >= self.min_finfo_read]
        logger.debug(f"{self} kept {reads.size} reads with informative "
                     f"fractions ≥ {self.min_finfo_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _filter_max_fmut_read(self, batch: RefseqMutsBatch):
        """ Filter out reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be in [0, 1], but got "
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
        reads = info.index[fmut_read <= self.max_fmut_read]
        logger.debug(f"{self} kept {reads.size} reads with mutated "
                     f"fractions ≤ {self.max_fmut_read} in {batch}")
        # Return a new batch of only those reads.
        return apply_mask(batch, reads)

    def _mask_min_mut_gap(self, batch: RefseqMutsBatch):
        """ Filter out reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_mut_gap == 0:
            # No read can have a pair of mutations that are too close.
            logger.debug(f"{self} skipped filtering reads with pairs of "
                         f"mutations too close in {batch}")
            return batch
        reads = batch.nonprox_muts(self.pattern, self.min_mut_gap)
        logger.debug(f"{self} kept {reads.size} reads with no two mutations "
                     f"separated by < {self.min_mut_gap} nt in {batch}")
        return apply_mask(batch, reads)

    def _exclude_positions(self):
        """ Exclude positions from the section. """
        self.section.mask_polya(self.exclude_polya)
        if self.exclude_gu:
            self.section.mask_gu()
        self.section.mask_pos(self.exclude_pos)

    def _filter_batch_reads(self, batch: RefseqMutsBatch):
        """ Remove the reads in the batch that do not pass the filters
        and return a new batch without those reads. """
        # Keep only the unmasked positions.
        batch = apply_mask(batch, positions=self.pos_kept)
        self._n_reads[self.MASK_READ_INIT] += (n := batch.num_reads)
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
                                      overwrite=True)
        self.checksums.append(checksum)
        return batch

    def _filter_positions(self, info: pd.Series, muts: pd.Series):
        """ Remove the positions that do not pass the filters. """
        # Mask the positions with insufficient informative reads.
        self.section.add_mask(
            self.MASK_POS_NINFO,
            index_to_pos(info.index[info < self.min_ninfo_pos]))
        # Mask the positions with excessive mutation fractions.
        self.section.add_mask(
            self.MASK_POS_FMUT,
            index_to_pos(info.index[(muts / info) > self.max_fmut_pos]))

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
            exclude_gu=self.exclude_gu,
            exclude_polya=self.exclude_polya,
            exclude_pos=self.exclude_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            n_pos_init=self.section.length,
            n_pos_gu=self.pos_gu.size,
            n_pos_polya=self.pos_polya.size,
            n_pos_user=self.pos_user.size,
            n_pos_min_ninfo=self.pos_min_ninfo.size,
            n_pos_max_fmut=self.pos_max_fmut.size,
            n_pos_kept=self.pos_kept.size,
            pos_gu=self.pos_gu,
            pos_polya=self.pos_polya,
            pos_user=self.pos_user,
            pos_min_ninfo=self.pos_min_ninfo,
            pos_max_fmut=self.pos_max_fmut,
            pos_kept=self.pos_kept,
            min_finfo_read=self.min_finfo_read,
            max_fmut_read=self.max_fmut_read,
            min_mut_gap=self.min_mut_gap,
            n_reads_init=self.n_reads_init,
            n_reads_min_finfo=self.n_reads_min_finfo,
            n_reads_max_fmut=self.n_reads_max_fmut,
            n_reads_min_gap=self.n_reads_min_gap,
            n_reads_kept=self.n_reads_kept,
            began=began,
            ended=ended,
        )

    def __str__(self):
        return f"Mask {self.dataset} over {self.section} with {self.pattern}"


def mask_section(dataset: RelateLoader,
                 section: Section,
                 count_del: bool,
                 count_ins: bool,
                 discount: Iterable[str], *,
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
        pattern = RelPattern.from_counts(count_del, count_ins, discount)
        masker = RelMasker(dataset, section, pattern, **kwargs)
        masker.mask()
        ended = datetime.now()
        report = masker.create_report(began, ended)
        report.save(dataset.top, overwrite=True)
    return report_file

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
