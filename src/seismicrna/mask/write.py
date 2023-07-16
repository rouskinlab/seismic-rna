"""
Mask -- Write Module
"""

from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .report import MaskReport
from ..core import path
from ..core.bitcall import BitCaller
from ..core.bitvect import BitBatch, BitCounter
from ..core.files import digest_file
from ..core.sect import Section, index_to_pos
from ..relate.load import RelateLoader

logger = getLogger(__name__)


def write_batch(read_names: Iterable[str],
                out_dir: Path, sample: str, ref: str, sect: str, batch: int):
    """ Write the names of the reads in one batch to a file. """
    # Determine the path to the batch file.
    batch_file = MaskReport.build_batch_path(out_dir, batch,
                                             sample=sample, ref=ref, sect=sect,
                                             ext=path.CSVZIP_EXT)
    # Write the read names to the batch file.
    read_data = pd.Series(read_names)
    read_data.to_csv(batch_file, index=False)
    logger.debug(f"Wrote {read_data.size} read names to {batch_file}:\n"
                 f"{read_data}")
    # Compute the checksum of the batch file.
    return digest_file(batch_file)


class BitMasker(object):
    """ Mask bit vectors. """

    MASK_READ_FINFO = "read-finfo"
    MASK_READ_FMUT = "read-fmut"
    MASK_READ_GAP = "read-gap"
    MASK_POS_NINFO = "pos-ninfo"
    MASK_POS_FMUT = "pos-fmut"

    def __init__(self,
                 loader: RelateLoader,
                 section: Section,
                 bit_caller: BitCaller, *,
                 exclude_polya: int = 0,
                 exclude_gu: bool = False,
                 exclude_pos: Iterable[tuple[str, int]] = (),
                 min_mut_gap: int = 0,
                 min_finfo_read: float = 0.,
                 max_fmut_read: float = 1.,
                 min_ninfo_pos: int = 0,
                 max_fmut_pos: float = 1.):
        """
        Parameters
        ----------
        loader: RelateLoader
            Relation vector loader
        section: Section
            Load this section from the reference
        bit_caller: BitCaller
            Bit caller
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
        self.out_dir = loader.out_dir
        self.sample = loader.sample
        self.ref = loader.ref
        # Create a new section to compute the masked positions.
        self.section = Section(loader.ref, loader.get_refseq(),
                               end5=section.end5, end3=section.end3,
                               name=section.name)
        # Set the parameters for excluding positions from the section.
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = np.array([pos for ref, pos in exclude_pos
                                     if ref == self.ref],
                                    dtype=int)
        # Set the parameters for filtering reads.
        self.min_mut_gap = min_mut_gap
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        # Set the parameters for filtering positions.
        if not min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {min_ninfo_pos}")
        self.min_ninfo_pos = min_ninfo_pos
        if not 0. <= max_fmut_pos < 1.:
            raise ValueError(
                f"max_fmut_pos must be in [0, 1), but got {max_fmut_pos}")
        self.max_fmut_pos = max_fmut_pos
        # Exclude positions based on the parameters.
        self.section.mask_polya(self.exclude_polya)
        self.section.mask_gu(self.exclude_gu)
        self.section.mask_pos(self.exclude_pos)
        # Create a new BitCaller with the masked section.
        self.bit_caller = BitCaller(self.section,
                                    bit_caller.affi_call,
                                    bit_caller.anti_call)
        # For each batch, load the positions that remain unmasked,
        # count the informative and mutated bits in every vector and
        # position, and drop reads that do not pass the filters.
        self.counter = BitCounter(self.section, self.bit_caller.iter(
            loader.iter_batches_processed(positions=self.pos_kept),
            masks={self.MASK_READ_FINFO: self._mask_min_finfo_read,
                   self.MASK_READ_FMUT: self._mask_max_fmut_read,
                   self.MASK_READ_GAP: self._mask_min_mut_gap}
        ))
        # Write batches of read names and record their checksums.
        self.checksums = [write_batch(read_names, self.out_dir, self.sample,
                                      self.ref, self.section.name, batch)
                          for batch, read_names
                          in enumerate(self.counter.read_batches)]
        # Warn if no reads were counted.
        if self.counter.nreads == 0:
            logger.warning(f"No reads passed through {self}")
        # Find the positions with insufficient informative reads.
        mask_min_ninfo = self.counter.n_info_per_pos < self.min_ninfo_pos
        # Mask the positions with insufficient informative reads.
        self.section.add_mask(self.MASK_POS_NINFO,
                              index_to_pos(mask_min_ninfo.index[mask_min_ninfo]))
        # Find the positions with excessive mutation fractions.
        mask_max_fmut = self.counter.f_affi_per_pos > self.max_fmut_pos
        # Mask the positions with excessive mutation fractions.
        self.section.add_mask(self.MASK_POS_FMUT,
                              index_to_pos(mask_max_fmut.index[mask_max_fmut]))
        # Check if any positions passed the filters.
        if not np.any(self.section.unmasked_bool):
            logger.warning(f"No positions passed through {self}")

    @property
    def n_reads_init(self):
        return self.counter.nreads_given

    @property
    def n_reads_min_finfo(self):
        return self.counter.nmasked[self.MASK_READ_FINFO]

    @property
    def n_reads_max_fmut(self):
        return self.counter.nmasked[self.MASK_READ_FMUT]

    @property
    def n_reads_min_gap(self):
        return self.counter.nmasked[self.MASK_READ_GAP]

    @property
    def n_reads_kept(self):
        """ Number of reads kept. """
        return self.counter.nreads

    @property
    def pos_gu(self):
        """ Positions masked for having a G or U base. """
        return self.section.get_mask(self.section.MASK_GU)

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

    def _mask_min_finfo_read(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with too few informative bits. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be in [0, 1], but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # Nothing to mask.
            return np.zeros_like(batch.reads, dtype=bool)
        return batch.f_info_per_read.values < self.min_finfo_read

    def _mask_max_fmut_read(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be in [0, 1], but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read == 1.:
            # Nothing to mask.
            return np.zeros_like(batch.reads, dtype=bool)
        return batch.f_affi_per_read.values > self.max_fmut_read

    def _mask_min_mut_gap(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        # Initially, mask no reads (set all to False).
        read_mask = np.zeros_like(batch.reads, dtype=bool)
        if self.min_mut_gap == 0:
            # Nothing to mask.
            return read_mask
        # Make a mask for the current window of positions.
        positions = batch.section.unmasked_int
        pos_mask = pd.Series(False, index=positions)
        # Loop through every 5' position in the bit vectors.
        for pos5 in positions:
            # Find the 3'-most position that must be checked.
            pos3 = pos5 + self.min_mut_gap
            # Mask all positions between pos5 and pos3, inclusive.
            # For example, if positions 1, 2, 4, 5, and 7 remain and
            # min_gap = 3, then the positions checked in each
            # iteration are [1, 2, 4], [2, 4, 5], [4, 5], [4, 5, 7].
            pos_mask.loc[pos5: pos3] = True
            # Sum over axis 1 to count the mutations in each read
            # between positions pos5 and pos3, inclusive. If a read
            # has > 1 mutation in this range, then it has mutations
            # that are too close. Set all such reads' flags to True.
            read_mask |= batch.affi.loc[:, pos_mask.values].sum(axis=1) > 1
            # Erase the mask over the current window of positions.
            pos_mask.loc[pos5: pos3] = False
            if pos3 >= positions[-1]:
                # The end of the section has been reached. Stop.
                break
        return read_mask

    def create_report(self):
        """ Create a FilterReport from a BitFilter object. """
        return MaskReport(
            out_dir=self.out_dir,
            sample=self.sample,
            ref=self.ref,
            sect=self.section.name,
            end5=self.section.end5,
            end3=self.section.end3,
            checksums=self.checksums,
            n_batches=self.n_batches,
            count_refs=self.bit_caller.anti_call,
            count_muts=self.bit_caller.affi_call,
            exclude_gu=self.exclude_gu,
            exclude_polya=self.exclude_polya,
            exclude_pos=self.exclude_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            n_pos_init=self.section.length,
            n_pos_gu=len(self.pos_gu),
            n_pos_polya=len(self.pos_polya),
            n_pos_user=len(self.pos_user),
            n_pos_min_ninfo=len(self.pos_min_ninfo),
            n_pos_max_fmut=len(self.pos_max_fmut),
            n_pos_kept=len(self.pos_kept),
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
        )

    def __str__(self):
        return f"Mask {self.sample} over {self.section} with {self.bit_caller}"


def mask_section(loader: RelateLoader,
                 section: Section,
                 count_del: bool,
                 count_ins: bool,
                 discount: Iterable[str], *,
                 rerun: bool,
                 **kwargs):
    """ Filter a section of a set of bit vectors. """
    # Check if the report file already exists.
    report_file = MaskReport.build_path(loader.out_dir,
                                        sample=loader.sample,
                                        ref=loader.ref,
                                        sect=section.name)
    if rerun or not report_file.is_file():
        # Create the bit caller.
        bit_caller = BitCaller.from_counts(section, count_del, count_ins,
                                           discount)
        # Create and apply the mask.
        bit_masker = BitMasker(loader, section, bit_caller, **kwargs)
        # Output the results of filtering.
        report = bit_masker.create_report()
        report.save()
    else:
        logger.warning(f"File exists: {report_file}")
    return report_file
