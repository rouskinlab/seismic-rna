"""

Relation Vector Writing Module
========================================================================

Given alignment map (BAM) files, split each file into batches of reads,
write the relation vectors for each batch to a compressed file, and
write a report summarizing the results.

"""

from datetime import datetime
from logging import getLogger
from pathlib import Path
from typing import Iterable

from .io import from_reads, QnamesBatchIO, RelateBatchIO
from .py.relate import find_rels_line
from .report import RelateReport
from .sam import XamViewer
from ..core import path
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.ngs import encode_phred
from ..core.io import RefseqIO
from ..core.seq import DNA, get_fasta_seq
from ..core.write import need_write

logger = getLogger(__name__)


def mib_to_bytes(batch_size: float):
    """
    Return the number of bytes per batch of a given size in mebibytes.

    Parameters
    ----------
    batch_size: float
        Size of the batch in mebibytes (MiB): 1 MiB = 2^20 bytes

    Return
    ------
    int
        Number of bytes per batch, to the nearest integer
    """
    if batch_size <= 0.:
        raise ValueError(f"batch_size must be > 0, but got {batch_size}")
    nbytes = round(batch_size * 2 ** 20)
    if nbytes <= 0:
        logger.warning(f"Using batch_size of {batch_size} MiB gave {nbytes} "
                       f"bytes per batch: defaulting to 1")
        nbytes = 1
    return nbytes


def get_reads_per_batch(bytes_per_batch: int, seq_len: int):
    """ Compute the number of reads per batch. """
    reads_per_batch = bytes_per_batch // seq_len
    if reads_per_batch < 1:
        logger.warning(f"Cannot have {bytes_per_batch} bytes per batch with a "
                       f"sequence of {seq_len} nt. Using 1 read per batch.")
        return 1
    return reads_per_batch


def generate_batch(batch: int, *,
                   xam_view: XamViewer,
                   out_dir: Path,
                   refseq: DNA,
                   min_mapq: int,
                   min_qual: str,
                   ambrel: bool,
                   brotli_level: int):
    """ Compute relation vectors for every SAM record in one batch,
    write the vectors to a batch file, and return its MD5 checksum
    and the number of vectors. """
    logger.info(f"Began computing relation vectors for batch {batch} "
                f"of {xam_view}")

    def relate_records(records: Iterable[tuple[str, str]]):
        for line1, line2 in records:
            try:
                yield find_rels_line(line1,
                                     line2,
                                     xam_view.ref,
                                     refseq,
                                     min_mapq,
                                     min_qual,
                                     ambrel)
            except Exception as error:
                logger.error(f"Failed to compute relation vector: {error}")

    names, relvecs = from_reads(relate_records(xam_view.iter_records(batch)),
                                xam_view.sample,
                                xam_view.ref,
                                len(refseq),
                                batch)
    logger.info(f"Ended computing relation vectors for batch {batch} "
                f"of {xam_view}")
    _, relv_check = relvecs.save(out_dir, brotli_level, overwrite=True)
    _, name_check = names.save(out_dir, brotli_level, overwrite=True)
    return relvecs.num_reads, relv_check, name_check


class RelationWriter(object):
    """
    Compute and write relation vectors for all reads from one sample
    mapped to one reference sequence.
    """

    def __init__(self, xam_view: XamViewer, seq: DNA):
        self.xam = xam_view
        self.seq = seq

    @property
    def sample(self):
        return self.xam.sample

    @property
    def ref(self):
        return self.xam.ref

    def _write_report(self, *, out_dir: Path, **kwargs):
        report = RelateReport(sample=self.sample,
                              ref=self.ref,
                              **kwargs)
        return report.save(out_dir, overwrite=True)

    def _write_refseq(self, out_dir: Path, brotli_level: int):
        """ Write the reference sequence to a file. """
        refseq_file = RefseqIO(sample=self.sample,
                               ref=self.ref,
                               refseq=self.seq)
        _, checksum = refseq_file.save(out_dir, brotli_level, overwrite=True)
        return checksum

    def _generate_batches(self, *,
                          out_dir: Path,
                          keep_temp: bool,
                          min_mapq: int,
                          phred_enc: int,
                          min_phred: int,
                          ambrel: bool,
                          brotli_level: int,
                          n_procs: int):
        """ Compute a relation vector for every record in a XAM file,
        split among one or more batches. For each batch, write a matrix
        of the vectors to one batch file, and compute its checksum. """
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        # The SAM file will remain open until exiting the with block.
        logger.info(f"Began {self}")
        try:
            # Collect the keyword arguments.
            disp_kwargs = dict(xam_view=self.xam,
                               out_dir=out_dir,
                               refseq=self.seq,
                               min_mapq=min_mapq,
                               min_qual=encode_phred(min_phred, phred_enc),
                               ambrel=ambrel,
                               brotli_level=brotli_level)
            # Generate and write relation vectors for each batch.
            results = dispatch(generate_batch,
                               n_procs,
                               parallel=True,
                               pass_n_procs=False,
                               args=as_list_of_tuples(self.xam.indexes),
                               kwargs=disp_kwargs)
            if results:
                nums_reads, relv_checks, name_checks = map(list,
                                                           zip(*results,
                                                               strict=True))
            else:
                nums_reads = list()
                relv_checks = list()
                name_checks = list()
            n_reads = sum(nums_reads)
            n_batches = len(nums_reads)
            checksums = {RelateBatchIO.btype(): relv_checks,
                         QnamesBatchIO.btype(): name_checks}
            logger.info(f"Ended {self}: {n_reads} reads in {n_batches} batches")
            return n_reads, n_batches, checksums
        finally:
            if not keep_temp:
                # Delete the temporary SAM file before exiting.
                self.xam.delete_temp_sam()

    def write(self, *,
              out_dir: Path,
              brotli_level: int,
              force: bool,
              **kwargs):
        """ Compute a relation vector for every record in a BAM file,
        write the vectors into one or more batch files, compute their
        checksums, and write a report summarizing the results. """
        report_file = RelateReport.build_path(top=out_dir,
                                              sample=self.sample,
                                              ref=self.ref)
        if need_write(report_file, force):
            began = datetime.now()
            # Write the reference sequence to a file.
            refcheck = self._write_refseq(out_dir, brotli_level)
            # Compute relation vectors and time how long it takes.
            (nreads,
             nbats,
             checks) = self._generate_batches(out_dir=out_dir,
                                              brotli_level=brotli_level,
                                              **kwargs)
            ended = datetime.now()
            # Write a report of the relation step.
            self._write_report(out_dir=out_dir,
                               n_reads_rel=nreads,
                               n_batches=nbats,
                               checksums=checks,
                               refseq_checksum=refcheck,
                               began=began,
                               ended=ended)
        return report_file

    def __str__(self):
        return f"Relate {self.xam}"


def write_one(xam_file: Path, *,
              fasta: Path,
              temp_dir: Path,
              batch_size: float,
              min_reads: int = 0,
              **kwargs):
    """ Write the batches of relation vectors for one XAM file. """
    # Get the reference sequence.
    ref = path.parse(xam_file, *path.XAM_SEGS)[path.REF]
    seq = get_fasta_seq(fasta, DNA, ref)
    # Compute the number of records per batch.
    rec_per_batch = get_reads_per_batch(mib_to_bytes(batch_size), len(seq))
    # Determine if there are enough reads.
    xam = XamViewer(xam_file, temp_dir, rec_per_batch)
    if min_reads > 0 and xam.n_reads < min_reads:
        raise ValueError(
            f"Insufficient reads in {xam}: {xam.n_reads} (< {min_reads})")
    # Write the batches.
    writer = RelationWriter(xam, seq)
    return writer.write(**kwargs)


def write_all(xam_files: Iterable[Path],
              max_procs: int,
              parallel: bool,
              **kwargs):
    """  """
    return dispatch(write_one,
                    max_procs,
                    parallel,
                    args=as_list_of_tuples(path.deduplicated(xam_files)),
                    kwargs=kwargs)

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
