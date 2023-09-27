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

import numpy as np
import pandas as pd

from .relate import find_rels_line
from .report import RelateReport
from .sam import XamViewer
from .seqpos import format_seq_pos
from ..core import path
from ..core.fasta import get_fasta_seq
from ..core.files import digest_file
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.qual import encode_phred
from ..core.array import from_reads
from ..core.seq import DNA

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
    return round(batch_size * 2 ** 20)


def get_reads_per_batch(bytes_per_batch: int, seq_len: int):
    """ Compute the number of reads per batch. """
    reads_per_batch = bytes_per_batch // seq_len
    if reads_per_batch < 1:
        logger.warning(f"Cannot have {bytes_per_batch} bytes per batch with a "
                       f"sequence of {seq_len} nt. Using 1 read per batch.")
        return 1
    return reads_per_batch


def _bytes_to_df(relbytes: tuple[bytearray, ...], reads: list[str], refseq: DNA):
    """ Convert the relation vectors from bytearrays to a DataFrame. """
    # Ideally, this step would use the NumPy unsigned 8-bit integer
    # (np.uint8) data type because the data must be read back from the
    # file as this type. But the PyArrow backend of to_parquet does not
    # currently support uint8, so we are using np.byte, which works.
    relarray = np.frombuffer(b"".join(relbytes), dtype=np.byte)
    relarray.shape = relarray.size // len(refseq), len(refseq)
    # Determine the numeric positions in the reference sequence.
    positions = np.arange(1, len(refseq) + 1)
    # Parquet format requires that the label of each column be a string.
    # This requirement provides a good opportunity to add the reference
    # sequence into the column labels themselves. Subsequent tasks can
    # then obtain the entire reference sequence without needing to read
    # the report file: relate.export.as_iter(), for example.
    columns = format_seq_pos(refseq, positions, 1)
    # Data must be converted to pd.DataFrame for PyArrow to write.
    # Set copy=False to prevent copying the relation vectors.
    return pd.DataFrame(data=relarray, index=reads, columns=columns, copy=False)


def relate_batch(batch: int, *,
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

    relvecs = from_reads(relate_records(xam_view.iter_records(batch)),
                         batch,
                         refseq)
    batch_file = RelateReport.build_batch_path(out_dir,
                                               batch,
                                               sample=xam_view.sample,
                                               ref=xam_view.ref,
                                               ext=path.PICKLE_BROTLI_EXT)
    relvecs.save(batch_file, brotli_level=brotli_level, overwrite=True)
    checksum = digest_file(batch_file)
    logger.info(f"Ended computing relation vectors for batch {batch} "
                f"of {xam_view}")
    return relvecs.num_reads, checksum


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
        report = RelateReport(out_dir=out_dir,
                              seq=self.seq,
                              sample=self.sample,
                              ref=self.ref,
                              **kwargs)
        report.save()
        return report.get_path()

    def _write_batches(self, *,
                       out_dir: Path,
                       save_temp: bool,
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
        logger.info(f"Began running {self}")
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
            results = dispatch(relate_batch,
                               n_procs,
                               parallel=True,
                               pass_n_procs=False,
                               args=as_list_of_tuples(self.xam.indexes),
                               kwargs=disp_kwargs)
            # The list of results contains, for each batch, a tuple of
            # the number of relation vectors in the batch and the MD5
            # checksum of the batch file. Compute the total number of
            # vectors and list all the checksums.
            nreads = sum(result[0] for result in results)
            checksums: list[str] = [result[1] for result in results]
            logger.info(f"Ended {self}: {nreads} reads")
            return nreads, checksums
        finally:
            if not save_temp:
                # Delete the temporary SAM file before exiting.
                self.xam.delete_temp_sam()

    def write(self, *, rerun: bool, out_dir: Path, **kwargs):
        """ Compute a relation vector for every record in a BAM file,
        write the vectors into one or more batch files, compute their
        checksums, and write a report summarizing the results. """
        report_file = RelateReport.build_path(out_dir,
                                              sample=self.sample,
                                              ref=self.ref)
        # Check if the report file already exists.
        if rerun or not report_file.is_file():
            # Compute relation vectors and time how long it takes.
            began = datetime.now()
            nreads, checksums = self._write_batches(out_dir=out_dir, **kwargs)
            ended = datetime.now()
            # Write a report of the relation step.
            self._write_report(out_dir=out_dir,
                               n_reads_rel=nreads,
                               checksums=checksums,
                               began=began,
                               ended=ended)
        else:
            logger.warning(f"File exists: {report_file}")
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
        raise ValueError(f"Insufficient reads in {xam}: {xam.n_reads}")
    # Write the batches.
    writer = RelationWriter(xam, seq)
    return writer.write(**kwargs)


def write_all(xam_files: Iterable[Path],
              max_procs: int,
              parallel: bool,
              **kwargs):
    """  """
    return dispatch(write_one, max_procs, parallel,
                    args=as_list_of_tuples(xam_files),
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
