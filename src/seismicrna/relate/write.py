"""

Relation Vector Writing Module
========================================================================

Given alignment map (BAM) files, split each file into batches of reads,
write the relation vectors for each batch to a compressed file, and
write a report summarizing the results.

"""

from datetime import datetime
from itertools import starmap
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .blank import blank_relvec
from .relate import relate_line, relate_pair
from .report import RelateReport
from .sam import XamViewer
from .seqpos import format_seq_pos
from ..core import path
from ..core.fasta import get_fasta_seq
from ..core.files import digest_file
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.qual import encode_phred
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


def get_records_per_batch(bytes_per_batch: int, seq_len: int):
    """ Compute the number of records per batch. """
    records_per_batch = bytes_per_batch // seq_len
    if records_per_batch < 1:
        logger.warning(f"Cannot have {bytes_per_batch} bytes per batch with a "
                       f"sequence of {seq_len} nt. Using 1 record per batch.")
        return 1
    return records_per_batch


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


def write_batch(batch: int,
                relvecs: pd.DataFrame,
                sample: str,
                ref: str,
                out_dir: Path):
    """ Write a batch of relation vectors to an ORC file. """
    logger.info(
        f"Began writing sample '{sample}' reference '{ref}' batch {batch}")
    # Build the path to the batch file.
    batch_path = RelateReport.build_batch_path(out_dir, batch, sample=sample,
                                               ref=ref, ext=path.PARQ_EXTS[0])
    relvecs.to_parquet(batch_path, index=True, engine="pyarrow",
                       compression="brotli")
    logger.info(f"Ended writing sample '{sample}' reference '{ref}' "
                f"batch {batch} to {batch_path}")
    return batch_path


def _relate_record(relvec: bytearray,
                   line1: str,
                   line2: str,
                   ref: str,
                   refseq: DNA,
                   min_qual: str,
                   ambrel: bool):
    """ Compute the relation vector of a record in a SAM file. """
    # Fill the relation vector with data from the SAM line(s).
    if line2:
        relate_pair(relvec, line1, line2,
                    ref, refseq, len(refseq), min_qual, ambrel)
    else:
        relate_line(relvec, line1, ref, refseq, len(refseq), min_qual, ambrel)


def _relate_batch(batch: int, *,
                  xam_view: XamViewer,
                  out_dir: Path,
                  refseq: DNA,
                  min_qual: str,
                  ambrel: bool):
    """ Compute relation vectors for every SAM record in one batch,
    write the vectors to a batch file, and return its MD5 checksum
    and the number of vectors. """
    logger.info(f"Began computing relation vectors for batch {batch}"
                f"of {xam_view}")
    # Cache a blank relation vector as a template.
    blank = bytearray(blank_relvec(len(refseq)).tobytes())

    # Wrap self._relate_record with keyword arguments and a
    # try-except block so that if one record fails to vectorize,
    # it does not crash all the others.

    def relate_record(read_name: str, line1: str, line2: str):
        # Copy the blank template to get a new relation vector.
        relvec = blank.copy()
        try:
            _relate_record(relvec, line1, line2, xam_view.ref, refseq,
                           min_qual, ambrel)
        except Exception as err:
            logger.error(f"Failed to relate read '{read_name}': {err}")
            # Return an empty read name and relation vector.
            return "", bytearray()
        else:
            return read_name, relvec

    # Vectorize every record in the batch.
    read_names, relbytes = zip(*starmap(relate_record,
                                        xam_view.iter_records(batch)))
    # For every read for which creating a relation vector failed, an
    # empty string was returned as the read name and an empty
    # bytearray as the relation vector. The empty read names must be
    # filtered out, while the empty relation vectors will not cause
    # problems because, being of length zero, they will disappear
    # when concatenated with the other vectors into a 1D array.
    read_names = list(filter(None, read_names))
    # Compute the number of reads that passed and failed.
    n_total = len(relbytes)  # has empty bytearray for each failed read
    n_pass = len(read_names)  # has no item for any failed read
    n_fail = n_total - n_pass  # difference between total and passed
    if not n_pass:
        logger.warning(f"Batch {batch} of {xam_view} yielded 0 vectors")
    # Write the names and vectors to a file.
    batch_file = write_batch(batch,
                             _bytes_to_df(relbytes, read_names, refseq),
                             sample=xam_view.sample,
                             ref=xam_view.ref,
                             out_dir=out_dir)
    # Compute the MD5 checksum of the file.
    checksum = digest_file(batch_file)
    logger.info(f"Ended computing relation vectors for batch {batch} "
                f"of {xam_view}")
    return n_pass, n_fail, checksum


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
            disp_kwargs = dict(xam_view=self.xam, out_dir=out_dir,
                               refseq=self.seq, ambrel=ambrel,
                               min_qual=encode_phred(min_phred, phred_enc))
            # Generate and write relation vectors for each batch.
            results = dispatch(_relate_batch,
                               n_procs,
                               parallel=True,
                               pass_n_procs=False,
                               args=as_list_of_tuples(self.xam.indexes),
                               kwargs=disp_kwargs)
            # The list of results contains, for each batch, a tuple of
            # the number of relation vectors in the batch and the MD5
            # checksum of the batch file. Compute the total number of
            # vectors and list all the checksums.
            n_pass = sum(result[0] for result in results)
            n_fail = sum(result[1] for result in results)
            checksums: list[str] = [result[2] for result in results]
            logger.info(f"Ended {self}: {n_pass} pass, {n_fail} fail")
            return n_pass, n_fail, checksums
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
            n_pass, n_fail, checksums = self._write_batches(out_dir=out_dir,
                                                            **kwargs)
            ended = datetime.now()
            # Write a report of the relation step.
            self._write_report(out_dir=out_dir,
                               n_reads_rel_pass=n_pass,
                               n_reads_rel_fail=n_fail,
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
    rec_per_batch = get_records_per_batch(mib_to_bytes(batch_size), len(seq))
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
