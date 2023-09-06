"""
Relate -- Write Module
======================
Auth: Matty

Given alignment map (BAM) files, split each file into batches of reads,
write the relation vectors for each batch to a compressed file, and
write a report summarizing the results.
"""

from __future__ import annotations

from datetime import datetime
from functools import cached_property
from itertools import starmap
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .relate import relate_line, relate_pair
from .report import RelateReport
from .sam import iter_batch_indexes, iter_records
from .seqpos import format_seq_pos
from ..core import path
from ..core.cmd import CMD_REL
from ..core.files import digest_file
from ..core.parallel import dispatch
from ..core.rel import blank_relvec
from ..core.seq import DNA, parse_fasta
from ..core.xam import count_total_reads, run_flagstat, run_view_xam

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
    return round(batch_size * 2**20)


def write_batch(batch: int,
                relvecs: tuple[bytearray, ...],
                read_names: list[str], *,
                sample: str,
                ref: str,
                seq: DNA,
                out_dir: Path):
    """ Write a batch of relation vectors to an ORC file. """
    logger.info(
        f"Began writing sample '{sample}' reference '{ref}' batch {batch}")
    # Process the relation vectors into a 2D NumPy array (matrix).
    # Ideally, this step would use the NumPy unsigned 8-bit integer
    # (np.uint8) data type because the data must be read back from the
    # file as this type. But the PyArrow backend of to_parquet does not
    # currently support uint8, so we are using np.byte, which works.
    relmatrix = np.frombuffer(b"".join(relvecs), dtype=np.byte)
    relmatrix.shape = relmatrix.size // len(seq), len(seq)
    # Determine the numeric positions in the reference sequence.
    positions = np.arange(1, len(seq) + 1)
    # Parquet format requires that the label of each column be a string.
    # This requirement provides a good opportunity to add the reference
    # sequence into the column labels themselves. Subsequent tasks can
    # then obtain the entire reference sequence without needing to read
    # the report file: relate.export.as_iter(), for example.
    columns = format_seq_pos(seq, positions, 1)
    # Data must be converted to pd.DataFrame for PyArrow to write.
    # Set copy=False to prevent copying the relation vectors.
    relframe = pd.DataFrame(data=relmatrix,
                            index=read_names,
                            columns=columns,
                            copy=False)
    # Build the path to the batch file.
    batch_path = RelateReport.build_batch_path(out_dir, batch, sample=sample,
                                               ref=ref, ext=path.PARQ_EXTS[0])
    relframe.to_parquet(batch_path, index=True, engine="pyarrow",
                        compression="brotli")
    logger.info(f"Ended writing sample '{sample}' reference '{ref}' "
                f"batch {batch} to {batch_path}")
    return batch_path


def _relate_record(relvec: bytearray, line1: str, line2: str, refseq: DNA, *,
                   min_qual: str, ambrel: bool):
    """ Compute the relation vector of a record in a SAM file. """
    # Fill the relation vector with data from the SAM line(s).
    if line2:
        relate_pair(relvec, line1, line2, refseq, len(refseq), min_qual, ambrel)
    else:
        relate_line(relvec, line1, refseq, len(refseq), min_qual, ambrel)


def _relate_batch(batch: int, start: int, stop: int, *,
                  temp_sam: Path, out_dir: Path,
                  sample: str, ref: str, refseq: DNA,
                  min_qual: str, ambrel: bool):
    """ Compute relation vectors for every SAM record in one batch,
    write the vectors to a batch file, and return its MD5 checksum
    and the number of vectors. """
    logger.info(f"Began computing relation vectors for batch {batch} of "
                f"{temp_sam} (file indexes {start} - {stop})")
    # Cache a blank relation vector as a template.
    blank = bytearray(blank_relvec(len(refseq)).tobytes())

    # Wrap self._relate_record with keyword arguments and a
    # try-except block so that if one record fails to vectorize,
    # it does not crash all the others.

    def relate_record(read_name: str, line1: str, line2: str):
        # Copy the blank template to get a new relation vector.
        relvec = blank.copy()
        try:
            _relate_record(relvec, line1, line2, refseq,
                           min_qual=min_qual, ambrel=ambrel)
        except Exception as err:
            logger.error(f"Failed to relate read '{read_name}': {err}")
            # Return an empty read name and relation vector.
            return "", bytearray()
        else:
            return read_name, relvec

    with open(temp_sam) as sam_file:
        # Vectorize every record in the batch.
        records = iter_records(sam_file, start, stop)
        read_names, relvecs = zip(*starmap(relate_record, records))
        # For every read for which creating a relation vector failed, an
        # empty string was returned as the read name and an empty
        # bytearray as the relation vector. The empty read names must be
        # filtered out, while the empty relation vectors will not cause
        # problems because, being of length zero, they will disappear
        # when concatenated with the other vectors into a 1D array.
        read_names = list(filter(None, read_names))
    # Compute the number of reads that passed and failed.
    n_total = len(relvecs)  # has an empty byte for each failed read
    n_pass = len(read_names)  # has no item for any failed read
    n_fail = n_total - n_pass  # difference between total and passed
    if not n_pass:
        logger.warning(f"Batch {batch} of {temp_sam} yielded 0 vectors")
    # Write the names and vectors to a file.
    batch_file = write_batch(batch, relvecs, read_names, sample=sample,
                             ref=ref, seq=refseq, out_dir=out_dir)
    # Compute the MD5 checksum of the file.
    checksum = digest_file(batch_file)
    logger.info(f"Ended computing relation vectors for batch {batch} of "
                f"{temp_sam} (file indexes {start} - {stop})")
    return n_pass, n_fail, checksum


class RelationWriter(object):
    """
    Compute and write relation vectors for all reads from one sample
    mapped to one reference sequence.
    """

    def __init__(self, bam_file: Path, seq: DNA):
        self.bam = bam_file
        self.seq = seq

    @cached_property
    def sample_ref(self):
        fields = path.parse(self.bam, *path.XAM_SEGS)
        return fields[path.SAMP], fields[path.REF]

    @property
    def sample(self):
        return self.sample_ref[0]

    @property
    def ref(self):
        return self.sample_ref[1]

    def _write_report(self, *, out_dir: Path, **kwargs):
        report = RelateReport(out_dir=out_dir,
                              seq=self.seq,
                              sample=self.sample,
                              ref=self.ref,
                              **kwargs)
        report.save()
        return report.get_path()

    def _relate_bam(self, *,
                    out_dir: Path,
                    temp_dir: Path,
                    save_temp: bool,
                    min_mapq: int,
                    batch_size: int,
                    phred_enc: int,
                    min_phred: int,
                    ambrel: bool,
                    n_procs: int):
        """ Compute a relation vector for every record in a BAM file,
        split among one or more batches. For each batch, write a matrix
        of the vectors to one batch file, and compute its checksum. """
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        # The SAM file will remain open until exiting the with block.
        logger.info(f"Began running {self}")
        # Determine the path of the temporary SAM file.
        temp_sam = path.build(*path.XAM_STEP_SEGS,
                              top=temp_dir, sample=self.sample,
                              cmd=CMD_REL, step=path.STEPS_VECT[0],
                              ref=self.ref, ext=path.SAM_EXT)
        # Create the temporary SAM file.
        run_view_xam(self.bam, temp_sam, min_mapq=min_mapq, n_procs=n_procs)
        try:
            with open(temp_sam) as sam_file:
                # Compute the number of records per batch.
                n_per_bat = max(1, mib_to_bytes(batch_size) // len(self.seq))
                # Compute the batch indexes.
                disp_args = list(iter_batch_indexes(sam_file, n_per_bat))
                # Collect the keyword arguments.
                disp_kwargs = dict(temp_sam=temp_sam, out_dir=out_dir,
                                   sample=self.sample, ref=self.ref,
                                   refseq=self.seq, ambrel=ambrel,
                                   min_qual=get_min_qual(min_phred, phred_enc))
                # Generate and write relation vectors for each batch.
                results = dispatch(_relate_batch, n_procs,
                                   parallel=True, pass_n_procs=False,
                                   args=disp_args, kwargs=disp_kwargs)
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
                temp_sam.unlink(missing_ok=True)

    def relate_sample_ref(self, *, rerun: bool, out_dir: Path, **kwargs):
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
            n_pass, n_fail, checksums = self._relate_bam(out_dir=out_dir,
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
        return f"Relate {self.bam}"


def get_min_qual(min_phred: int, phred_enc: int):
    """
    Return the minimum quality for a base in a read to be considered
    informative, as the character in the FASTQ file encoding that would
    be the minimum valid quality.

    Parameters
    ----------
    min_phred: int
        The minimum Phred score needed to use the value of a base call.
    phred_enc: int
        The encoding offset for Phred scores. A Phred score is encoded
        as the character whose ASCII value is the sum of the phred score
        and the encoding offset.

    Return
    ------
    str
        The character whose ASCII code, in the encoding scheme of the
        FASTQ file, represents the minimum valid quality.

    Examples
    --------
    For example, if the minimum Phred score (`min_phred`) that
    is accepted as informative is 20, and the Phred encoding of the
    FASTQ file (`phred_enc`) is 33 (i.e. ASCII+33), then the
    minimum quality as an ASCII integer (`min_qual`) is 20 + 33
    = 53, which is character '5'. If `min_phred` were 37, then
    `min_qual` would be 37 + 33 = 70, which is character 'F'.
    """
    return chr(min_phred + phred_enc)


def get_sufficient_xam_files(xam_files: Iterable[Path],
                             min_reads: int,
                             n_procs: int):
    """ Yield the XAM files with enough reads. """
    for xam_file in xam_files:
        n_reads = count_total_reads(run_flagstat(xam_file, n_procs=n_procs))
        if n_reads >= min_reads:
            yield xam_file
        else:
            logger.warning(f"Skipping {xam_file} with {n_reads} reads, which "
                           f"is less than the minimum ({min_reads})")


def get_relaters(xam_files: Iterable[Path], fasta: Path, *,
                 min_reads: int, max_procs: int):
    """ Return a RelationWriter for every XAM file with a reference in
    the given FASTA file and with a sufficient number of reads. """
    logger.info("Began creating relation writers")
    # Filter out the XAM files with insufficient reads and index them by
    # the name of the reference.
    xam_files = {path.parse(xam, *path.XAM_SEGS)[path.REF]: xam for xam in
                 get_sufficient_xam_files(xam_files, min_reads, max_procs)}
    # Cache the sequences of the references for those XAM files.
    seqs = {ref: seq for ref, seq in parse_fasta(fasta) if ref in xam_files}
    # Create a RelationWriter for each XAM file.
    relaters = [RelationWriter(xam_file, seqs[ref])
                for ref, xam_file in xam_files.items()]
    logger.info(f"Ended creating {len(relaters)} relation writer(s)")
    return relaters


def relate_all(relaters: list[RelationWriter], *,
               max_procs: int, parallel: bool, **kwargs) -> list[Path]:
    """ Run one or more RelationWriters in series or parallel. """
    logger.info("Began generating relation vector sets")
    # Determine method of parallelization. Do not use hybrid mode, which
    # would try to process multiple SAM files in parallel and use more
    # than one processe for each file.
    return dispatch([relater.relate_sample_ref for relater in relaters],
                    max_procs, parallel, kwargs=kwargs)
