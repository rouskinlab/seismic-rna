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
from ..core.io import RefseqIO
from ..core.ngs import encode_phred
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.seq import DNA, get_fasta_seq
from ..core.write import need_write

logger = getLogger(__name__)


def generate_batch(batch: int, *,
                   xam_view: XamViewer,
                   top: Path,
                   refseq: DNA,
                   min_mapq: int,
                   min_qual: str,
                   ambindel: bool,
                   overhangs: bool,
                   clip_end5: int,
                   clip_end3: int,
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
                                     ambindel,
                                     overhangs,
                                     clip_end5,
                                     clip_end3)
            except Exception as error:
                logger.error(f"Failed to compute relation vector: {error}")

    names, relvecs = from_reads(relate_records(xam_view.iter_records(batch)),
                                xam_view.sample,
                                xam_view.ref,
                                refseq,
                                batch)
    logger.info(f"Ended computing relation vectors for batch {batch} "
                f"of {xam_view}")
    _, relv_check = relvecs.save(top, brotli_level)
    _, name_check = names.save(top, brotli_level)
    return relvecs.num_reads, relv_check, name_check


class RelationWriter(object):
    """
    Compute and write relation vectors for all reads from one sample
    mapped to one reference sequence.
    """

    def __init__(self, xam_view: XamViewer, seq: DNA):
        self._xam = xam_view
        self.seq = seq

    @property
    def sample(self):
        return self._xam.sample

    @property
    def ref(self):
        return self._xam.ref

    @property
    def num_reads(self):
        return self._xam.n_reads

    def _write_report(self, *, top: Path, **kwargs):
        report = RelateReport(sample=self.sample, ref=self.ref, **kwargs)
        return report.save(top)

    def _write_refseq(self, top: Path, brotli_level: int):
        """ Write the reference sequence to a file. """
        refseq_file = RefseqIO(sample=self.sample,
                               ref=self.ref,
                               refseq=self.seq)
        _, checksum = refseq_file.save(top, brotli_level)
        return checksum

    def _generate_batches(self, *,
                          top: Path,
                          keep_tmp: bool,
                          min_mapq: int,
                          phred_enc: int,
                          min_phred: int,
                          ambindel: bool,
                          overhangs: bool,
                          clip_end5: int,
                          clip_end3: int,
                          brotli_level: int,
                          n_procs: int):
        """ Compute a relation vector for every record in a XAM file,
        split among one or more batches. For each batch, write a matrix
        of the vectors to one batch file, and compute its checksum. """
        logger.info(f"Began {self}")
        try:
            # Collect the keyword arguments.
            kwargs = dict(xam_view=self._xam,
                          top=top,
                          refseq=self.seq,
                          min_mapq=min_mapq,
                          min_qual=encode_phred(min_phred, phred_enc),
                          ambindel=ambindel,
                          overhangs=overhangs,
                          clip_end5=clip_end5,
                          clip_end3=clip_end3,
                          brotli_level=brotli_level)
            # Generate and write relation vectors for each batch.
            results = dispatch(generate_batch,
                               n_procs,
                               parallel=True,
                               pass_n_procs=False,
                               args=as_list_of_tuples(self._xam.indexes),
                               kwargs=kwargs)
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
            if not keep_tmp:
                # Delete the temporary SAM file before exiting.
                self._xam.delete_tmp_sam()

    def write(self, *,
              out_dir: Path,
              release_dir: Path,
              min_mapq: int,
              min_reads: int,
              brotli_level: int,
              force: bool,
              overhangs: bool,
              min_phred: int,
              phred_enc: int,
              ambindel: bool,
              clip_end5: int,
              clip_end3: int,
              **kwargs):
        """ Compute a relation vector for every record in a BAM file,
        write the vectors into one or more batch files, compute their
        checksums, and write a report summarizing the results. """
        report_file = RelateReport.build_path(top=out_dir,
                                              sample=self.sample,
                                              ref=self.ref)
        if need_write(report_file, force):
            began = datetime.now()
            # Determine if there are enough reads.
            if self.num_reads < min_reads:
                raise ValueError(f"Insufficient reads in {self._xam}: "
                                 f"{self.num_reads} < {min_reads}")
            # Write the reference sequence to a file.
            refseq_checksum = self._write_refseq(release_dir, brotli_level)
            # Compute relation vectors and time how long it takes.
            (nreads,
             nbats,
             checks) = self._generate_batches(top=release_dir,
                                              brotli_level=brotli_level,
                                              min_mapq=min_mapq,
                                              overhangs=overhangs,
                                              min_phred=min_phred,
                                              phred_enc=phred_enc,
                                              ambindel=ambindel,
                                              clip_end5=clip_end5,
                                              clip_end3=clip_end3,
                                              **kwargs)
            ended = datetime.now()
            # Write a report of the relation step.
            report_saved = self._write_report(top=release_dir,
                                              min_mapq=min_mapq,
                                              min_phred=min_phred,
                                              phred_enc=phred_enc,
                                              overhangs=overhangs,
                                              ambindel=ambindel,
                                              clip_end5=clip_end5,
                                              clip_end3=clip_end3,
                                              min_reads=min_reads,
                                              n_reads_xam=self.num_reads,
                                              n_reads_rel=nreads,
                                              n_batches=nbats,
                                              checksums=checks,
                                              refseq_checksum=refseq_checksum,
                                              began=began,
                                              ended=ended)
            release_to_out(out_dir, release_dir, report_saved.parent)
        return report_file

    def __str__(self):
        return f"Relate {self._xam}"


def write_one(xam_file: Path, *,
              fasta: Path,
              tmp_dir: Path,
              batch_size: int,
              n_procs: int,
              **kwargs):
    """ Write the batches of relation vectors for one XAM file. """
    release_dir, working_dir = get_release_working_dirs(tmp_dir)
    ref = path.parse(xam_file, *path.XAM_SEGS)[path.REF]
    writer = RelationWriter(XamViewer(xam_file,
                                      working_dir,
                                      batch_size,
                                      n_procs=n_procs),
                            get_fasta_seq(fasta, DNA, ref))
    return writer.write(**kwargs, n_procs=n_procs, release_dir=release_dir)


def write_all(xam_files: Iterable[Path],
              max_procs: int,
              parallel: bool,
              **kwargs):
    """  """
    return dispatch(write_one,
                    max_procs,
                    parallel,
                    args=as_list_of_tuples(path.deduplicate(xam_files)),
                    kwargs=kwargs)

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
