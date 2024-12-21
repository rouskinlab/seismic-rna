from datetime import datetime
from pathlib import Path
from typing import Iterable

from .io import from_reads, ReadNamesBatchIO, RelateBatchIO
from .report import RelateReport
from .sam import XamViewer
from .table import RelateCountTabulator
from ..core import path
from ..core.io import RefseqIO
from ..core.logs import logger
from ..core.ngs import encode_phred
from ..core.seq import DNA, get_fasta_seq
from ..core.table import all_patterns
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.write import need_write


def relate_records(records: Iterable[tuple[str, str, str]],
                   ref: str,
                   refseq: str,
                   min_mapq: int,
                   min_qual: int,
                   insert3: bool,
                   ambindel: bool,
                   overhangs: bool,
                   clip_end5: int,
                   clip_end3: int,
                   relate_cx: bool):
    # Load the module.
    if relate_cx:
        # Load the C extension module.
        try:
            from .cx.relate import RelateError, calc_rels_lines
        except ImportError:
            logger.warning(
                "Failed to import the C extension for the relate algorithm; "
                "defaulting to the Python version, which is much slower"
            )
            from .py.relate import RelateError, calc_rels_lines
    else:
        # Load the Python module.
        from .py.relate import RelateError, calc_rels_lines
    # Process the records.
    for name, line1, line2 in records:
        try:
            yield name, calc_rels_lines(line1,
                                        line2,
                                        ref,
                                        refseq,
                                        min_mapq,
                                        min_qual,
                                        insert3,
                                        ambindel,
                                        overhangs,
                                        clip_end5,
                                        clip_end3)
        except RelateError as error:
            logger.error(error)


def generate_batch(batch: int, *,
                   xam_view: XamViewer,
                   top: Path,
                   refseq: DNA,
                   brotli_level: int,
                   count_pos: bool,
                   count_read: bool,
                   **kwargs):
    """ Compute relation vectors for every SAM record in one batch,
    write the vectors to a batch file, and return its MD5 checksum
    and the number of vectors. """
    logger.routine("Began computing read-reference relationships "
                   f"for batch {batch} of {xam_view}")
    names, relvecs = from_reads(relate_records(xam_view.iter_records(batch),
                                               ref=xam_view.ref,
                                               refseq=str(refseq),
                                               **kwargs),
                                xam_view.sample,
                                xam_view.ref,
                                refseq,
                                batch)
    logger.routine("Ended computing read-reference relationships "
                   f"for batch {batch} of {xam_view}")
    _, relv_check = relvecs.save(top, brotli_level)
    _, name_check = names.save(top, brotli_level)
    return (relvecs.count_all(all_patterns(),
                              count_pos=count_pos,
                              count_read=count_read,
                              count_ends=False),
            relv_check,
            name_check)


class RelationWriter(object):
    """ Compute and write relationships for all reads from one sample
    aligned to one reference sequence. """

    def __init__(self, xam_view: XamViewer, refseq: DNA):
        self._xam = xam_view
        self.refseq = refseq

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
                               refseq=self.refseq)
        _, checksum = refseq_file.save(top, brotli_level)
        return checksum

    def _generate_batches(self, *,
                          top: Path,
                          keep_tmp: bool,
                          phred_enc: int,
                          min_phred: int,
                          n_procs: int,
                          **kwargs):
        """ Compute a relation vector for every record in a XAM file,
        split among one or more batches. For each batch, write a matrix
        of the vectors to one batch file, and compute its checksum. """
        logger.routine(f"Began generating batches for {self}")
        try:
            # Collect the keyword arguments.
            kwargs = dict(xam_view=self._xam,
                          top=top,
                          refseq=self.refseq,
                          min_qual=ord(encode_phred(min_phred, phred_enc)),
                          **kwargs)
            # Generate and write relation vectors for each batch.
            num_batches = len(self._xam.indexes)
            results = dispatch(generate_batch,
                               n_procs,
                               pass_n_procs=False,
                               args=as_list_of_tuples(self._xam.indexes),
                               kwargs=kwargs)
            if results:
                batch_counts, relv_checks, name_checks = map(list,
                                                             zip(*results,
                                                                 strict=True))
            else:
                batch_counts = list()
                relv_checks = list()
                name_checks = list()
            if len(batch_counts) != num_batches:
                raise ValueError(f"Expected {num_batches} batch(es), "
                                 f"but generated {len(batch_counts)}")
            checksums = {RelateBatchIO.btype(): relv_checks,
                         ReadNamesBatchIO.btype(): name_checks}
            logger.routine(f"Ended generating batches {num_batches} for {self}")
            return batch_counts, num_batches, checksums
        finally:
            if not keep_tmp:
                # Delete the temporary SAM file before exiting.
                self._xam.delete_tmp_sam()

    def write(self, *,
              out_dir: Path,
              release_dir: Path,
              min_mapq: int,
              min_reads: int,
              min_phred: int,
              phred_enc: int,
              insert3: bool,
              ambindel: bool,
              overhangs: bool,
              clip_end5: int,
              clip_end3: int,
              relate_pos_table: bool,
              relate_read_table: bool,
              brotli_level: int,
              force: bool,
              n_procs: int,
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
            # Compute relationships and time how long it takes.
            (batch_counts,
             n_batches,
             checks) = self._generate_batches(
                top=release_dir,
                brotli_level=brotli_level,
                min_mapq=min_mapq,
                min_phred=min_phred,
                phred_enc=phred_enc,
                insert3=insert3,
                ambindel=ambindel,
                overhangs=overhangs,
                clip_end5=clip_end5,
                clip_end3=clip_end3,
                count_pos=relate_pos_table,
                count_read=relate_read_table,
                n_procs=n_procs,
                **kwargs
            )
            # Tabulate the data.
            tabulator = RelateCountTabulator(batch_counts=batch_counts,
                                             top=release_dir,
                                             sample=self.sample,
                                             ref=self.ref,
                                             refseq=self.refseq,
                                             count_pos=relate_pos_table,
                                             count_read=relate_read_table,
                                             validate=False)
            tabulator.write_tables(pos=relate_pos_table, read=relate_read_table)
            ended = datetime.now()
            # Write a report of the relation step.
            report_saved = self._write_report(
                top=release_dir,
                min_mapq=min_mapq,
                min_phred=min_phred,
                phred_enc=phred_enc,
                insert3=insert3,
                ambindel=ambindel,
                overhangs=overhangs,
                clip_end5=clip_end5,
                clip_end3=clip_end3,
                min_reads=min_reads,
                n_reads_xam=self.num_reads,
                n_reads_rel=tabulator.num_reads,
                n_batches=n_batches,
                checksums=checks,
                refseq_checksum=refseq_checksum,
                began=began,
                ended=ended
            )
            release_to_out(out_dir, release_dir, report_saved.parent)
        return report_file.parent

    def __str__(self):
        return f"Relate {self._xam}"


def write_one(xam_file: Path, *,
              fasta: Path,
              tmp_dir: Path,
              batch_size: int,
              n_procs: int,
              **kwargs):
    """ Write the batches of relationships for one XAM file. """
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
              **kwargs):
    """ Write the batches of relationships for all XAM files. """
    return dispatch(write_one,
                    max_procs,
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
