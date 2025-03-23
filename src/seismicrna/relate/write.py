from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import Iterable

from .io import from_reads, ReadNamesBatchIO, RelateBatchIO, RefseqIO
from .report import RelateReport
from .sam import XamViewer
from .table import RelateCountTabulator
from ..core import path
from ..core.logs import logger
from ..core.ngs import encode_phred
from ..core.seq import DNA, Region, get_fasta_seq
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
        # Try to load the C extension module.
        try:
            from .cx.relate import RelateError, calc_rels_lines
        except ImportError:
            logger.warning(
                "Failed to import the C extension for the relate algorithm; "
                "defaulting to the Python implementation, which is much slower"
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
            logger.error(RelateError(f"Read {repr(name)}: {error}"))


def generate_batch(batch: int, *,
                   xam_view: XamViewer,
                   top: Path,
                   refseq: DNA,
                   brotli_level: int,
                   count_pos: bool,
                   count_read: bool,
                   write_read_names: bool,
                   **kwargs):
    """ Compute relationships for every SAM record in one batch. """
    logger.routine(f"Began calculating batch {batch} of {xam_view}")
    relate_batch, name_batch = from_reads(
        relate_records(xam_view.iter_records(batch),
                       ref=xam_view.ref,
                       refseq=str(refseq),
                       **kwargs),
        branches=xam_view.branches,
        sample=xam_view.sample,
        ref=xam_view.ref,
        refseq=refseq,
        batch=batch,
        write_read_names=write_read_names
    )
    logger.routine(f"Ended calculating batch {batch} of {xam_view}")
    # Save the RelateBatchIO instance, which has as few attributes as
    # possible to make the file as small as possible.
    _, relate_checksum = relate_batch.save(top, brotli_level)
    if write_read_names:
        assert isinstance(name_batch, ReadNamesBatchIO)
        _, name_checksum = name_batch.save(top, brotli_level)
    else:
        assert name_batch is None
        name_checksum = None
    # Generate a RelateRegionMutsBatch in order to count the mutations,
    # which the RelateBatchIO instance cannot do.
    relate_region_batch = relate_batch.to_region_batch(Region(xam_view.ref,
                                                              refseq))
    return (relate_region_batch.count_all(all_patterns(),
                                          count_pos=count_pos,
                                          count_read=count_read,
                                          count_ends=False),
            relate_checksum,
            name_checksum)


class RelationWriter(object):
    """ Compute and write relationships for all reads from one sample
    aligned to one reference sequence. """

    def __init__(self, xam_view: XamViewer, fasta_file: str | Path):
        self._xam = xam_view
        self._fasta = fasta_file

    @property
    def sample(self):
        return self._xam.sample

    @property
    def ref(self):
        return self._xam.ref

    @cached_property
    def refseq(self):
        return get_fasta_seq(self._fasta, DNA, self.ref)

    @property
    def num_reads(self):
        return self._xam.n_reads

    @property
    def branches(self):
        return self._xam.branches

    def _write_report(self, *, top: Path, **kwargs):
        report = RelateReport(sample=self.sample,
                              ref=self.ref,
                              branches=self.branches,
                              n_reads_xam=self.num_reads,
                              **kwargs)
        return report.save(top)

    def _write_refseq(self, top: Path, brotli_level: int):
        """ Write the reference sequence to a file. """
        refseq_file = RefseqIO(branches=self.branches,
                               sample=self.sample,
                               ref=self.ref,
                               refseq=self.refseq)
        _, checksum = refseq_file.save(top, brotli_level)
        return checksum

    def _generate_batches(self, *,
                          top: Path,
                          write_read_names: bool,
                          keep_tmp: bool,
                          phred_enc: int,
                          min_phred: int,
                          num_cpus: int,
                          **kwargs):
        """ Compute the relationships for every read in a XAM file,
        split among one or more batches. """
        logger.routine(f"Began generating batches for {self._xam}")
        try:
            kwargs = dict(xam_view=self._xam,
                          top=top,
                          refseq=self.refseq,
                          min_qual=ord(encode_phred(min_phred, phred_enc)),
                          write_read_names=write_read_names,
                          **kwargs)
            results = dispatch(generate_batch,
                               num_cpus=num_cpus,
                               pass_num_cpus=False,
                               as_list=True,
                               ordered=True,
                               raise_on_error=True,
                               args=as_list_of_tuples(self._xam.indexes),
                               kwargs=kwargs)
            if results:
                (batch_counts,
                 relate_checksums,
                 name_checksums) = map(list, zip(*results, strict=True))
            else:
                batch_counts = list()
                relate_checksums = list()
                name_checksums = list()
            assert (len(self._xam.indexes)
                    == len(batch_counts)
                    == len(relate_checksums)
                    == len(name_checksums))
            checksums = {RelateBatchIO.btype(): relate_checksums}
            if write_read_names:
                assert all(name_checksums)
                checksums[ReadNamesBatchIO.btype()] = name_checksums
            else:
                assert not any(name_checksums)
            logger.routine(f"Ended generating batches for {self._xam}")
            return batch_counts, checksums
        finally:
            if not keep_tmp:
                # Delete the temporary SAM file (which can be massive)
                # before exiting.
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
              num_cpus: int,
              **kwargs):
        """ Compute relationships for every record in a XAM file. """
        report_file = RelateReport.build_path({path.TOP: out_dir,
                                               path.SAMPLE: self.sample,
                                               path.BRANCHES: self.branches,
                                               path.REF: self.ref})
        if need_write(report_file, force):
            began = datetime.now()
            # Determine if there are enough reads.
            if self.num_reads < min_reads:
                raise ValueError(f"Insufficient reads in {self._xam}: "
                                 f"{self.num_reads} < {min_reads}")
            # Write the reference sequence to a file.
            refseq_checksum = self._write_refseq(release_dir, brotli_level)
            # Compute relationships and time how long it takes.
            batch_counts, checks = self._generate_batches(
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
                num_cpus=num_cpus,
                **kwargs
            )
            # Tabulate the data.
            tabulator = RelateCountTabulator(batch_counts=batch_counts,
                                             top=release_dir,
                                             branches=self.branches,
                                             sample=self.sample,
                                             ref=self.ref,
                                             refseq=self.refseq,
                                             count_pos=relate_pos_table,
                                             count_read=relate_read_table,
                                             validate=False)
            tabulator.write_tables(pos=relate_pos_table, read=relate_read_table)
            ended = datetime.now()
            # Write the report.
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
                n_reads_rel=tabulator.num_reads,
                n_batches=len(batch_counts),
                checksums=checks,
                refseq_checksum=refseq_checksum,
                began=began,
                ended=ended
            )
            release_to_out(out_dir, release_dir, report_saved.parent)
        return report_file.parent

    def __str__(self):
        return f"{type(self).__name__}:{self._xam}"


def relate_xam(xam_file: Path, *,
               fasta: Path,
               tmp_dir: Path,
               branch: str,
               batch_size: int,
               num_cpus: int,
               **kwargs):
    """ Write the batches of relationships for one XAM file. """
    release_dir, working_dir = get_release_working_dirs(tmp_dir)
    writer = RelationWriter(XamViewer(xam_file,
                                      working_dir,
                                      branch,
                                      batch_size,
                                      num_cpus=num_cpus),
                            fasta)
    return writer.write(**kwargs, num_cpus=num_cpus, release_dir=release_dir)
