from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import Iterable

from .io import from_reads, ReadNamesBatchIO, IDmutBatchIO, RefseqIO
from .report import IDmutReport
from .sam import SamFileViewer
from .table import IDmutCountTabulator
from ..core import path
from ..core.logs import logger, format_sample_reference_region
from ..core.ngs import encode_phred
from ..core.seq import DNA, Region, get_fasta_seq
from ..core.table import all_patterns
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.write import need_write


def idmut_records(
    records: Iterable[tuple[str, str, str]],
    ref: str,
    refseq: str,
    min_mapq: int,
    min_qual: int,
    insert3: bool,
    ambindel: bool,
    overhangs: bool,
    clip_end5: int,
    clip_end3: int,
    idmut_cx: bool,
):
    """
    Yield relationships for each SAM record in an iterable.

    Parameters
    ----------
    records: Iterable[tuple[str, str, str]]
        Iterable of (name, line1, line2) tuples from a SAM file, where
        `line2` is empty for single-end reads.
    ref: str
        Reference name expected in each SAM record.
    refseq: str
        Full reference sequence string.
    min_mapq: int
        Minimum acceptable mapping quality score.
    min_qual: int
        Minimum Phred quality score (as an integer) for base calls.
    insert3: bool
        Whether to mark insertions on the 3' flanking reference
        position (True) or the 5' position (False).
    ambindel: bool
        Whether to find and label ambiguous indel positions.
    overhangs: bool
        Whether to allow paired-end mates to overhang one another.
    clip_end5: int
        Number of bases to clip from the 5' end of each read.
    clip_end3: int
        Number of bases to clip from the 3' end of each read.
    idmut_cx: bool
        Whether to use the C extension for the idmut algorithm;
        falls back to the Python implementation if import fails.

    Yields
    ------
    tuple[str, tuple]
        Read name and the result of `id_muts_lines` for that read.
    """
    # Load the module.
    if idmut_cx:
        # Try to load the C extension module.
        try:
            from .cx.idmut import IDmutError, id_muts_lines
        except ImportError:
            logger.warning(
                "Failed to import the C extension for the idmut algorithm; "
                "defaulting to the Python implementation, which is much slower"
            )
            from .py.idmut import IDmutError, id_muts_lines
    else:
        # Load the Python module.
        from .py.idmut import IDmutError, id_muts_lines
    # Process the records.
    for name, line1, line2 in records:
        try:
            yield (
                name,
                id_muts_lines(
                    line1,
                    line2,
                    ref,
                    refseq,
                    min_mapq,
                    min_qual,
                    insert3,
                    ambindel,
                    overhangs,
                    clip_end5,
                    clip_end3,
                ),
            )
        except IDmutError as error:
            logger.error(IDmutError(f"Read {repr(name)}: {error}"))


def generate_batch(
    batch: int,
    *,
    sam_view: SamFileViewer,
    top: Path,
    refseq: DNA,
    brotli_level: int,
    count_pos: bool,
    count_read: bool,
    write_read_names: bool,
    **kwargs,
):
    """Compute relationships for every SAM record in one batch."""
    with logger.debug.single_context("calculating batch {} of {}", batch, sam_view):
        idmut_batch, name_batch = from_reads(
            idmut_records(
                sam_view.iter_records(batch),
                ref=sam_view.ref,
                refseq=str(refseq),
                **kwargs,
            ),
            branches=sam_view.branches,
            sample=sam_view.sample,
            ref=sam_view.ref,
            refseq=refseq,
            batch=batch,
            write_read_names=write_read_names,
        )
    # Save the IDmutBatchIO instance, which has as few attributes as
    # possible to make the file as small as possible.
    _, idmut_checksum = idmut_batch.save(top, brotli_level)
    if write_read_names:
        assert isinstance(name_batch, ReadNamesBatchIO)
        _, name_checksum = name_batch.save(top, brotli_level)
    else:
        assert name_batch is None
        name_checksum = None
    # Generate an IDmutRegionMutsBatch in order to count the mutations,
    # which the IDmutBatchIO instance cannot do.
    idmut_region_batch = idmut_batch.to_region_batch(Region(sam_view.ref, refseq))
    return (
        idmut_region_batch.count_all(
            all_patterns(), count_pos=count_pos, count_read=count_read, count_ends=False
        ),
        idmut_checksum,
        name_checksum,
    )


class RelationWriter(object):
    """Compute and write relationships for all reads from one sample
    aligned to one reference sequence."""

    def __init__(self, sam_view: SamFileViewer, fasta_file: str | Path):
        """
        Initialize a RelationWriter.

        Parameters
        ----------
        sam_view: SamFileViewer
            Viewer for the input XAM file, providing access to reads
            and batch indexes.
        fasta_file: str | Path
            Path to the FASTA file containing the reference sequence.
        """
        self._xam = sam_view
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
        """
        Build and save an IDmutReport to disk.

        Parameters
        ----------
        top: Path
            Top-level directory in which to write the report file.
        **kwargs
            Additional fields forwarded to `IDmutReport.__init__`.

        Returns
        -------
        Path
            Path of the saved report file.
        """
        report = IDmutReport(
            sample=self.sample,
            ref=self.ref,
            branches=self.branches,
            n_reads_xam=self.num_reads,
            **kwargs,
        )
        return report.save(top)

    def _write_refseq(self, top: Path, brotli_level: int):
        """Write the reference sequence to a file."""
        refseq_file = RefseqIO(
            branches=self.branches, sample=self.sample, ref=self.ref, refseq=self.refseq
        )
        _, checksum = refseq_file.save(top, brotli_level)
        return checksum

    def _generate_batches(
        self,
        *,
        top: Path,
        write_read_names: bool,
        keep_tmp: bool,
        phred_enc: int,
        min_phred: int,
        num_cpus: int,
        **kwargs,
    ):
        """Compute the relationships for every read in a XAM file,
        split among one or more batches."""
        logger.debug("Began generating batches for {}", self._xam)
        try:
            kwargs = dict(
                sam_view=self._xam,
                top=top,
                refseq=self.refseq,
                min_qual=ord(encode_phred(min_phred, phred_enc)),
                write_read_names=write_read_names,
                **kwargs,
            )
            results = dispatch(
                generate_batch,
                num_cpus=num_cpus,
                pass_num_cpus=False,
                as_list=True,
                ordered=True,
                raise_on_error=True,
                args=as_list_of_tuples(self._xam.indexes),
                kwargs=kwargs,
            )
            if results:
                (batch_counts, idmut_checksums, name_checksums) = map(
                    list, zip(*results, strict=True)
                )
            else:
                batch_counts = list()
                idmut_checksums = list()
                name_checksums = list()
            assert (
                len(self._xam.indexes)
                == len(batch_counts)
                == len(idmut_checksums)
                == len(name_checksums)
            )
            checksums = {IDmutBatchIO.btype(): idmut_checksums}
            if write_read_names:
                assert all(name_checksums)
                checksums[ReadNamesBatchIO.btype()] = name_checksums
            else:
                assert not any(name_checksums)
            logger.debug("Ended generating batches for {}", self._xam)
            return batch_counts, checksums
        finally:
            if not keep_tmp:
                # Delete the temporary SAM file (which can be massive)
                # before exiting.
                self._xam.delete_tmp_sam()

    def write(
        self,
        *,
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
        idmut_pos_table: bool,
        idmut_read_table: bool,
        brotli_level: int,
        force: bool,
        num_cpus: int,
        **kwargs,
    ):
        """Compute relationships for every record in a XAM file."""
        report_file = IDmutReport.build_path(
            {
                path.TOP: out_dir,
                path.SAMPLE: self.sample,
                path.BRANCHES: self.branches,
                path.REF: self.ref,
            }
        )
        if need_write(report_file, force):
            began = datetime.now()
            # Determine if there are enough reads.
            if self.num_reads < min_reads:
                raise ValueError(
                    f"Insufficient reads in {self._xam}: {self.num_reads} < {min_reads}"
                )
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
                count_pos=idmut_pos_table,
                count_read=idmut_read_table,
                num_cpus=num_cpus,
                **kwargs,
            )
            # Tabulate the data.
            tabulator = IDmutCountTabulator(
                batch_counts=batch_counts,
                top=release_dir,
                branches=self.branches,
                sample=self.sample,
                ref=self.ref,
                refseq=self.refseq,
                count_pos=idmut_pos_table,
                count_read=idmut_read_table,
                validate=False,
            )
            tabulator.write_tables(pos=idmut_pos_table, read=idmut_read_table)
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
                ended=ended,
            )
            release_to_out(out_dir, release_dir, report_saved.parent)
        return report_file.parent

    def __str__(self):
        srr = format_sample_reference_region(self.sample, self.ref)
        return f"{type(self).__name__} of {srr}"


def idmut_xam(
    xam_file: Path,
    *,
    fasta: Path,
    tmp_dir: Path,
    branch: str,
    batch_size: int,
    num_cpus: int,
    **kwargs,
):
    """Write the batches of relationships for one XAM file."""
    release_dir, working_dir = get_release_working_dirs(tmp_dir)
    writer = RelationWriter(
        SamFileViewer(xam_file, working_dir, branch, batch_size, num_cpus=num_cpus), fasta
    )
    return writer.write(**kwargs, num_cpus=num_cpus, release_dir=release_dir)
