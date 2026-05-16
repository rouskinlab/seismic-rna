from datetime import datetime
from pathlib import Path

import numpy as np

from .mm import iter_mm_file
from ..core import path
from ..core.logs import logger
from ..core.rel.code import DELET, INS_3, INS_5, SUB_A, SUB_C, SUB_G, SUB_T, SUB_N
from ..core.seq import DNA
from ..core.tmp import get_release_working_dirs, release_to_out
from ..core.write import need_write
from ..relate.io import from_reads, ReadNamesBatchIO, RelateBatchIO, RefseqIO
from ..relate.report import RelateReport

# All substitution codes except the one that would match the reference base.
_SUB_FOR_BASE: dict[str, int] = {
    'A': SUB_C | SUB_G | SUB_T,
    'C': SUB_A | SUB_G | SUB_T,
    'G': SUB_A | SUB_C | SUB_T,
    'T': SUB_A | SUB_C | SUB_G,
    'N': SUB_N
}


def _build_mut_codes(refseq: DNA, insert3: bool) -> np.ndarray:
    """ Build a 1-indexed array mapping each reference position to the
    mutation code to assign when that position is mutated in an MM read.

    The code combines DELET, one insertion side (INS_3 or INS_5 per
    insert3), and all substitutions except the one matching the reference
    base (which would be a match, not a mutation).
    """
    ins_code = INS_3 if insert3 else INS_5
    codes = np.zeros(len(refseq) + 1, dtype=np.uint8)
    for i, base in enumerate(str(refseq), start=1):
        codes[i] = DELET | ins_code | _SUB_FOR_BASE.get(base, SUB_N)
    return codes


def _iter_batch_reads(reads: list[tuple[int, int, list[int]]],
                      batch_start: int,
                      batch_size: int,
                      mut_codes: np.ndarray):
    """ Yield (name, ((end5s, end3s), poss)) for one batch of MM reads.

    Coordinates are converted from 0-based (MM) to 1-based (seismic-rna).
    """
    for i, (start, end, mut_positions) in enumerate(
            reads[batch_start: batch_start + batch_size]):
        poss = {pos + 1: int(mut_codes[pos + 1]) for pos in mut_positions}
        yield (f"read_{batch_start + i}",
               (([start + 1], [end + 1]), poss))


def _import_one_ref(ref_id: str,
                    refseq: DNA,
                    reads: list[tuple[int, int, list[int]]],
                    *,
                    sample: str,
                    branches: dict[str, str],
                    out_dir: Path,
                    release_dir: Path,
                    min_reads: int,
                    batch_size: int,
                    insert3: bool,
                    write_read_names: bool,
                    brotli_level: int,
                    force: bool) -> Path | None:
    """ Convert one MM transcript block into relate batches and a report. """
    report_path = RelateReport.build_path(
        {path.TOP: out_dir,
         path.SAMPLE: sample,
         path.BRANCHES: branches,
         path.REF: ref_id}
    )
    if not need_write(report_path, force):
        return report_path.parent

    n_reads_xam = len(reads)
    if n_reads_xam < min_reads:
        logger.warning(
            f"Skipping reference {repr(ref_id)}: "
            f"{n_reads_xam} reads < min_reads ({min_reads})"
        )
        return None

    logger.routine(f"Began importing MM reference {repr(ref_id)}")
    began = datetime.now()

    # Write the reference sequence.
    refseq_io = RefseqIO(branches=branches,
                         sample=sample,
                         ref=ref_id,
                         refseq=refseq)
    _, refseq_checksum = refseq_io.save(release_dir, brotli_level)

    # Build per-position mutation codes.
    mut_codes = _build_mut_codes(refseq, insert3)

    # Write batches.
    relate_checksums: list[str] = []
    name_checksums: list[str | None] = []
    n_reads_rel = 0
    batch_num = 0
    batch_start = 0
    while batch_start < n_reads_xam:
        relate_batch, name_batch = from_reads(
            _iter_batch_reads(reads, batch_start, batch_size, mut_codes),
            sample=sample,
            branches=branches,
            ref=ref_id,
            refseq=refseq,
            batch=batch_num,
            write_read_names=write_read_names,
        )
        _, rc = relate_batch.save(release_dir, brotli_level)
        relate_checksums.append(rc)
        n_reads_rel += relate_batch.num_reads
        if write_read_names:
            assert isinstance(name_batch, ReadNamesBatchIO)
            _, nc = name_batch.save(release_dir, brotli_level)
            name_checksums.append(nc)
        else:
            name_checksums.append(None)
        batch_num += 1
        batch_start += batch_size

    checksums = {RelateBatchIO.btype(): relate_checksums}
    if write_read_names:
        checksums[ReadNamesBatchIO.btype()] = name_checksums

    ended = datetime.now()

    # Write the report.
    report = RelateReport(
        sample=sample,
        ref=ref_id,
        branches=branches,
        min_mapq=0,
        phred_enc=0,
        min_phred=0,
        insert3=insert3,
        ambindel=False,
        overhangs=True,
        clip_end5=0,
        clip_end3=0,
        min_reads=min_reads,
        n_reads_xam=n_reads_xam,
        n_reads_rel=n_reads_rel,
        n_batches=len(relate_checksums),
        checksums=checksums,
        refseq_checksum=refseq_checksum,
        began=began,
        ended=ended,
    )
    report_saved = report.save(release_dir)
    release_to_out(out_dir, release_dir, report_saved.parent)
    logger.routine(f"Ended importing MM reference {repr(ref_id)}")
    return report_path.parent


def import_mm(mm_path: Path, *,
              sample: str,
              out_dir: Path,
              tmp_dir: Path,
              branch: str,
              min_reads: int,
              batch_size: int,
              insert3: bool,
              write_read_names: bool,
              brotli_level: int,
              force: bool) -> list[Path]:
    """ Convert all transcript blocks in one MM file into relate outputs.

    Returns a list of output directories (one per reference that was
    successfully imported).
    """
    release_dir, _ = get_release_working_dirs(tmp_dir)
    branches = path.add_branch(path.RELATE_STEP, branch, {})
    results = []
    for ref_id, refseq, reads in iter_mm_file(mm_path):
        result = _import_one_ref(
            ref_id, refseq, reads,
            sample=sample,
            branches=branches,
            out_dir=out_dir,
            release_dir=release_dir,
            min_reads=min_reads,
            batch_size=batch_size,
            insert3=insert3,
            write_read_names=write_read_names,
            brotli_level=brotli_level,
            force=force,
        )
        if result is not None:
            results.append(result)
    return results
