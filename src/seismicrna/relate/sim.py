from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .batch import format_read_name
from .io import ReadNamesBatchIO, RelateBatchIO, RefseqIO
from .report import RelateReport
from ..core import path
from ..core.seq import DNA
from ..core.tmp import release_to_out
from ..core.write import need_write

rng = np.random.default_rng()


def _update_checksums(current_checksums: dict[str, list[str]],
                      new_checksums: dict[str, list[str]]):
    for key, values in new_checksums.items():
        try:
            current_checksums[key].extend(values)
        except KeyError:
            current_checksums[key] = values


def simulate_batch(sample: str,
                   branches: dict[str, str],
                   ref: str,
                   batch: int,
                   write_read_names: bool,
                   pmut: pd.DataFrame,
                   uniq_end5s: np.ndarray,
                   uniq_end3s: np.ndarray,
                   pends: np.ndarray,
                   paired: bool,
                   read_length: int,
                   p_rev: float,
                   min_mut_gap: int,
                   num_reads: int,
                   formatter: Callable[[int, int], str] = format_read_name):
    """ Simulate a pair of RelateBatchIO and ReadNamesBatchIO. """
    relate_batch = RelateBatchIO.simulate(sample=sample,
                                          branches=branches,
                                          ref=ref,
                                          batch=batch,
                                          pmut=pmut,
                                          uniq_end5s=uniq_end5s,
                                          uniq_end3s=uniq_end3s,
                                          pends=pends,
                                          paired=paired,
                                          read_length=read_length,
                                          p_rev=p_rev,
                                          min_mut_gap=min_mut_gap,
                                          num_reads=num_reads)
    if write_read_names:
        name_batch = ReadNamesBatchIO.simulate(sample=sample,
                                               branches=branches,
                                               ref=ref,
                                               batch=batch,
                                               num_reads=relate_batch.num_reads,
                                               formatter=formatter)
    else:
        name_batch = None
    return relate_batch, name_batch


def simulate_cluster(first_batch: int,
                     batch_size: int,
                     num_reads: int,
                     **kwargs):
    """ Simulate all batches for one cluster. """
    # Determine the numbers of batches and reads per batch.
    num_full_batches, last_batch_size = divmod(int(num_reads), int(batch_size))
    last_batch = first_batch + num_full_batches
    # Simulate every full-size batch.
    for batch in range(first_batch, last_batch):
        yield simulate_batch(batch=batch,
                             num_reads=batch_size,
                             **kwargs)
    # Simulate the last batch, which may have fewer reads.
    if last_batch_size > 0:
        yield simulate_batch(batch=last_batch,
                             num_reads=last_batch_size,
                             **kwargs)


def simulate_batches(batch_size: int,
                     pmut: pd.DataFrame,
                     pclust: pd.Series,
                     num_reads: int,
                     **kwargs):
    # Simulate the number of reads per cluster.
    num_reads_clusters = pd.Series(rng.multinomial(num_reads, pclust),
                                   pclust.index)
    # Simulate batches for each cluster.
    first_batch = 0
    for k, clust in pclust.index:
        num_reads_cluster = num_reads_clusters.at[(k, clust)]
        pmut_cluster = pmut.loc[:, (slice(None), k, clust)]
        for batch in simulate_cluster(first_batch,
                                      batch_size,
                                      num_reads_cluster,
                                      pmut=pmut_cluster,
                                      **kwargs):
            yield batch
            first_batch += 1


def simulate_relate(*,
                    out_dir: Path,
                    tmp_dir: Path,
                    branch: str,
                    sample: str,
                    ref: str,
                    refseq: DNA,
                    batch_size: int,
                    num_reads: int,
                    write_read_names: bool,
                    pmut: pd.DataFrame,
                    uniq_end5s: np.ndarray,
                    uniq_end3s: np.ndarray,
                    pends: np.ndarray,
                    pclust: pd.Series,
                    brotli_level: int,
                    force: bool,
                    **kwargs):
    """ Simulate an entire relate step. """
    branches = path.add_branch(path.RELATE_STEP, branch, dict())
    report_file = RelateReport.build_path({path.TOP: out_dir,
                                           path.SAMPLE: sample,
                                           path.BRANCHES: branches,
                                           path.REF: ref})
    if need_write(report_file, force):
        began = datetime.now()
        # Write the reference sequence to a file.
        refseq_file = RefseqIO(sample=sample,
                               branches=branches,
                               ref=ref,
                               refseq=refseq)
        _, refseq_checksum = refseq_file.save(tmp_dir,
                                              brotli_level=brotli_level,
                                              force=True)
        # Simulate and write the batches.
        checksums = {RelateBatchIO.btype(): list()}
        if write_read_names:
            checksums[ReadNamesBatchIO.btype()] = list()
        read_count = 0
        for relate_batch, name_batch in simulate_batches(
                sample=sample,
                branches=branches,
                ref=ref,
                batch_size=batch_size,
                num_reads=num_reads,
                write_read_names=write_read_names,
                pmut=pmut,
                uniq_end5s=uniq_end5s,
                uniq_end3s=uniq_end3s,
                pends=pends,
                pclust=pclust,
                **kwargs
        ):
            _, relate_checksum = relate_batch.save(tmp_dir,
                                                   brotli_level=brotli_level,
                                                   force=True)
            checksums[RelateBatchIO.btype()].append(relate_checksum)
            if write_read_names:
                assert isinstance(name_batch, ReadNamesBatchIO)
                _, name_checksum = name_batch.save(tmp_dir,
                                                   brotli_level=brotli_level,
                                                   force=True)
                checksums[ReadNamesBatchIO.btype()].append(name_checksum)
            else:
                assert name_batch is None
            read_count += relate_batch.num_reads
        if write_read_names:
            assert ReadNamesBatchIO.btype() in checksums
            assert (len(checksums[ReadNamesBatchIO.btype()])
                    == len(checksums[RelateBatchIO.btype()]))
        else:
            assert ReadNamesBatchIO.btype() not in checksums
        ended = datetime.now()
        # Write the report.
        report = RelateReport(sample=sample,
                              branches=branches,
                              ref=ref,
                              min_mapq=0,
                              phred_enc=0,
                              min_phred=0,
                              insert3=False,
                              ambindel=False,
                              clip_end5=0,
                              clip_end3=0,
                              min_reads=0,
                              n_reads_xam=0,
                              n_reads_rel=read_count,
                              n_batches=len(checksums[RelateBatchIO.btype()]),
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=ended,
                              overhangs=True)
        report_saved = report.save(tmp_dir, force=True)
        release_to_out(out_dir, tmp_dir, report_saved.parent)
    return report_file
