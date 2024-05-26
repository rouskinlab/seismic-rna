from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .batch import format_read_name
from .io import QnamesBatchIO, RelateBatchIO
from .report import RelateReport
from ..core.io import RefseqIO
from ..core.seq import DNA
from ..core.write import need_write

rng = np.random.default_rng()


def _update_checksums(current_checksums: dict[str, list[str]],
                      new_checksums: dict[str, list[str]]):
    for key, values in new_checksums.items():
        try:
            current_checksums[key].extend(values)
        except KeyError:
            current_checksums[key] = values


def simulate_batch(out_dir: Path,
                   brotli_level: int,
                   sample: str,
                   ref: str,
                   batch: int,
                   pmut: pd.DataFrame,
                   uniq_end5s: np.ndarray,
                   uniq_end3s: np.ndarray,
                   pends: np.ndarray,
                   num_reads: int,
                   formatter: Callable[[int, int], str] = format_read_name):
    """ Simulate a pair of a QnamesBatchIO and RelateBatchIO, save them,
    and return their checksums. """
    qnames_batch = QnamesBatchIO.simulate(sample=sample,
                                          ref=ref,
                                          batch=batch,
                                          num_reads=num_reads,
                                          formatter=formatter)
    relate_batch = RelateBatchIO.simulate(sample=sample,
                                          ref=ref,
                                          batch=batch,
                                          pmut=pmut,
                                          uniq_end5s=uniq_end5s,
                                          uniq_end3s=uniq_end3s,
                                          pends=pends,
                                          num_reads=num_reads)
    _, relate_check = relate_batch.save(out_dir, brotli_level, force=True)
    _, qnames_check = qnames_batch.save(out_dir, brotli_level, force=True)
    return {RelateBatchIO.btype(): [relate_check],
            QnamesBatchIO.btype(): [qnames_check]}


def simulate_cluster(first_batch: int,
                     batch_size: int,
                     num_reads: int,
                     **kwargs):
    """ Simulate all batches for one cluster. """
    checksums = dict()
    # Determine the numbers of batches and reads per batch.
    num_full_batches, last_batch_size = divmod(int(num_reads), int(batch_size))
    last_batch = first_batch + num_full_batches
    # Simulate every full-size batch.
    for batch in range(first_batch, last_batch):
        _update_checksums(checksums,
                          simulate_batch(batch=batch,
                                         num_reads=batch_size,
                                         **kwargs))
    # Simulate the last batch, which may have fewer reads.
    if last_batch_size > 0:
        _update_checksums(checksums,
                          simulate_batch(batch=last_batch,
                                         num_reads=last_batch_size,
                                         **kwargs))
    return checksums


def simulate_batches(batch_size: int,
                     pmut: list[pd.DataFrame],
                     pclust: np.ndarray,
                     num_reads: int,
                     **kwargs):
    # Simulate the number of reads per cluster.
    num_reads_clusters = rng.multinomial(num_reads, pclust)
    # Simulate batches for each cluster.
    first_batch = 0
    checksums = dict()
    for num_reads_cluster, pmut_cluster in zip(num_reads_clusters,
                                               pmut,
                                               strict=True):
        _update_checksums(checksums,
                          simulate_cluster(first_batch,
                                           batch_size,
                                           num_reads_cluster,
                                           pmut=pmut_cluster,
                                           **kwargs))
        first_batch = len(checksums[RelateBatchIO.btype()])
    return checksums


def simulate_relate(out_dir: Path,
                    sample: str,
                    ref: str,
                    refseq: DNA,
                    batch_size: int,
                    num_reads: int,
                    pmut: list[pd.DataFrame],
                    uniq_end5s: np.ndarray,
                    uniq_end3s: np.ndarray,
                    pends: np.ndarray,
                    pclust: np.ndarray,
                    brotli_level: int,
                    force: bool,
                    **kwargs):
    """ Simulate an entire relate step. """
    report_file = RelateReport.build_path(top=out_dir,
                                          sample=sample,
                                          ref=ref)
    if need_write(report_file, force):
        began = datetime.now()
        # Write the reference sequence to a file.
        refseq_file = RefseqIO(sample=sample, ref=ref, refseq=refseq)
        _, refseq_checksum = refseq_file.save(out_dir, brotli_level, force=True)
        # Simulate and write the batches.
        checksums = simulate_batches(out_dir=out_dir,
                                     sample=sample,
                                     ref=ref,
                                     batch_size=batch_size,
                                     num_reads=num_reads,
                                     pmut=pmut,
                                     uniq_end5s=uniq_end5s,
                                     uniq_end3s=uniq_end3s,
                                     pends=pends,
                                     pclust=pclust,
                                     brotli_level=brotli_level,
                                     **kwargs)
        ns_batches = sorted(set(map(len, checksums.values())))
        if len(ns_batches) != 1:
            raise ValueError("Expected exactly 1 number of batches, "
                             f"but got {ns_batches}")
        n_batches = ns_batches[0]
        ended = datetime.now()
        # Write the report.
        report = RelateReport(sample=sample,
                              ref=ref,
                              min_mapq=-1,
                              phred_enc=0,
                              min_phred=0,
                              ambindel=False,
                              min_reads=num_reads,
                              n_reads_xam=num_reads,
                              n_reads_rel=num_reads,
                              n_batches=n_batches,
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=ended,
                              overhangs=True)
        report.save(out_dir, force=True)
    return report_file
