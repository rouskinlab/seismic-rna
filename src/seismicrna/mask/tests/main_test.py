import unittest as ut
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

from seismicrna.core import path
from seismicrna.core.io.seq import RefseqIO
from seismicrna.core.seq.section import Section
from seismicrna.core.seq.xna import DNA
from seismicrna.relate.data import RelateDataset, ReadNamesDataset
from seismicrna.relate.io import RelateBatchIO, ReadNamesBatchIO
from seismicrna.relate.report import RelateReport

SAMPLE = "sample"
REF = "ref"
REF_SEQ = DNA("CGCAAATC")
SECTION = Section(REF, REF_SEQ)

# Read CCGAAATC
#    0 ========
#    1 ?=======
#    2 ==?=?===
#    3 ===D====
#    4 ====II==
#    5 ..===G=G
#    6 .===C===
#    7 .A==T===
#    8 .===.===
READ_NAMES = [f"Read{i}" for i in range(8)]
SINGLE_END5S = [[1], [1], [1], [1], [3], [2], [1], [2]]
SINGLE_END3S = [[8], [8], [8], [8], [8], [8], [7], [8]]
PAIRED_END5S = [[1, 2], [3, 1], [1, 4], [5, 1], [3, 6], [7, 2], [1, 8], [6, 2]]
PAIRED_END3S = [[1, 8], [8, 2], [3, 8], [8, 4], [5, 8], [8, 6], [7, 8], [8, 4]]
MUTS = {1: {209: [1]},
        2: {16: [7]},
        3: {177: [2]},
        4: {2: [3]},
        5: {5: [4], 32: [6], 128: [7], 225: [2]},
        6: {9: [4], 64: [5]},
        7: {},
        8: {64: [5]}}


def write_dataset(out_dir: Path,
                  paired: bool,
                  reads_per_sample: dict[str, int],
                  reads_per_batch: int):
    if sum(reads_per_sample.values()) != len(READ_NAMES):
        raise ValueError(f"reads_per_sample must sum to {len(READ_NAMES)}, "
                         f"but got {sum(reads_per_sample.values())}")
    read_num = 0
    for sample, num_reads in reads_per_sample.items():
        began = datetime.now()
        # Write the reference sequence.
        refseq = RefseqIO(sample=sample, ref=REF, refseq=REF_SEQ)
        _, refseq_checksum = refseq.save(out_dir)
        # Assign read numbers to each batch.
        read_nums = iter(range(read_num, read_num + num_reads))
        i_to_batches = list()
        # Map the read numbers (i) to the read numbers in the batch.
        while batch_read_nums := list(zip(range(reads_per_batch), read_nums)):
            i_to_batches.append({i: batch_i for batch_i, i in batch_read_nums})
        # Take the read numbers for this sample in groups equal to the
        # number of reads per batch.
        checksums = {ReadNamesBatchIO.btype(): list(),
                     RelateBatchIO.btype(): list()}
        for batch, i_to_batch in enumerate(i_to_batches):
            names = np.array([READ_NAMES[i] for i in i_to_batch])
            end5s = np.array([PAIRED_END5S[i] if paired else SINGLE_END5S[i]
                              for i in i_to_batch])
            end3s = np.array([PAIRED_END3S[i] if paired else SINGLE_END3S[i]
                              for i in i_to_batch])
            # When assembling the mutations for this batch, map the read
            # numbers to the corresponding numbers within the batch.
            muts = dict()
            for j, rels in MUTS.items():
                muts[j] = defaultdict(list)
                for rel, reads in rels.items():
                    for i in reads:
                        batch_i = i_to_batch.get(i)
                        if batch_i is not None:
                            muts[j][rel].append(batch_i)
            relate_batch = RelateBatchIO(
                sample=sample,
                section=SECTION,
                batch=batch,
                seg_end5s=end5s,
                seg_end3s=end3s,
                muts={j: {rel: np.array(reads)
                          for rel, reads in rels.items()}
                      for j, rels in muts.items()}
            )
            checksums[RelateBatchIO.btype()].append(
                relate_batch.save(out_dir)[1]
            )
            name_batch = ReadNamesBatchIO(sample=sample,
                                          ref=REF,
                                          batch=batch,
                                          names=names)
            checksums[ReadNamesBatchIO.btype()].append(
                name_batch.save(out_dir)[1]
            )
        report = RelateReport(top=out_dir,
                              sample=sample,
                              ref=REF,
                              min_mapq=0,
                              min_phred=0,
                              phred_enc=33,
                              overhangs=True,
                              ambindel=False,
                              clip_end5=0,
                              clip_end3=0,
                              min_reads=0,
                              n_reads_xam=num_reads,
                              n_reads_rel=num_reads,
                              n_batches=len(i_to_batches),
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=datetime.now())
        report.save(out_dir)


if __name__ == "__main__":
    write_dataset(Path("dummy"), True, {"sample1": 4, "sample2": 4}, 2)
    ut.main()
