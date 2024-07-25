from pathlib import Path

import pandas as pd

from .compare import EMRunsK
from .io import ClusterBatchIO
from ..mask.data import MaskMutsDataset


def write_batches(dataset: MaskMutsDataset,
                  ks: list[EMRunsK],
                  brotli_level: int,
                  top: Path):
    """ Write the cluster memberships to batch files. """
    checksums = list()
    for batch_num in dataset.batch_nums:
        resps = pd.concat((runs.best.get_resps(batch_num) for runs in ks),
                          axis=1)
        batch = ClusterBatchIO(sample=dataset.sample,
                               ref=dataset.ref,
                               sect=dataset.sect,
                               batch=batch_num,
                               resps=resps)
        _, checksum = batch.save(top, brotli_level=brotli_level)
        checksums.append(checksum)
    return checksums
