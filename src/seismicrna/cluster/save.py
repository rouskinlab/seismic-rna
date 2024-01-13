import pandas as pd

from .compare import RunOrderResults
from .data import ClusterMutsDataset
from .io import ClustBatchIO
from ..mask.data import MaskMutsDataset


def write_batches(dataset: MaskMutsDataset,
                  orders: list[RunOrderResults],
                  brotli_level: int):
    """ Write the cluster memberships to batch files. """
    checksums = list()
    for batch_num in dataset.batch_nums:
        resps = pd.concat((runs.best.get_resps(batch_num) for runs in orders),
                          axis=1)
        batch = ClustBatchIO(sample=dataset.sample,
                             ref=dataset.ref,
                             sect=dataset.sect,
                             batch=batch_num,
                             resps=resps)
        _, checksum = batch.save(top=dataset.top,
                                 brotli_level=brotli_level,
                                 force=True)
        checksums.append(checksum)
    return checksums


def update_batches(dataset: ClusterMutsDataset,
                   orders: list[RunOrderResults],
                   brotli_level: int):
    """ Update the cluster memberships in batches. """
    checksums = list()
    for batch in dataset.iter_batches():
        # Merge the original responsibilities with the new ones.
        resps = pd.concat([batch.resps] + [runs.best.get_resps(batch.batch)
                                           for runs in orders],
                          axis=1,
                          verify_integrity=True)
        batch = ClustBatchIO(sample=dataset.sample,
                             ref=dataset.ref,
                             sect=dataset.sect,
                             batch=batch.batch,
                             resps=resps)
        _, checksum = batch.save(top=dataset.top,
                                 brotli_level=brotli_level,
                                 force=True)
        checksums.append(checksum)
    return checksums
