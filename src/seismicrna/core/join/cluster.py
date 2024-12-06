from pathlib import Path

import pandas as pd

from ..header import ClustHeader, parse_header


def parse_join_clusts_file(file: str | Path):
    """ Parse a file of joined clusters. """
    n_cols = len(ClustHeader.level_names())
    clusts_df = pd.read_csv(file, index_col=list(range(n_cols)))
    header = parse_header(clusts_df.index)
    # Verify the index: use type() is not rather than isinstance() so
    # that subclasses of ClustHeader will yield False, not True.
    if type(header) is not ClustHeader:
        raise TypeError(f"Expected first {n_cols} of {file} to be a valid "
                        f"{ClustHeader.__name__}, but got {header}")
    # Rearrange the DataFrame into a dict.
    clusts_dict = {reg: {k: dict() for k in header.ks}
                   for reg in clusts_df.columns}
    for reg, clusts in clusts_df.items():
        for (k, clust), reg_clust in clusts.items():
            if not 1 <= reg_clust <= k:
                raise ValueError(f"Region {repr(reg)} k {k} got a "
                                 f"cluster number out of range: {reg_clust}")
            if reg_clust in clusts_dict[reg][k].values():
                raise ValueError(f"Region {repr(reg)} k={k} got a "
                                 f"repeated cluster number: {reg_clust}")
            clusts_dict[reg][k][clust] = reg_clust
    return clusts_dict
