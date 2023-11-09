from logging import getLogger
from pathlib import Path

import pandas as pd

logger = getLogger(__name__)

SAMPLE_INDEX = "Sample"
REFERENCE_INDEX = "Reference"


def _parse_metadata(file: Path, index_col: str):
    """ Parse a CSV file of metadata.

    Parameters
    ----------
    file: Path
        CSV file of metadata

    Returns
    -------
    dict
        Parsed metadata
    """
    metadata = dict()
    df = pd.read_csv(file, index_col=index_col)
    if df.columns.has_duplicates:
        raise ValueError(f"Got duplicate metadata fields in file {file}: "
                         f"{df.columns[df.columns.duplicated()]}")
    for item in df.index:
        if item in metadata:
            logger.warning(f"Metadata for {index_col.lower()} {repr(item)} "
                           f"were given multiple times in file {file}")
            continue
        metadata[item] = {k: v for k, v in df.loc[item].to_dict().items()
                          if not pd.isnull(v)}
    return metadata


def parse_samples_metadata(file: Path):
    """ Parse a CSV file of metadata for each sample.

    Parameters
    ----------
    file: Path
        CSV file of metadata for each sample
    
    Returns
    -------
    dict
        Parsed metadata for each sample
    """
    return _parse_metadata(file, SAMPLE_INDEX)


def parse_refs_metadata(file: Path):
    """ Parse a CSV file of metadata for each reference.

    Parameters
    ----------
    file: Path
        CSV file of metadata for each reference

    Returns
    -------
    dict
        Parsed metadata for each reference
    """
    return _parse_metadata(file, REFERENCE_INDEX)
