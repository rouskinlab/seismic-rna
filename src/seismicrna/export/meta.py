from logging import getLogger
from pathlib import Path
from typing import Any

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


def combine_metadata(special_metadata: dict[str, Any],
                     parsed_metadata: dict[Any, dict],
                     item: Any,
                     what: str = "item"):
    try:
        item_metadata = parsed_metadata[item]
    except KeyError:
        logger.debug(f"No metadata were given for {what} {repr(item)}")
        return special_metadata
    # Check for any keys in the parsed metadata that match those in the
    # special metadata.
    for field in set(special_metadata) & set(item_metadata):
        if (s := special_metadata[field]) != (p := item_metadata[field]):
            raise ValueError(f"Metadata field {repr(field)} of {what} "
                             f"{repr(item)} is {repr(s)}, but was also given "
                             f"as {repr(p)} in the metadata file")
    return special_metadata | item_metadata

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
