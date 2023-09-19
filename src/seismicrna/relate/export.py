from logging import getLogger

import pandas as pd

from .human import humanize_relvec
from .seqpos import parse_seq


logger = getLogger(__name__)


def as_iter(vectors: pd.DataFrame, reference: bool = False):
    """ For each vector (row) in vectors, yield a string of its name and
    its mutations as human-readable text.

    Parameters
    ----------
    vectors: DataFrame
        The vectors to display as a DataFrame of bytes, typically from
        the method get_batch(), get_all_batches(), or get_all_vectors()
        of a .load.VectorLoader instance.
    reference: bool = False
        Whether to yield the reference sequence as the first item, prior
        to any relation vectors. The reference sequence will be the same
        length as every relation vector.
    """
    if reference:
        # Prepend the reference sequence to the lines of vectors.
        yield f"Reference\t{parse_seq(vectors.columns).decode()}"
    for index, row in zip(vectors.index, vectors.values, strict=True):
        yield f"{index}\t{humanize_relvec(row).decode()}"


def as_block(vectors: pd.DataFrame, reference: bool = False):
    """ Display all relation vectors as a block of human-readable text,
    with each vector on a new line and all the positions in the vector
    aligned vertically. Parameters are the same as those of as_iter. """
    return "\n".join(as_iter(vectors, reference=reference))

########################################################################
#                                                                      #
# ©2023, the Rouskin Lab.                                              #
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
