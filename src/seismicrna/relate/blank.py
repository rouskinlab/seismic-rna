import numpy as np
import pandas as pd

from ..core.rel import NOCOV, NP_TYPE
from ..core.sect import seq_pos_to_index
from ..core.seq import DNA


def blank_relvec(bases: int | DNA,
                 reads: int | list | np.ndarray | pd.Index | None = None):
    """
    Return blank relation vector(s) of a given length.

    Parameters
    ----------
    bases: int | DNA
        The reference sequence (if DNA) or just its length (if int).
        If the sequence itself is given, then return a Series/DataFrame
        whose main index (index for Series, columns for DataFrame) is a
        MultiIndex of 1-indexed positions and bases in the sequence.
        If only the sequence length is given, then return a NumPy array.
    reads: int | list | np.ndarray | pd.Index | None = None
        Optional names of the relation vectors. If given, then return a
        DataFrame whose indexes are the read names if bases is DNA,
        otherwise a 2D NumPy array with one row for each name in reads.
        If None, then return a Series if bases is DNA, otherwise a 1D
        NumPy array.

    Returns
    -------
    np.ndarray | pd.Series | pd.DataFrame
        The blank relation vector(s).
    """
    # Determine whether to return a Pandas or NumPy object.
    if isinstance(bases, DNA):
        # Make a Series or DataFrame with the sequence as its main axis.
        sequence = seq_pos_to_index(bases, np.arange(1, len(bases) + 1), 1)
        if reads is None:
            # Return a Series representing just one relation vector.
            return pd.Series(NOCOV, index=sequence, dtype=NP_TYPE)
        # Ensure that names is a sequence of read names as str objects.
        if isinstance(reads, int):
            names = [f"Read_{i}" for i in range(1, reads + 1)]
        else:
            names = list(map(str, reads))
        # Return a DataFrame with one row per relation vector.
        return pd.DataFrame(NOCOV, index=names, columns=sequence, dtype=NP_TYPE)
    # Determine the size of the NumPy array.
    size = (bases if reads is None
            else ((reads, bases) if isinstance(reads, int)
                  else (len(reads), bases)))
    return np.full(size, fill_value=NOCOV, dtype=NP_TYPE)

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
