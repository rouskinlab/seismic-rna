from logging import getLogger
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .batch import RelateBatch
from .report import RelateReport
from .seqpos import format_seq_pos, parse_pos
from ..core.data import BatchLoader
from ..core.rel import REL_TYPE
from ..core.sect import seq_pos_to_index
from ..core.seq import DNA

logger = getLogger(__name__)

POS = "by_position"
VEC = "by_vector"
READ = "Read Name"


class RelateLoader(BatchLoader):
    """ Load batches of relation vectors. """

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @classmethod
    def get_batch_type(cls):
        return RelateBatch

    @classmethod
    def get_btype_name(cls):
        return

    def get_refseq(self):
        return DNA.load(self._report.refseq_file_path(self.top))

    def load_data_personal(self, batch_file: Path, *,
                           positions: Sequence[int] | None = None):
        """
        Return the relation vectors from one batch. Optionally, return
        a subset of the positions (columns) in the relation vectors.

        Parameters
        ----------
        batch_file: Path
            File of the batch
        positions: Sequence[int] | None = None
            Positions of the sequence to load (1-indexed).

        Return
        ------
        DataFrame
            Relation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number.
        """
        # Determine which columns to read from the file.
        if positions is None:
            # Load all columns.
            columns = None
        else:
            # Load the columns corresponding to the given positions.
            columns = format_seq_pos(self.seq, positions, self.end5)
        # Read the batch file using the selected positions.
        vectors = pd.read_parquet(batch_file, columns=columns)
        # Convert the columns to a MultiIndex of positions and bases.
        vectors.columns = seq_pos_to_index(self.seq, parse_pos(vectors.columns),
                                           self.end5)
        # Name the index.
        vectors.index.rename(READ, inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work.
        return vectors.astype(REL_TYPE, copy=False)

    def iter_batches_personal(self, *, positions: Sequence[int] | None = None):
        yield from super().iter_batches_personal(positions=positions)

    def iter_batches_processed(self, *, positions: Sequence[int] | None = None):
        yield from super().iter_batches_processed(positions=positions)


def open_reports(report_files: Iterable[Path]):
    """ Load an arbitrary number of vector reports. """
    reports = dict()
    for report_file in report_files:
        try:
            # Load the report and collect basic information.
            report = RelateLoader.open(report_file)
            key = report.sample, report.ref
            if key in reports:
                logger.warning(f"Got multiple reports for {key}")
            else:
                reports[key] = report
        except Exception as error:
            logger.error(f"Failed to open {report_file}: {error}")
    return list(reports.values())

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
