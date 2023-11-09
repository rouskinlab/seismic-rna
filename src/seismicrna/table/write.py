from abc import ABC
from logging import getLogger
from pathlib import Path

from .base import (Table,
                   PosTable,
                   ReadTable,
                   RelPosTable,
                   RelReadTable,
                   MaskPosTable,
                   MaskReadTable,
                   ClustPosTable,
                   ClustReadTable,
                   ClustFreqTable)
from .calc import (Tabulator,
                   AvgTabulator,
                   RelateTabulator,
                   MaskTabulator,
                   ClustTabulator,
                   tabulate_loader)
from ..clust.data import ClustMerger
from ..core import path
from ..core.write import need_write
from ..mask.data import MaskMerger
from ..relate.data import RelateLoader

logger = getLogger(__name__)

PRECISION = 1


# Table Writer Base Classes ############################################

class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: AvgTabulator | ClustTabulator):
        self.tabulator = tabulator

    @property
    def top(self):
        return self.tabulator.top

    @property
    def sample(self):
        return self.tabulator.sample

    @property
    def ref(self):
        return self.tabulator.ref

    @property
    def sect(self):
        return self.tabulator.section.name

    @property
    def columns(self):
        return self.header.index

    def write(self, force: bool):
        """ Write the table's rounded data to the table's CSV file. """
        if need_write(self.path, force):
            # Write self._data instead of self.data because the former
            # includes any positions that were masked out, while these
            # positions are not present in the latter.
            data = self._data.T if self.transposed() else self._data
            data.round(decimals=PRECISION).to_csv(self.path)
        return self.path


# Write by Index (position/read/cluster) ###############################

class PosTableWriter(TableWriter, PosTable, ABC):

    @property
    def _data(self):
        return self.tabulator.table_per_pos


class ReadTableWriter(TableWriter, ReadTable, ABC):

    @property
    def _data(self):
        return self.tabulator.table_per_read


# Instantiable Table Writers ###########################################

class RelPosTableWriter(PosTableWriter, RelPosTable):
    pass


class RelReadTableWriter(ReadTableWriter, RelReadTable):
    pass


class MaskPosTableWriter(PosTableWriter, MaskPosTable):
    pass


class MaskReadTableWriter(ReadTableWriter, MaskReadTable):
    pass


class ClustPosTableWriter(PosTableWriter, ClustPosTable):
    pass


class ClustReadTableWriter(ReadTableWriter, ClustReadTable):
    pass


class ClustFreqTableWriter(TableWriter, ClustFreqTable):

    @property
    def _data(self):
        return self.tabulator.table_per_clust


# Helper Functions #####################################################

def infer_report_loader_type(report_file: Path):
    """ Given a report file path, infer the type of Loader it needs. """
    if path.RelateRepSeg.ptrn.match(report_file.name):
        return RelateLoader
    if path.MaskRepSeg.ptrn.match(report_file.name):
        return MaskMerger
    if path.ClustRepSeg.ptrn.match(report_file.name):
        return ClustMerger
    raise ValueError(f"Failed to infer loader type for {report_file}")


def get_tabulator_writer_types(tabulator: Tabulator):
    if isinstance(tabulator, RelateTabulator):
        return RelPosTableWriter, RelReadTableWriter
    if isinstance(tabulator, MaskTabulator):
        return MaskPosTableWriter, MaskReadTableWriter
    if isinstance(tabulator, ClustTabulator):
        return ClustPosTableWriter, ClustReadTableWriter, ClustFreqTableWriter
    raise TypeError(f"Invalid tabulator type: {type(tabulator).__name__}")


def get_tabulator_writers(tabulator: AvgTabulator | ClustTabulator):
    for writer_type in get_tabulator_writer_types(tabulator):
        yield writer_type(tabulator)


def write(report_file: Path, force: bool):
    """ Helper function to write a table from a report file. """
    # Determine the needed type of report loader.
    report_loader_type = infer_report_loader_type(report_file)
    # Load the report.
    report_loader = report_loader_type.load(report_file)
    # Create the tabulator for the report's data.
    tabulator = tabulate_loader(report_loader)
    # For each table associated with this tabulator, create the table,
    # write it, and return the path to the table output file.
    return [table.write(force) for table in get_tabulator_writers(tabulator)]

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
