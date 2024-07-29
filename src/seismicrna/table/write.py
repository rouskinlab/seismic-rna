from abc import ABC
from functools import cached_property
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
                   ClustFreqTable)
from .calc import (Tabulator,
                   AvgTabulator,
                   RelateTabulator,
                   MaskTabulator,
                   ClustTabulator,
                   tabulate_loader)
from ..core.data import LoadFunction
from ..core.write import need_write
from ..cluster.data import ClusterMutsDataset, JoinClusterMutsDataset
from ..mask.data import MaskMutsDataset, JoinMaskMutsDataset
from ..pool.data import PoolDataset
from ..relate.data import RelateDataset

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
            self.data.round(decimals=PRECISION).to_csv(self.path)
        return self.path


# Write by Index (position/read/cluster) ###############################

class PosTableWriter(TableWriter, PosTable, ABC):

    @cached_property
    def data(self):
        return self.tabulator.table_per_pos


class ReadTableWriter(TableWriter, ReadTable, ABC):

    @cached_property
    def data(self):
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


class ClustFreqTableWriter(TableWriter, ClustFreqTable):

    @cached_property
    def data(self):
        return self.tabulator.table_per_clust


# Helper Functions #####################################################

def get_tabulator_writer_types(tabulator: Tabulator):
    if isinstance(tabulator, RelateTabulator):
        return RelPosTableWriter, RelReadTableWriter
    if isinstance(tabulator, MaskTabulator):
        return MaskPosTableWriter, MaskReadTableWriter
    if isinstance(tabulator, ClustTabulator):
        return ClustPosTableWriter, ClustFreqTableWriter
    raise TypeError(f"Invalid tabulator type: {type(tabulator).__name__}")


def get_tabulator_writers(tabulator: AvgTabulator | ClustTabulator, *,
                          table_pos: bool = True,
                          table_read: bool = True,
                          table_clust: bool = True):
    types_write = {PosTableWriter: table_pos,
                   ReadTableWriter: table_read,
                   ClustFreqTableWriter: table_clust}
    for writer_type in get_tabulator_writer_types(tabulator):
        for table_type, type_write in types_write.items():
            if issubclass(writer_type, table_type):
                if type_write:
                    yield writer_type(tabulator)
                else:
                    logger.debug(f"Skipped {writer_type} for {tabulator}")
                break
        else:
            raise TypeError(f"Invalid writer type: {writer_type.__name__}")


load_any_dataset = LoadFunction(RelateDataset,
                                PoolDataset,
                                MaskMutsDataset,
                                JoinMaskMutsDataset,
                                ClusterMutsDataset,
                                JoinClusterMutsDataset)


def write(report_file: Path, *, force: bool, **kwargs):
    """ Helper function to write a table from a report file. """
    # Load the dataset.
    dataset = load_any_dataset(report_file)
    # Create the tabulator for the dataset.
    tabulator = tabulate_loader(dataset)
    # For each table associated with this tabulator, create the table,
    # write it, and return the path to the table output file.
    return [table.write(force)
            for table in get_tabulator_writers(tabulator, **kwargs)]

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
