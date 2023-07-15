from abc import ABC, abstractmethod
from logging import getLogger
from pathlib import Path

import pandas as pd

from .base import (Table, PosTable, ReadTable,
                   RelPosTable, RelReadTable, MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustFreqTable)
from .tabulate import (Tabulator, RelTabulator, MaskTabulator, ClusterTabulator,
                       tabulate_loader)
from ..cluster.load import ClustLoader
from ..core import path
from ..mask.load import MaskLoader
from ..relate.load import RelateLoader

logger = getLogger(__name__)

PRECISION = 1


# Table Writer Base Classes ############################################

class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: Tabulator | ClusterTabulator):
        self._tab = tabulator

    @property
    def out_dir(self):
        return self._tab.out_dir

    @property
    def sample(self):
        return self._tab.sample

    @property
    def ref(self):
        return self._tab.ref

    @property
    def sect(self):
        return self._tab.sect

    @abstractmethod
    def load_data(self):
        """ Load the table's data from a DataLoader. """
        return pd.DataFrame()

    @property
    def data(self):
        return self.load_data()

    def write(self, rerun: bool):
        """ Write the table's rounded data to the table's CSV file. """
        # Check if the File exists.
        if rerun or not self.path.is_file():
            self.data.round(decimals=PRECISION).to_csv(self.path)
        else:
            logger.warning(f"File exists: {self.path}")
        return self.path


# Write by Index (position/read/cluster) ###############################

class PosTableWriter(TableWriter, PosTable, ABC):

    def load_data(self):
        # Load the data for each position, including excluded positions.
        all_positions = self._tab.section.range
        return self._tab.tabulate_by_pos().reindex(index=all_positions)


class ReadTableWriter(TableWriter, ReadTable, ABC):

    def load_data(self):
        return self._tab.tabulate_by_read()


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

    @classmethod
    def clusters_on_columns(cls):
        return True


class ClustReadTableWriter(ReadTableWriter, ClustReadTable):

    @classmethod
    def clusters_on_columns(cls):
        return True


class ClustFreqTableWriter(TableWriter, ClustFreqTable):

    @classmethod
    def clusters_on_columns(cls):
        return False

    def load_data(self):
        return self._tab.tabulate_by_clust()


# Helper Functions #####################################################

def infer_report_loader_type(report_file: Path):
    """ Given a report file path, infer the type of Loader it needs. """
    if path.RelateRepSeg.ptrn.match(report_file.name):
        return RelateLoader
    if path.MaskRepSeg.ptrn.match(report_file.name):
        return MaskLoader
    if path.ClustRepSeg.ptrn.match(report_file.name):
        return ClustLoader
    raise ValueError(f"Failed to infer loader for {report_file}")


def get_tabulator_writer_types(tabulator: Tabulator):
    if isinstance(tabulator, RelTabulator):
        return RelPosTableWriter, RelReadTableWriter
    if isinstance(tabulator, MaskTabulator):
        return MaskPosTableWriter, MaskReadTableWriter
    if isinstance(tabulator, ClusterTabulator):
        return ClustPosTableWriter, ClustReadTableWriter, ClustFreqTableWriter
    raise TypeError(f"Invalid tabulator type: {type(tabulator).__name__}")


def get_tabulator_writers(tabulator: Tabulator):
    for writer_type in get_tabulator_writer_types(tabulator):
        yield writer_type(tabulator)


def write(report_file: Path, rels: str, rerun: bool):
    """ Helper function to write a table from a report file. """
    # Determine the needed type of report loader.
    report_loader_type = infer_report_loader_type(report_file)
    # Load the report.
    report_loader = report_loader_type.open(report_file)
    # Create the tabulator for the report's data.
    tabulator = tabulate_loader(report_loader, rels)
    # For each table associated with this tabulator, create the table,
    # write it, and return the path to the table output file.
    return [table.write(rerun) for table in get_tabulator_writers(tabulator)]
