from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Generator, Iterable

from .base import GraphBase, GraphRunner, GraphWriter
from ..cluster.data import load_cluster_dataset
from ..core import path
from ..core.data import MutsDataset
from ..mask.data import load_mask_dataset
from ..relate.data import load_relate_dataset


class DatasetGraph(GraphBase, ABC):
    """ Graph based on one or more tables. """

    def __init__(self, *, dataset: MutsDataset, **kwargs):
        super().__init__(**kwargs)
        self.dataset = dataset

    @property
    def top(self):
        return self.dataset.top

    @property
    def sample(self):
        return self.dataset.sample

    @property
    def ref(self):
        return self.dataset.ref

    @property
    def reg(self):
        return self.dataset.region.name

    @property
    def seq(self):
        return self.dataset.region.seq

    @cached_property
    def _title_main(self):
        return [f"{self.what()} between "
                f"{self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over reference {repr(self.ref)} "
                f"region {repr(self.reg)}"]

    @property
    def details(self):
        return []

    @cached_property
    def predicate(self):
        return self.codestring


class DatasetGraphWriter(GraphWriter, ABC):

    def __init__(self, *datasets: MutsDataset):
        self.datasets = list(datasets)

    @abstractmethod
    def iter_graphs(self,
                    *args,
                    **kwargs) -> Generator[DatasetGraph, None, None]:
        pass


def load_datasets(input_path: Iterable[str | Path]):
    for load_func in (load_relate_dataset,
                      load_mask_dataset,
                      load_cluster_dataset):
        yield from map(load_func,
                       path.find_files_chain(input_path,
                                             load_func.report_path_seg_types))


class DatasetGraphRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[DatasetGraphWriter]:
        pass

    @classmethod
    def get_input_loader(cls):
        return load_datasets

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
