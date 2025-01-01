from abc import ABC, abstractmethod
from functools import cached_property
from typing import Callable, Generator

from .base import (GraphBase,
                   GraphRunner,
                   GraphWriter,
                   get_action_name,
                   make_path_subject,
                   make_title_action_sample)
from .onedata import OneDataGraph
from ..core.data import MutsDataset


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


class OneDatasetGraph(DatasetGraph, OneDataGraph, ABC):

    @cached_property
    def action(self):
        return get_action_name(self.dataset)

    @cached_property
    def path_subject(self):
        return make_path_subject(self.action, self.k, self.clust)

    @cached_property
    def title_action_sample(self):
        return make_title_action_sample(self.action, self.sample)


class DatasetGraphWriter(GraphWriter, ABC):

    def __init__(self, *datasets: MutsDataset):
        self.datasets = list(datasets)

    @abstractmethod
    def iter_graphs(self,
                    *args,
                    **kwargs) -> Generator[DatasetGraph, None, None]:
        """ Yield every graph. """


class DatasetGraphRunner(GraphRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[DatasetGraphWriter]:
        """ Type of GraphWriter. """

    @classmethod
    @abstractmethod
    def get_table_loader(cls) -> Callable[[tuple[str, ...]], Generator]:
        """ Function to find and filter table files. """

    @classmethod
    def list_table_files(cls, input_path: tuple[str, ...]):
        """ Find, filter, and list all table files from input files. """
        finder = cls.get_table_loader()
        return list(finder(input_path))
