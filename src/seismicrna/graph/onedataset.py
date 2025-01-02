from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain

from .base import (cgroup_table,
                   get_action_name,
                   make_path_subject,
                   make_title_action_sample)
from .dataset import DatasetGraph, DatasetGraphWriter, DatasetGraphRunner
from .onesource import OneSourceGraph
from ..core.data import MutsDataset
from ..core.task import dispatch


class OneDatasetGraph(DatasetGraph, OneSourceGraph, ABC):

    @cached_property
    def action(self):
        return get_action_name(self.dataset)

    @cached_property
    def path_subject(self):
        return make_path_subject(self.action, self.k, self.clust)

    @cached_property
    def title_action_sample(self):
        return make_title_action_sample(self.action, self.sample)


class OneDatasetWriter(DatasetGraphWriter, ABC):

    def __init__(self, dataset: MutsDataset):
        super().__init__(dataset)

    @cached_property
    def dataset(self):
        """ The dataset providing the data for the graph(s). """
        return self.datasets[0]

    @abstractmethod
    def get_graph(self, *args, **kwargs) -> OneDatasetGraph:
        """ Return a graph instance. """

    def iter_graphs(self,
                    rels: tuple[str, ...],
                    cgroup: str,
                    **kwargs):
        for cparams in cgroup_table(self.dataset, cgroup):
            for rel in rels:
                yield self.get_graph(rel, **kwargs | cparams)


class OneDatasetRunner(DatasetGraphRunner, ABC):

    @classmethod
    def run(cls,
            input_path: tuple[str, ...], *,
            max_procs: int,
            **kwargs):
        # Generate a table writer for each table.
        writer_type = cls.get_writer_type()
        writers = [writer_type(table_file)
                   for table_file in cls.list_input_files(input_path)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    pass_n_procs=False,
                                    kwargs=kwargs)))
