from abc import ABC, abstractmethod
from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import Any, Callable, Iterable

from . import path
from .batch import MutsBatch, RegionMutsBatch, ReadBatch, list_batch_nums
from .error import InconsistentValueError, NoDataError
from .header import NO_K
from .io import MutsBatchIO, ReadBatchIO
from .logs import logger
from .rel import RelPattern
from .report import (BranchesF,
                     SampleF,
                     RefF,
                     RegF,
                     TimeBeganF,
                     TimeEndedF,
                     NumBatchesF,
                     ChecksumsF,
                     PooledSamplesF,
                     JoinedRegionsF,
                     Report,
                     BatchedReport,
                     MissingFieldWithNoDefaultError)
from .seq import DNA, Region, unite
from .validate import (require_isinstance,
                       require_issubclass,
                       require_atleast,
                       require_atmost)


class BadTimeStampError(RuntimeError):
    """ A dataset has a timestamp that is earlier than a dataset that
    should have been written before it. """


class MissingBatchError(RuntimeError):
    """ A dataset does not have a batch of a given type and number. """


class MissingBatchTypeError(MissingBatchError):
    """ A dataset does not have a batch of a given type. """


class FailedToLoadDatasetError(RuntimeError):
    """ A batch failed to load. """


class Dataset(ABC):
    """ Dataset comprising batches of data. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[Report]:
        """ Type of report. """

    def __init__(self, report_file: str | Path, verify_times: bool = True):
        self.report_file = Path(report_file)
        self.verify_times = verify_times
        # Load the report here so that if the report does not have the
        # expected fields, an error will be raised inside __init__.
        # Doing so is essential for LoadFunction to determine which type
        # of dataset to load, since it bases that decision on whether
        # __init__ succeeds or raises an error.
        self.report = self.get_report_type().load(self.report_file)

    @cached_property
    def top(self) -> Path:
        """ Top-level directory of the dataset. """
        top, fields = self.get_report_type().parse_path(self.report_file)
        return top

    @property
    def branches(self):
        """ Branches of the workflow. """
        return self.report.get_field(BranchesF)

    @cached_property
    def sample(self) -> str:
        """ Name of the sample. """
        return self.report.get_field(SampleF)

    @cached_property
    def ref(self) -> str:
        """ Name of the reference. """
        return self.report.get_field(RefF)

    @property
    @abstractmethod
    def pattern(self) -> RelPattern | None:
        """ Pattern of mutations to count. """

    @cached_property
    @abstractmethod
    def ks(self) -> list[int]:
        """ Numbers of clusters. """

    @cached_property
    @abstractmethod
    def best_k(self) -> int:
        """ Best number of clusters. """

    @cached_property
    def is_clustered(self):
        """ Whether the dataset is clustered. """
        return any(k != NO_K for k in self.ks)

    @cached_property
    def time_began(self) -> datetime:
        """ Time at which the data were written. """
        return self.report.get_field(TimeBeganF)

    @cached_property
    def time_ended(self) -> datetime:
        """ Time at which the data were written. """
        return self.report.get_field(TimeEndedF)

    @property
    def dir(self) -> Path:
        """ Directory containing the dataset. """
        return self.report_file.parent

    @property
    @abstractmethod
    def data_dirs(self) -> list[Path]:
        """ All directories containing data for the dataset. """

    @property
    @abstractmethod
    def num_batches(self) -> int:
        """ Number of batches. """

    @property
    def batch_nums(self):
        """ Numbers of the batches. """
        return list_batch_nums(self.num_batches)

    @abstractmethod
    def get_batch(self, batch_num: int) -> ReadBatch:
        """ Get a specific batch of data. """

    def iter_batches(self):
        """ Yield each batch. """
        for batch_num in self.batch_nums:
            yield self.get_batch(batch_num)

    @cached_property
    def num_reads(self):
        """ Number of reads in the dataset. """
        return sum(batch.num_reads for batch in self.iter_batches())

    def link_data_dirs_to_tmp(self, tmp_dir: Path):
        """ Make links to a dataset in a temporary directory. """
        for data_dir in self.data_dirs:
            tmp_data_dir = path.transpath(tmp_dir,
                                          self.top,
                                          data_dir,
                                          strict=True)
            path.mkdir_if_needed(tmp_data_dir.parent)
            path.symlink_if_needed(tmp_data_dir, data_dir)

    def __str__(self):
        return f"{type(self).__name__}({self.report_file})"

    def __repr__(self):
        return str(self)


class UnbiasDataset(Dataset, ABC):
    """ Dataset with attributes for correcting observer bias. """

    def __init__(self,
                 *args,
                 masked_read_nums: dict[[int, list]] | None = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        if masked_read_nums is not None:
            require_isinstance("masked_read_nums", masked_read_nums, dict)
        self.masked_read_nums = masked_read_nums

    @property
    @abstractmethod
    def min_mut_gap(self) -> int:
        """ Minimum gap between two mutations. """

    @property
    @abstractmethod
    def quick_unbias(self) -> bool:
        """ Use the quick heuristic for unbiasing. """

    @property
    @abstractmethod
    def quick_unbias_thresh(self) -> float:
        """ Consider mutation rates less than or equal to this threshold
        to be 0 when using the quick heuristic for unbiasing. """


class RegionDataset(Dataset, ABC):
    """ Dataset with a known reference sequence and region. """

    @property
    @abstractmethod
    def refseq(self) -> DNA:
        """ Sequence of the reference. """

    @property
    def reflen(self):
        """ Length of the reference sequence. """
        return len(self.refseq)

    @cached_property
    @abstractmethod
    def region(self) -> Region:
        """ Region of the dataset. """


class MutsDataset(RegionDataset, ABC):
    """ Dataset with a known region and explicit mutational data. """

    @abstractmethod
    def get_batch(self, batch_num: int) -> RegionMutsBatch:
        # Redeclare this method (already declared by Dataset) for the
        # sole purpose of enabling type linters to determine that
        # MutsDataset.get_batch() returns a RegionMutsBatch object
        # instead of plain ReadBatch objects.
        pass

    def get_batch_count_all(self, batch_num: int, **kwargs):
        """ Calculate the counts for a specific batch of data. """
        return self.get_batch(batch_num).count_all(**kwargs)

    def iter_batches(self):
        # Reimplement this method (already implemented by Dataset) for
        # the sole purpose of enabling type linters to determine that
        # MutsDataset.iter_batches() yields RegionMutsBatch objects
        # instead of plain ReadBatch objects.
        for batch_num in self.batch_nums:
            yield self.get_batch(batch_num)


class LoadFunction(object):
    """ Function to load a dataset. """

    def __init__(self,
                 dataset_type: type[Dataset], /,
                 *more_types: type[Dataset]):
        require_issubclass("dataset_type", dataset_type, Dataset)
        for i in range(len(more_types)):
            require_issubclass(f"more_types[{i}]", more_types[i], Dataset)
        self.dataset_types = (dataset_type,) + more_types

    def _get_dataset_type_consensus(self,
                                    method: Callable[[type[Dataset]], Any]):
        """ Get the consensus value among all types of dataset. """
        value0 = method(self.dataset_types[0])
        for dataset_type in self.dataset_types[1:]:
            assert method(dataset_type) == value0
        return value0

    @cached_property
    def report_path_seg_types(self):
        """ Segment types of the report file path. """
        return self._get_dataset_type_consensus(
            lambda dt: dt.get_report_type().get_path_seg_types()
        )

    @cached_property
    def report_path_auto_fields(self):
        """ Automatic field values of the report file path. """
        return self._get_dataset_type_consensus(
            lambda dt: dt.get_report_type().get_auto_path_fields()
        )

    def build_report_path(self, path_fields: dict[str, Any]):
        """ Build the path of a report file. """
        return path.build(self.report_path_seg_types,
                          {**self.report_path_auto_fields, **path_fields})

    def __call__(self, report_file: str | Path, **kwargs):
        """ Load a dataset from the report file. """
        # Try to load the report file using each type of dataset.
        errors = dict()
        for dataset_type in self.dataset_types:
            try:
                # Return the first dataset type that works.
                return dataset_type(report_file, **kwargs)
            except FileNotFoundError:
                # Re-raise FileNotFoundError because if the report file
                # does not exist, then no dataset type can load it.
                raise
            except BadTimeStampError:
                # Re-raise BadTimeStampError because if the report file
                # has a timestamp that is earlier than that of one of
                # its constituents, then no dataset type can load it.
                raise
            except MissingFieldWithNoDefaultError as error:
                # If a dataset type fails because of missing a default
                # value, then it was probably the wrong datasettype, so
                # record the error silently.
                errors[dataset_type.__name__] = error
        # If all dataset types failed, then raise an error.
        errmsg = "\n".join(f"{type_name}: {error}"
                           for type_name, error in errors.items())
        raise FailedToLoadDatasetError(
            f"{self} failed to load {report_file}:\n{errmsg}"
        )

    def iterate(self, input_path: Iterable[str | Path], *, raise_on_error : bool = False, **kwargs):
        """ Yield a Dataset from each report file in `input_path`. """
        for report_file in path.find_files_chain(input_path,
                                                 self.report_path_seg_types):
            try:
                yield self(report_file, **kwargs)
            except Exception as error:
                if raise_on_error:
                    raise error
                logger.error(error)

    def __str__(self):
        names = ", ".join(dataset_type.__name__
                          for dataset_type in self.dataset_types)
        return f"{type(self).__name__}({names})"


class LoadedDataset(Dataset, ABC):
    """ Dataset created by loading directly from a Report. """

    @classmethod
    @abstractmethod
    def get_batch_type(cls) -> type[ReadBatchIO | MutsBatchIO]:
        """ Type of batch. """

    @classmethod
    @abstractmethod
    def get_report_type(cls) -> type[BatchedReport]:
        """ Type of report. """

    @classmethod
    def get_btype_name(cls):
        """ Name of the type of batch. """
        return cls.get_batch_type().btype()

    @property
    def num_batches(self):
        return self.report.get_field(NumBatchesF)

    @property
    def data_dirs(self):
        return [self.dir]

    def get_batch_path(self, batch_num: int):
        """ Get the path to a batch of a specific number. """
        require_atleast("batch_num", batch_num, 0, classes=int)
        require_atmost("batch_num",
                       batch_num,
                       self.num_batches - 1,
                       "num_batches",
                       classes=int)
        return path.cast_path(self.report_file,
                              self.report.get_path_seg_types(),
                              self.get_batch_type().get_path_seg_types(),
                              {path.BATCH: batch_num,
                               path.EXT: path.BRICKLE_EXT})

    def get_batch_checksum(self, batch_num: int):
        """ Get the checksum of a specific batch from the report. """
        require_atleast("batch_num", batch_num, 0, classes=int)
        require_atmost("batch_num",
                       batch_num,
                       self.num_batches - 1,
                       "num_batches",
                       classes=int)
        checksums = self.report.get_field(ChecksumsF)
        try:
            return checksums[self.get_btype_name()][batch_num]
        except KeyError:
            # Report does not have checksums for this type of batch.
            raise MissingBatchTypeError(self.get_batch_type())
        except IndexError:
            # Report does not have a checksum for this batch number.
            raise MissingBatchError(batch_num)

    def get_batch(self, batch_num: int) -> ReadBatchIO | MutsBatchIO:
        batch = self.get_batch_type().load(
            self.get_batch_path(batch_num),
            checksum=self.get_batch_checksum(batch_num)
        )
        assert batch.batch == batch_num
        return batch


class MergedDataset(Dataset, ABC):
    """ Dataset made by merging one or more constituent datasets. """

    @classmethod
    @abstractmethod
    def get_dataset_load_func(cls) -> LoadFunction:
        """ Function to load one constituent dataset. """

    @cached_property
    @abstractmethod
    def datasets(self) -> list[Dataset]:
        """ Constituent datasets that were merged. """

    @cached_property
    def data_dirs(self):
        return [self.dir] + [dataset_dir
                             for dataset in self.datasets
                             for dataset_dir in dataset.data_dirs]

    def _list_dataset_attr(self, name: str, *subnames: str):
        """ Get a list of an attribute for each dataset. """
        values = list()
        for dataset in self.datasets:
            value = getattr(dataset, name)
            for subname in subnames:
                value = getattr(value, subname)
            values.append(value)
        return values

    def _get_common_attr(self, name: str, *subnames):
        """ Get a common attribute among datasets. """
        values = list()
        for value in self._list_dataset_attr(name, *subnames):
            if value not in values:
                values.append(value)
        if len(values) != 1:
            raise InconsistentValueError(
                f"Attribute {repr(name)} has multiple values: {values}"
            )
        return values[0]

    @cached_property
    def pattern(self):
        return self._get_common_attr("pattern")


class MergedRegionDataset(MergedDataset, RegionDataset, ABC):

    @cached_property
    def refseq(self):
        return self._get_common_attr("refseq")


class MergedUnbiasDataset(MergedDataset, UnbiasDataset, ABC):
    """ MergedDataset with attributes for correcting observer bias. """

    @property
    def min_mut_gap(self):
        return self._get_common_attr("min_mut_gap")

    @property
    def quick_unbias(self):
        return self._get_common_attr("quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return self._get_common_attr("quick_unbias_thresh")


class TallDataset(MergedDataset, ABC):
    """ Dataset made by vertically pooling other datasets from one or
    more samples aligned to the same reference sequence. """

    @cached_property
    def datasets(self):
        # Determine the sample name and pooled samples from the report.
        pooled_samples = self.report.get_field(PooledSamplesF)
        # Determine the report file for each sample in the pool.
        load_func = self.get_dataset_load_func()
        sample_report_files = [path.cast_path(
            self.report_file,
            self.get_report_type().get_path_seg_types(),
            load_func.report_path_seg_types,
            {**load_func.report_path_auto_fields, path.SAMPLE: sample}
        ) for sample in pooled_samples]
        if not sample_report_files:
            raise NoDataError(f"{self} has no datasets")
        return list(map(self.get_dataset_load_func(),
                        sample_report_files))

    @cached_property
    def samples(self) -> list[str]:
        """ Names of all samples in the pool. """
        return self._list_dataset_attr("sample")

    @cached_property
    def nums_batches(self) -> list[int]:
        """ Number of batches in each dataset in the pool. """
        return self._list_dataset_attr("num_batches")

    @cached_property
    def num_batches(self):
        return sum(self.nums_batches)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._batch_nums_table = dict()

    def _translate_batch_num(self, batch_num: int):
        """ Translate a batch number into the numbers of the dataset and
        of the batch within the dataset. """
        require_atleast("batch_num", batch_num, 0, classes=int)
        require_atmost("batch_num",
                       batch_num,
                       self.num_batches - 1,
                       "num_batches",
                       classes=int)
        try:
            # Check if the translated batch number has been cached.
            dataset_num, dataset_batch_num = self._batch_nums_table[batch_num]
            return dataset_num, dataset_batch_num
        except KeyError:
            # If not, then calculate it.
            dataset_batch_num = batch_num
            for dataset_num, num_batches in enumerate(self.nums_batches):
                assert dataset_batch_num >= 0
                if dataset_batch_num < num_batches:
                    # Cache the translated batch number.
                    self._batch_nums_table[batch_num] = (dataset_num,
                                                         dataset_batch_num)
                    return dataset_num, dataset_batch_num
                dataset_batch_num -= num_batches
        raise ValueError(f"{self} has no batch with number {batch_num}")

    def get_batch(self, batch_num: int):
        # Determine the dataset and the batch number in that dataset.
        dataset_num, dataset_batch_num = self._translate_batch_num(batch_num)
        # Load the batch.
        batch = self.datasets[dataset_num].get_batch(dataset_batch_num)
        # Renumber the batch from the numbering in its original dataset
        # to the numbering in the pooled dataset.
        batch.batch = batch_num
        return batch


class WideDataset(MergedRegionDataset, ABC):
    """ Dataset made by horizontally joining other datasets from one or
    more regions of the same reference sequence. """

    @cached_property
    def datasets(self):
        # Determine the name of the joined region and the individual
        # regions from the report.
        joined_regs = self.report.get_field(JoinedRegionsF)
        # Determine the report file for each joined region.
        load_func = self.get_dataset_load_func()
        region_report_files = [path.cast_path(
            self.report_file,
            self.get_report_type().get_path_seg_types(),
            load_func.report_path_seg_types,
            {**load_func.report_path_auto_fields, path.REG: reg}
        ) for reg in joined_regs]
        if not region_report_files:
            raise NoDataError(f"{self} has no datasets")
        return list(map(self.get_dataset_load_func(),
                        region_report_files))

    @cached_property
    def num_batches(self):
        return self._get_common_attr("num_batches")

    @cached_property
    def region_names(self):
        """ Names of all joined regions. """
        return self._list_dataset_attr("region", "name")

    @cached_property
    def region(self):
        return unite(self._list_dataset_attr("region"),
                     name=self.report.get_field(RegF),
                     refseq=self.refseq)

    @abstractmethod
    def _join(self, batches: Iterable[tuple[str, ReadBatch]]) -> ReadBatch:
        """ Join corresponding batches of data. """

    def get_batch(self, batch_num: int):
        # Join the batch with that number from every dataset.
        return self._join((dataset.region.name, dataset.get_batch(batch_num))
                          for dataset in self.datasets)


class WideMutsDataset(WideDataset, MutsDataset, ABC):
    """ WideDataset with mutation data. """


class MultistepDataset(MutsDataset, ABC):
    """ Dataset made by integrating two datasets from different steps of
    the workflow. """

    @classmethod
    @abstractmethod
    def get_dataset1_load_func(cls) -> LoadFunction:
        """ Function to load Dataset 1. """

    @classmethod
    @abstractmethod
    def get_dataset2_type(cls) -> type[RegionDataset]:
        """ Type of Dataset 2. """

    @classmethod
    def get_dataset2_load_func(cls):
        """ Function to load Dataset 2. """
        return LoadFunction(cls.get_dataset2_type())

    @classmethod
    def get_report_type(cls):
        # The report is for dataset 2.
        return cls.get_dataset2_type().get_report_type()

    @classmethod
    def get_dataset1_report_file(cls,
                                 dataset2_report_file: Path,
                                 verify_times: bool):
        """ Given the report file for Dataset 2, determine the report
        file for Dataset 1. """
        load_func = cls.get_dataset1_load_func()
        dataset2 = cls.load_dataset2(dataset2_report_file, verify_times)
        return path.cast_path(
            dataset2_report_file,
            cls.get_report_type().get_path_seg_types(),
            load_func.report_path_seg_types,
            # Replace the default fields and branches of report path 2
            # with those of report path 1.
            {**load_func.report_path_auto_fields,
             path.BRANCHES: path.get_ancestors(dataset2.branches)}
        )

    @classmethod
    def load_dataset1(cls, dataset2_report_file: Path, verify_times: bool):
        """ Load Dataset 1. """
        load_func = cls.get_dataset1_load_func()
        return load_func(cls.get_dataset1_report_file(dataset2_report_file,
                                                      verify_times),
                         verify_times=verify_times)

    @classmethod
    def load_dataset2(cls, dataset2_report_file: Path, verify_times: bool):
        """ Load Dataset 2. """
        load_func = cls.get_dataset2_load_func()
        return load_func(dataset2_report_file,
                         verify_times=verify_times)

    def __init__(self, dataset2_report_file: Path, **kwargs):
        super().__init__(dataset2_report_file, **kwargs)
        dataset1 = self.load_dataset1(dataset2_report_file, self.verify_times)
        dataset2 = self.load_dataset2(dataset2_report_file, self.verify_times)
        time1 = dataset1.time_ended
        time2 = dataset2.time_began
        if self.verify_times and time1 > time2:
            raise BadTimeStampError(
                f"{dataset1.report_file} should have existed before "
                f"{dataset2.report_file} but was actually created afterwards. "
                f"If you are sure this inconsistency is not a problem, "
                f"then you can suppress this error using --no-verify-times"
            )
        self.dataset1 = dataset1
        self.dataset2 = dataset2

    @property
    def refseq(self):
        return self.dataset1.refseq

    @cached_property
    def data_dirs(self):
        return self.dataset1.data_dirs + self.dataset2.data_dirs

    @abstractmethod
    def _integrate(self, batch1: MutsBatch, batch2: ReadBatch) -> MutsBatch:
        """ Integrate corresponding batches of data. """

    @property
    def num_batches(self):
        return self.report.get_field(NumBatchesF)

    def get_batch(self, batch_num: int):
        return self._integrate(self.dataset1.get_batch(batch_num),
                               self.dataset2.get_batch(batch_num))
