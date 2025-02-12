import gzip
from abc import ABC, abstractmethod
from functools import cache
from pathlib import Path
from typing import Iterable, Self

import numpy as np
import pandas as pd

from . import path
from .arg import opt_max_fmut_pos, opt_min_ninfo_pos
from .io import RefFileIO
from .logs import logger
from .seq import FIELD_REF, POS_NAME
from .table import (INFOR_REL,
                    MUTAT_REL,
                    RelTypeTable,
                    PositionTableLoader,
                    RelTypeTableLoader)


class List(RefFileIO, ABC):
    """ List base class. """

    NAMES = [FIELD_REF, POS_NAME]

    @classmethod
    @abstractmethod
    def get_table_type(cls) -> type[RelTypeTableLoader]:
        """ Type of table that this type of list can process. """

    @classmethod
    def get_by_read(cls) -> bool:
        """ Whether the list is of reads. """
        return cls.get_table_type().get_by_read()

    @classmethod
    def get_is_gzip(cls):
        """ Whether the file is compressed with gzip. """
        return cls.get_by_read()

    @classmethod
    def get_ext(cls):
        return path.CSVZIP_EXT if cls.get_is_gzip() else path.CSV_EXT

    @classmethod
    def get_open_func(cls):
        return gzip.open if cls.get_is_gzip() else open

    @classmethod
    @cache
    def get_auto_path_fields(cls):
        """ Default values of the path fields. """
        return {path.STEP: cls.get_step(),
                path.LIST: cls.get_step(),
                **super().get_auto_path_fields()}

    @classmethod
    def get_path_from_table(cls, table: RelTypeTable, branch: str):
        """ Get the path of a list given a table. """
        return path.cast_path(table.path,
                              table.get_path_seg_types(),
                              cls.get_path_seg_types(),
                              {**cls.get_auto_path_fields(),
                               path.BRANCHES: path.add_branch(path.LIST_STEP,
                                                              branch,
                                                              table.branches)})

    @classmethod
    def load(cls, file: Path):
        top, path_field_values = cls.parse_path(file, exclude_auto=True)
        return cls(data=pd.read_csv(file, index_col=cls.NAMES).index.to_frame(),
                   **path_field_values)

    @classmethod
    @abstractmethod
    def from_table(cls, table: RelTypeTable, branch: str, **kwargs) -> Self:
        """ Create a list from a table. """

    def __init__(self, *,
                 sample: str,
                 branches: Iterable[str],
                 ref: str,
                 data: pd.DataFrame,
                 **kwargs):
        super().__init__(**kwargs)
        self.sample = sample
        if isinstance(branches, dict):
            branches = path.flatten_branches(branches)
        else:
            branches = list(branches)
            path.validate_branches_flat(branches)
        self.branches = branches
        self.ref = ref
        if not isinstance(data, pd.DataFrame):
            raise TypeError(
                f"data must be {pd.DataFrame}, but got {type(data)}"
            )
        if data.columns.to_list() != self.NAMES:
            raise ValueError(f"data must have columns {self.NAMES}, "
                             f"but got {data.columns.to_list()}")
        self.data = data

    def save(self, top: Path, force: bool = False):
        file = self.get_path(top)
        if not force and file.exists():
            raise FileExistsError(file)
        self.data.to_csv(file, header=True, index=False)
        return file


class PositionList(List, ABC):
    """ List of positions. """

    MASK_NINFO = "list_min_ninfo_pos"
    MASK_FMUT = "list_max_fmut_pos"

    @classmethod
    def get_data_type(cls):
        return int

    @classmethod
    def get_file_seg_type(cls):
        return path.PositionListSeg

    @classmethod
    def from_table(cls,
                   table: PositionTableLoader,
                   branch: str, *,
                   min_ninfo_pos: int = opt_min_ninfo_pos.default,
                   max_fmut_pos: float = opt_max_fmut_pos.default):
        logger.routine(f"Began making {cls} from {table}")
        if not isinstance(table, cls.get_table_type()):
            raise TypeError(
                f"table must be {cls.get_table_type()}, but got {type(table)}"
            )
        region = table.region.copy(masks=False)
        logger.detail(f"{table} has {region.length} positions of which "
                      f"{region.unmasked_int.size} are in use "
                      f"and {region.masked_int.size} are masked")
        region.add_mask(cls.MASK_NINFO,
                        region.range_int[table.fetch_count(rel=INFOR_REL,
                                                           squeeze=True)
                                         < min_ninfo_pos])
        logger.detail(f"{table} has {region.get_mask(cls.MASK_NINFO).size} "
                      f"positions with < {min_ninfo_pos} informative bases")
        region.add_mask(cls.MASK_FMUT,
                        region.range_int[table.fetch_ratio(rel=MUTAT_REL,
                                                           squeeze=True)
                                         > max_fmut_pos])
        logger.detail(f"{table} has {region.get_mask(cls.MASK_FMUT).size} "
                      f"positions with mutation rates > {max_fmut_pos}")
        positions = region.masked_int
        data = pd.DataFrame.from_dict({FIELD_REF: np.repeat(table.ref,
                                                            positions.size),
                                       POS_NAME: positions})
        new_list = cls(data=data,
                       sample=table.sample,
                       branches=path.add_branch(path.LIST_STEP,
                                                branch,
                                                table.branches),
                       ref=table.ref)
        logger.detail(f"{table} was used to produce a list of {positions.size} "
                      "positions to mask")
        logger.routine(f"Ended making {cls} from {table}")
        return new_list


class ReadList(List, ABC):
    """ List of reads. """

    @classmethod
    def get_data_type(cls):
        return str

    @classmethod
    def get_file_seg_type(cls):
        return path.ReadListSeg
