from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from inspect import getmembers
import json
from functools import cache
from logging import getLogger
from math import isclose, isnan, inf, nan
from numbers import Integral
from pathlib import Path
import re
import sys
from typing import Any, Hashable, Callable, Iterable

import numpy as np

from . import path
from .bitcall import SemiBitCaller
from .files import digest_file
from .seq import DNA

logger = getLogger(__name__)


# Field class

class Field(object):
    __slots__ = ["key", "title", "dtype", "iconv", "oconv",
                 "check_val", "check_rep_val"]

    def __init__(self, key: str, title: str, dtype: type, *,
                 iconv: Callable[[Any], Any] | None = None,
                 oconv: Callable[[Any], Any] | None = None,
                 check_val: Callable[[Any], bool] | None = None,
                 check_rep_val: Callable[[Report, Any], bool] | None = None):
        """
        Parameters
        ----------
        key: str
            Key under which the field is stored as an attribute of the
            Report instance.
        title: str
            Title by which the field is identified in the report file.
        dtype: type
            Data type of the field
        iconv: Callable[[Any], Any] | None = None
            Convert an input value to the right type for instantiation.
            If omitted, just cast the value to the proper data type.
        oconv: Callable[[Any], Any] | None = None
            Convert an output value to the right type for exporting.
            If omitted, export the value with no change.
        check_val: Callable[[Any], bool] | None = None
            Validate the value of the field (after validating its type)
            upon assigning the value to an attribute of a Report object.
            Must accept one argument (the value) and return True if the
            value is valid and False otherwise. If None, the value is
            not validated -- only its data type.
        check_rep_val: Callable[[Report, Any], bool] | None = None
            Validate the value of the field (after validating its type)
            upon assigning the value to an attribute of a Report object.
            Must accept two arguments (the Report object and the value)
            and return True if the value is valid and False otherwise.
            If None, the value is not validated -- only its data type.
        """
        self.key = key
        self.title = title
        self.dtype = dtype
        self.iconv = self.dtype if iconv is None else iconv
        self.oconv = oconv
        self.check_val = check_val
        self.check_rep_val = check_rep_val

    def validate(self, report: Report, value: Any):
        # Validate the type.
        if not isinstance(value, self.dtype):
            raise TypeError(f"{self} expected value to be {self.dtype}, "
                            f"but got {type(value)}")
        # Validate the value.
        if self.check_val is not None:
            if not self.check_val(value):
                raise ValueError(f"{self} got invalid value: {value}")
        if self.check_rep_val is not None:
            if not self.check_rep_val(report, value):
                raise ValueError(f"{self} got invalid value: {value}")

    def __str__(self):
        return f"Report Field '{self.title}' ({self.key})"


# Field calculation functions

# Note that each function takes a single argument: a Report instance.
# So these functions could be implemented as Report instance methods.
# But this implementation could cause confusion because no one Report
# class can work with all these methods.

def calc_time_taken(report: Report) -> float:
    """ Calculate the time taken in minutes. """
    delta = report.get_field(TimeEndedF) - report.get_field(TimeBeganF)
    minutes = (delta.seconds + 1e-6 * delta.microseconds) / 60.
    if minutes < 0.:
        raise ValueError(f"Time taken must be positive, but got {minutes} min")
    return minutes


def calc_seqlen(report: Report):
    return len(report.get_field(SeqF))


def calc_n_batches(report: Report):
    return len(report.get_field(ChecksumsF))


def calc_speed(report: Report) -> float:
    nvecs = report.get_field(NumVecF)
    taken = report.get_field(TimeTakenF)
    try:
        return nvecs / taken
    except ZeroDivisionError:
        return inf if nvecs > 0 else nan


# Field value checking functions

def agrees(num1: float, num2: float, precision: int | float):
    """ Check if two floats agree within a given decimal precision. """
    return isclose(num1, num2, abs_tol=10. ** -precision, rel_tol=0.)


def nanagrees(num1: float, num2: float, precision: int | float):
    """ Like agrees, but also return True if both floats are NaN. """
    return (isnan(num1) and isnan(num2)) or agrees(num1, num2, precision)


def check_name(name: str):
    return bool(name) and not bool(set(name) - path.STR_CHARS_SET)


def check_nonneg_int(num: int):
    return isinstance(num, Integral) and num >= 0


def check_pos_int(num: int):
    return isinstance(num, Integral) and num > 0


def check_array_ints(arr: np.ndarray):
    return isinstance(arr, np.ndarray) and issubclass(arr.dtype.type, Integral)


def check_array_pos_ints(arr: np.ndarray):
    return check_array_ints(arr) and np.all(arr > 0)


def check_nonneg_float(num: float):
    return isinstance(num, float) and num >= 0.


def check_nannonneg_float(num: float):
    return isinstance(num, float) and num >= 0. or isnan(num)


def check_bool(x: bool):
    return isinstance(x, bool)


def check_float(num: float):
    return isinstance(num, float) and not isnan(num)


def check_pos_float(num: float):
    return isinstance(num, float) and num > 0.


def check_probability(probability: float):
    return isinstance(probability, float) and 0. <= probability <= 1.


def check_nanprobability(probability: float):
    return isinstance(probability, float) and (0. <= probability <= 1.
                                               or isnan(probability))


def check_percentage(percentage: float):
    return isinstance(percentage, float) and 0. <= percentage <= 100.


def check_nanpercentage(percentage: float):
    return isinstance(percentage, float) and (0. <= percentage <= 100.
                                              or isnan(percentage))


def check_seqlen(report: Report, seqlen: int):
    return check_pos_int(seqlen) and seqlen == len(report.get_field(SeqF))


def check_checksums(checksums: list[str]):
    checksum_pattern = re.compile("^[0-9A-Fa-f]{32}$")
    return (isinstance(checksums, list)
            and all(isinstance(checksum, str) for checksum in checksums)
            and all(map(checksum_pattern.match, checksums)))


def check_n_batches(report: Report, n_batches):
    return check_nonneg_int(n_batches) and n_batches == calc_n_batches(report)


def check_time_ended(report: Report, ended: datetime):
    return (isinstance(ended, datetime)
            and ended >= report.get_field(TimeBeganF))


def check_ints_range(ints: Iterable[int]):
    if not all(isinstance(num, int) for num in ints):
        return False
    ints_sort = sorted(ints)
    return ints_sort == list(range(ints_sort[0], ints_sort[-1] + 1))


def check_cluster_dict(cdict: dict[int, Any]):
    return (isinstance(cdict, dict)
            and check_ints_range(cdict.keys())
            and all(map(check_pos_int, cdict.keys())))


def check_clusts_list_nonneg_int(lis: dict[int, list[int]]):
    return (check_cluster_dict(lis)
            and all(isinstance(li, list) for li in lis.values())
            and all(all(map(check_nonneg_int, li)) for li in lis.values()))


def check_clusts_floats(floats: dict[int, float]):
    return (check_cluster_dict(floats)
            and all(map(check_float, floats.values())))


def check_clusts_list_floats(lfs: dict[int, list[float]]):
    return (check_cluster_dict(lfs)
            and all(isinstance(lf, list) for lf in lfs.values())
            and all(all(map(check_float, lf)) for lf in lfs.values()))


def check_list_str(lstr: list[str]):
    return isinstance(lstr, list) and all(isinstance(x, str) for x in lstr)


def check_dir(d: Path):
    return isinstance(d, Path) and d.is_dir()


def check_file(f: Path):
    return isinstance(f, Path) and f.is_file()


# Field definitions

DATETIME_FORMAT = "%Y-%m-%d at %H:%M:%S"
DECIMAL_PRECISION = 3  # general precision for decimals
PERC_VEC_PRECISION = 1
TIME_TAKEN_PRECISION = 2
SPEED_PRECISION = 0


def iconv_int_keys(mapping: dict[Any, Any]):
    return {int(key): value for key, value in mapping.items()}


def iconv_array_int(nums: list[int]):
    return np.asarray(nums, dtype=int)


def iconv_array_float(nums: list[float]):
    return np.asarray(nums, dtype=float)


def oconv_array_int(nums: np.ndarray):
    return list(map(int, nums))


@cache
def get_oconv_float(precision: int = DECIMAL_PRECISION):
    def oconv_float(num: float):
        return round(num, precision)

    return oconv_float


@cache
def get_oconv_list_float(precision: int = DECIMAL_PRECISION):
    oconv_float = get_oconv_float(precision)

    def oconv_list_float(nums: list[float]) -> list[float]:
        return list(map(oconv_float, nums))

    return oconv_list_float


@cache
def get_oconv_array_float(precision: int = DECIMAL_PRECISION):
    oconv_list_float = get_oconv_list_float(precision)

    def oconv_array_float(nums: np.array):
        return oconv_list_float(nums.tolist())

    return oconv_array_float


@cache
def get_oconv_dict_float(precision: int = DECIMAL_PRECISION):
    oconv_float = get_oconv_float(precision)

    def oconv_dict_float(dnum: dict[Hashable, float]) -> dict[Hashable, float]:
        return {d: oconv_float(num) for d, num in dnum.items()}

    return oconv_dict_float


@cache
def get_oconv_dict_list_float(precision: int = DECIMAL_PRECISION):
    oconv_list_float = get_oconv_list_float(precision)

    def oconv_dict_list_float(dnums: dict[Hashable, list[float]]
                              ) -> dict[Hashable, list[float]]:
        return {d: oconv_list_float(nums) for d, nums in dnums.items()}

    return oconv_dict_list_float


def iconv_datetime(text: str):
    return datetime.strptime(text, DATETIME_FORMAT)


def oconv_datetime(dtime: datetime):
    return dtime.strftime(DATETIME_FORMAT)


# General
OutDirF = Field("out_dir", "", Path, check_val=check_dir)
SampleF = Field("sample", "Name of Sample", str, check_val=check_name)
RefF = Field("ref", "Name of Reference", str, check_val=check_name)
SeqF = Field("seq", "Sequence of Reference", DNA,
             iconv=DNA.parse, oconv=DNA.__str__)
SectF = Field("sect", "Name of Section", str, check_val=check_name)
End5F = Field("end5", "5' end of Section", int, check_val=check_pos_int)
End3F = Field("end3", "3' end of Section", int, check_val=check_pos_int)
SeqLenF = Field("length", "Length of Sequence (nt)", int,
                check_rep_val=check_seqlen)
TimeBeganF = Field("began", "Time Began", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime)
TimeEndedF = Field("ended", "Time Ended", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime,
                   check_rep_val=check_time_ended)
TimeTakenF = Field("taken", "Time Taken (minutes)", float,
                   oconv=get_oconv_float(TIME_TAKEN_PRECISION))

# Relation vector generation
NumVecF = Field("n_reads_rel_pass", "Number of Reads Passed", int,
                check_val=check_nonneg_int)
NumErrF = Field("n_reads_rel_fail", "Number of Reads Failed", int,
                check_val=check_nonneg_int)
NumBatchF = Field("n_batches", "Number of Batches", int,
                  check_val=check_nonneg_int)
ChecksumsF = Field("checksums", "MD5 Checksums of Batches", list,
                   check_val=check_checksums)
SpeedF = Field("speed", "Speed (reads per minute)", float,
               oconv=get_oconv_float(SPEED_PRECISION))

# Mutation calling
CountMutsF = Field("count_muts", "Count the Following as Mutations",
                   SemiBitCaller,
                   iconv=SemiBitCaller.from_report_format,
                   oconv=SemiBitCaller.to_report_format)
CountRefsF = Field("count_refs", "Count the Following as Matches",
                   SemiBitCaller,
                   iconv=SemiBitCaller.from_report_format,
                   oconv=SemiBitCaller.to_report_format)

# Positional filtering
ExclPolyAF = Field("exclude_polya",
                   "Exclude Poly(A) Sequences of at Least This Length (nt)",
                   int, check_val=check_nonneg_int)
ExclGUF = Field("exclude_gu", "Exclude G/U Bases", bool)
ExclUserPosF = Field("exclude_pos", "Exclude User-Defined Positions",
                     np.ndarray,
                     iconv=iconv_array_int, oconv=oconv_array_int,
                     check_val=check_array_pos_ints)

# Bit vector filtering
MinInfoPosF = Field("min_ninfo_pos",
                    "Minimum Number of Informative Reads per Position",
                    int, check_val=check_nonneg_int)
MinMutPosF = Field("min_fmut_pos",
                   "Minimum Fraction of Mutations per Position",
                   float, oconv=get_oconv_float(),
                   check_val=check_probability)
MaxMutPosF = Field("max_fmut_pos",
                   "Maximum Fraction of Mutations per Position",
                   float, oconv=get_oconv_float(),
                   check_val=check_probability)
MinMutGapF = Field("min_mut_gap",
                   "Minimum Gap Between Mutations (nt)",
                   int, check_val=check_nonneg_int)
MinInfoReadF = Field("min_finfo_read",
                     "Minimum Fraction of Informative Positions per Read",
                     float, oconv=get_oconv_float(),
                     check_val=check_probability)
MaxMutReadF = Field("max_fmut_read",
                    "Maximum Fraction of Mutations per Read",
                    float, oconv=get_oconv_float(),
                    check_val=check_probability)
PosInitF = Field("pos_init",
                 "Positions Initially Given",
                 np.ndarray,
                 iconv=iconv_array_int, oconv=oconv_array_int,
                 check_val=check_array_pos_ints)
PosCutPolyAF = Field("pos_polya",
                     "Positions Cut -- Poly(A) Sequence",
                     np.ndarray,
                     iconv=iconv_array_int, oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosCutGUF = Field("pos_gu",
                  "Positions Cut -- G/U Base",
                  np.ndarray,
                  iconv=iconv_array_int, oconv=oconv_array_int,
                  check_val=check_array_pos_ints)
PosCutUserF = Field("pos_user",
                    "Positions Cut -- User-Specified",
                    np.ndarray,
                    iconv=iconv_array_int, oconv=oconv_array_int,
                    check_val=check_array_pos_ints)
PosCutLoInfoF = Field("pos_min_ninfo",
                      "Positions Cut -- Too Few Informative Reads",
                      np.ndarray,
                      iconv=iconv_array_int, oconv=oconv_array_int,
                      check_val=check_array_pos_ints)
PosCutLoMutF = Field("pos_min_fmut",
                     "Positions Cut -- Too Few Mutations",
                     np.ndarray,
                     iconv=iconv_array_int, oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosCutHiMutF = Field("pos_max_fmut",
                     "Positions Cut -- Too Many Mutations",
                     np.ndarray,
                     iconv=iconv_array_int, oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosKeptF = Field("pos_kept",
                 "Positions Ultimately Kept",
                 np.ndarray,
                 iconv=iconv_array_int, oconv=oconv_array_int,
                 check_val=check_array_pos_ints)
NumPosInitF = Field("n_pos_init",
                    "Number of Positions Initially Given",
                    int, check_val=check_nonneg_int)
NumPosCutPolyAF = Field("n_pos_polya",
                        "Number of Positions Cut -- Poly(A) Sequence",
                        int, check_val=check_nonneg_int)
NumPosCutGUF = Field("n_pos_gu",
                     "Number of Positions Cut -- G/U Base",
                     int, check_val=check_nonneg_int)
NumPosCutUserF = Field("n_pos_user",
                       "Number of Positions Cut -- User-Specified",
                       int, check_val=check_nonneg_int)
NumPosCutLoInfoF = Field("n_pos_min_ninfo",
                         "Number of Positions Cut -- Too Few Informative Reads",
                         int, check_val=check_nonneg_int)
NumPosCutLoMutF = Field("n_pos_min_fmut",
                        "Number of Positions Cut -- Too Few Mutations",
                        int, check_val=check_nonneg_int)
NumPosCutHiMutF = Field("n_pos_max_fmut",
                        "Number of Positions Cut -- Too Many Mutations",
                        int, check_val=check_nonneg_int)
NumPosKeptF = Field("n_pos_kept",
                    "Number of Positions Ultimately Kept",
                    int, check_val=check_nonneg_int)
NumReadsInitF = Field("n_reads_init",
                      "Number of Reads Initially Given",
                      int, check_val=check_nonneg_int)
NumReadsLoInfoF = Field("n_reads_min_finfo",
                        "Number of Reads Cut -- Too Few Informative Positions",
                        int, check_val=check_nonneg_int)
NumReadsHiMutF = Field("n_reads_max_fmut",
                       "Number of Reads Cut -- Too Many Mutations",
                       int, check_val=check_nonneg_int)
NumReadsCloseMutF = Field("n_reads_min_gap",
                          "Number of Reads Cut -- Mutations Too Close Together",
                          int, check_val=check_nonneg_int)
NumReadsKeptF = Field("n_reads_kept",
                      "Number of Reads Ultimately Kept",
                      int, check_val=check_nonneg_int)
NumUniqReadKeptF = Field("n_uniq_reads",
                         "Number of Unique Bit Vectors",
                         int, check_val=check_nonneg_int)

# EM clustering
MinIterClustF = Field("min_iter",
                      "Minimum EM Iterations per Cluster",
                      int, check_val=check_nonneg_int)
MaxIterClustF = Field("max_iter",
                      "Maximum EM Iterations per Cluster",
                      int, check_val=check_pos_int)
ClustConvThreshF = Field("conv_thresh",
                         "Convergence Threshold for Log Likelihood",
                         float, oconv=get_oconv_float(),
                         check_val=check_pos_float)
MaxClustsF = Field("max_order",
                   "Maximum Number of Clusters",
                   int, check_val=check_pos_int)
ClustNumRunsF = Field("num_runs",
                      "Number of Independent EM Runs",
                      int, check_val=check_pos_int)
NumClustsF = Field("best_order",
                   "Optimal Number of Clusters",
                   int, check_val=check_pos_int)
ClustsBicF = Field("bic",
                   "Bayesian Information Criterion per Order",
                   dict, iconv=iconv_int_keys, oconv=get_oconv_dict_float(),
                   check_val=check_clusts_floats)
ClustsConvF = Field("converged",
                    "Iterations Until Convergence per Run",
                    dict, iconv=iconv_int_keys,
                    check_val=check_clusts_list_nonneg_int)
ClustsLogLikesF = Field("log_likes",
                        "Log Likelihood per Run",
                        dict, iconv=iconv_int_keys,
                        oconv=get_oconv_dict_list_float(),
                        check_val=check_clusts_list_floats)
ClustsLikeMeanF = Field("log_like_mean",
                        "Mean Log Likelihood per Order",
                        dict, iconv=iconv_int_keys,
                        oconv=get_oconv_dict_float(),
                        check_val=check_clusts_floats)
ClustsLikeStdF = Field("log_like_std",
                       "Std. Dev. Log Likelihood per Order",
                       dict, iconv=iconv_int_keys, oconv=get_oconv_dict_float(),
                       check_val=check_clusts_floats)
ClustsVarInfoF = Field("var_info",
                       "Variation of Information per Order",
                       dict, iconv=iconv_int_keys, oconv=get_oconv_dict_float(),
                       check_val=check_clusts_floats)


# Field managing functions

@cache
def fields() -> list[Field]:
    return [member for _, member in getmembers(sys.modules[__name__])
            if isinstance(member, Field)]


@cache
def field_keys() -> dict[str, Field]:
    return {field.key: field for field in fields()}


@cache
def field_titles() -> dict[str, Field]:
    return {field.title: field for field in fields() if field.title}


def lookup_key(key: str) -> Field:
    try:
        return field_keys()[key]
    except KeyError:
        # Suppress exception chaining.
        pass
    raise ValueError(f"No Report Field is keyed '{key}'")


def lookup_title(title: str) -> Field:
    if not title:
        raise ValueError("Got blank title for field")
    try:
        return field_titles()[title]
    except KeyError:
        # Suppress exception chaining.
        pass
    raise ValueError(f"No Report Field is titled '{title}'")


# Report classes

class Report(ABC):
    """ Abstract base class for a report from a step. """
    __slots__ = "out_dir",

    def __init__(self, **kwargs: Any | Callable[[Report], Any]):
        # Get the values of all attributes (specified in __slots__) from
        # the keyword arguments.
        for key, val in kwargs.items():
            if callable(val):
                # If the value of the keyword argument is callable, then
                # it must accept one argument -- self -- and return the
                # value of the attribute.
                val = val(self)
            self.__setattr__(key, val)
        logger.debug(f"Created new {self.__class__.__name__} with fields "
                     f"from keywords {kwargs}")

    @classmethod
    @abstractmethod
    def path_segs(cls):
        """ Return a tuple of the segments of the path. """
        return path.ModSeg, path.SampSeg, path.RefSeg

    @classmethod
    def auto_fields(cls) -> dict[str, Any]:
        """ Return the fields that are filled automatically. """
        return {path.EXT: path.JSON_EXT}

    @classmethod
    def build_path(cls, out_dir: Path, **path_fields):
        """ Build a path for a report from the given fields. """
        return path.buildpar(*cls.path_segs(), top=out_dir,
                             **{**cls.auto_fields(), **path_fields})

    def path_fields(self) -> dict[str, Any]:
        """ Return a dict of the fields of the path that are attributes
        of the report. """
        return {key: self.__getattribute__(key)
                for segment in self.path_segs()
                for key, field in segment.field_types.items()
                if hasattr(self, key)}

    def get_path(self):
        """ Return the path of the report. """
        return self.build_path(self.out_dir, **self.path_fields())

    def get_field(self, field: Field):
        """ Return the value of a field of the report using the field
        instance directly, not its key. """
        return self.__getattribute__(field.key)

    def to_dict(self):
        """ Return a dict of raw values of the fields, keyed by the
        titles of their fields. """
        odata = dict()
        for key in self.__slots__:
            field = lookup_key(key)
            if field.title:
                # Output only the fields with non-blank titles.
                value = self.get_field(field)
                if field.oconv is not None:
                    # Convert the value to the proper output value.
                    value = field.oconv(value)
                odata[field.title] = value
        return odata

    def save(self):
        with open(self.get_path(), "w") as f:
            json.dump(self.to_dict(), f, indent=4)
        logger.info(f"Wrote {self} to {self.get_path()}")

    @classmethod
    def from_dict(cls, out_dir: Path, odata: dict[str, Any]):
        """ Convert a dict of raw values (keyed by the titles of their
        fields) into a dict of encoded values (keyed by the keys of
        their fields), from which a new Report is instantiated. """
        if not isinstance(odata, dict):
            raise TypeError("Report classmethod from_data expected 'dict', "
                            f"but got '{type(odata).__name__}'")
        # Read every raw value, keyed by the title of its field.
        idata = {OutDirF.key: out_dir}
        for title, value in odata.items():
            # Get the field corresponding to the title.
            field = lookup_title(title)
            # Cast the value to the input type and key it by the field.
            idata[field.key] = field.iconv(value)
        # Instantiate and return a new Report from the values.
        return cls(**idata)

    @classmethod
    def open(cls, json_file: Path):
        """ Create a new Report instance from a JSON file. """
        path_fields = path.parse(json_file, *cls.path_segs())
        out_dir = path_fields[path.TOP]
        with open(json_file) as f:
            report: Report | BatchReport = cls.from_dict(out_dir, json.load(f))
        # Ensure that the path-related fields in the JSON data match the
        # actual path of the JSON file.
        for key, value in report.path_fields().items():
            if value != path_fields[key]:
                raise ValueError(f"Got different values for field '{key}' from "
                                 f"path ({path_fields[key]}) and contents "
                                 f"({value}) of report {json_file}")
        return report

    def __setattr__(self, key, value):
        """ Validate the attribute name and value before setting it. """
        lookup_key(key).validate(self, value)
        super().__setattr__(key, value)

    def __str__(self):
        descript = ", ".join(f"{key} = {repr(val)}"
                             for key, val in self.to_dict().items()
                             if isinstance(val, str))
        return f"{self.__class__.__name__}: {descript}"

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.to_dict() == other.to_dict()


class BatchReport(Report, ABC):

    @classmethod
    def open(cls, json_file: Path):
        # Open the report.
        report: BatchReport = super().open(json_file)
        # Verify that the batches exist and have the right checksums.
        missing, badsum = report.find_invalid_batches()
        if missing:
            files = list(map(report.get_batch_path, missing))
            raise FileNotFoundError(f"Missing batches: {files}")
        if badsum:
            files = list(map(report.get_batch_path, badsum))
            raise ValueError(f"Invalid MD5 checksums: {files}")
        return report

    @property
    def num_batches(self):
        """ Number of batches associated with the report. """
        return self.get_field(NumBatchF)

    @classmethod
    @abstractmethod
    def get_batch_seg(cls):
        """ Type of the path segment of each batch file. """
        return path.Segment("", dict())

    @classmethod
    def build_batch_path(cls, out_dir: Path, batch: int, ext: str | None = None,
                         **path_fields):
        """ Get the path to a batch of a specific number. """
        seg = cls.get_batch_seg()
        if ext is None:
            # Determine if a file with a valid extension already exists.
            for valid_ext in seg.exts:
                valid_file = cls.build_batch_path(out_dir, batch, valid_ext,
                                                  **path_fields)
                if valid_file.is_file():
                    # If such a file already exists, then return it.
                    return valid_file
            # If not, then return a path with the first file extension
            # listed in seg.exts.
            ext = seg.exts[0]
        # Determine the path to the parent directory of the batch.
        dir_path = cls.build_path(out_dir, **path_fields).parent
        # Return the batch segment inside that parent directory.
        return dir_path.joinpath(seg.build(batch=batch, ext=ext))

    def get_batch_path(self, batch: int):
        """ Get the path to a batch of a specific number. """
        return self.build_batch_path(self.out_dir, batch, **self.path_fields())

    def digest_batch(self, batch: int):
        """ Compute the MD5 digest of a specific batch. """
        return digest_file(self.get_batch_path(batch))

    def iter_batch_paths(self):
        """ Iterate through all batch paths. """
        yield from map(self.get_batch_path, range(self.get_field(NumBatchF)))

    def find_invalid_batches(self):
        """ Return all the batches of mutation vectors that either do
        not exist or do not match their expected checksums. """
        missing: list[int] = list()
        badsum: list[int] = list()
        for batch, checksum in enumerate(self.get_field(ChecksumsF)):
            try:
                if checksum != self.digest_batch(batch):
                    # The file exists does not match the checksum.
                    badsum.append(batch)
            except FileNotFoundError:
                # The batch file does not exist.
                missing.append(batch)
        return missing, badsum
