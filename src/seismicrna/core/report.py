from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from inspect import getmembers
import json
from functools import cache, partial
from logging import getLogger
from math import isclose, isnan, inf, nan
from numbers import Integral
from pathlib import Path
import re
import sys
from typing import Any, Hashable, Callable, Iterable

import numpy as np

from . import path
from .batch import BaseBatch
from .bitcall import SemiBitCaller
from .files import digest_file
from .output import Output
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

def calc_time_taken(report: Report):
    """ Calculate the time taken in minutes. """
    delta = report.get_field(TimeEndedF) - report.get_field(TimeBeganF)
    minutes = (delta.seconds + 1e-6 * delta.microseconds) / 60.
    if minutes < 0.:
        raise ValueError(f"Time taken must be positive, but got {minutes} min")
    return minutes


def calc_seqlen(report: Report):
    return len(report.get_field(SeqF))


def calc_speed(report: Report):
    nvecs = report.get_field(NumReadsRel)
    taken = report.get_field(TimeTakenF)
    try:
        return nvecs / taken
    except ZeroDivisionError:
        return inf if nvecs > 0 else nan


def calc_checksums(report: BatchReport, top: Path):
    """ Calculate the checksums. """
    return report.compute_checksums(top)


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


def check_dict_vals_nonneg_ints(mapping: dict[Any, int]):
    return isinstance(mapping, dict) and all(map(check_nonneg_int,
                                                 mapping.values()))


def check_nonneg_float(num: float):
    return isinstance(num, float) and num >= 0.


def check_nannonneg_float(num: float):
    return isinstance(num, float) and num >= 0. or isnan(num)


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


def check_checksums_btype(checksums: list[str]):
    checksum_pattern = re.compile("^[0-9A-Fa-f]{32}$")
    return (isinstance(checksums, list)
            and all(isinstance(checksum, str) for checksum in checksums)
            and all(map(checksum_pattern.match, checksums)))


def check_checksums(checksums: dict[str, list[str]]):
    return (isinstance(checksums, dict)
            and all(isinstance(btype, str) for btype in checksums.keys())
            and all(map(check_checksums_btype, checksums.values())))


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


def iconv_dict_str_int(mapping: dict[Any, Any]) -> dict[str, int]:
    return {str(key): int(value) for key, value in mapping.items()}


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
SampleF = Field("sample", "Name of Sample", str, check_val=check_name)
RefF = Field("ref", "Name of Reference", str, check_val=check_name)
SeqF = Field("seq", "Sequence of Reference", DNA, oconv=str)
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

# Alignment
IsDemultiplexed = Field("demultiplexed", "Use demultiplexed mode", bool)
IsPairedEnd = Field("paired_end", "Use paired-end mode", bool)
PhredEnc = Field("phred_enc", "Phred score encoding", int)
FastqcRun = Field("fastqc", "Check quality with FastQC", bool)
CutadaptRun = Field("cut", "Trim with Cutadapt", bool)
CutadaptQ1 = Field("cut_q1", "Minimum Phred score for read 1",
                   int, check_val=check_nonneg_int)
CutadaptQ2 = Field("cut_q2", "Minimum Phred score for read 2",
                   int, check_val=check_nonneg_int)
CutadaptG1 = Field("cut_g1", "5' adapter for read 1",
                   list, check_val=check_list_str)
CutadaptA1 = Field("cut_a1", "3' adapter for read 1",
                   list, check_val=check_list_str)
CutadaptG2 = Field("cut_g2", "5' adapter for read 2",
                   list, check_val=check_list_str)
CutadaptA2 = Field("cut_a2", "3' adapter for read 2",
                   list, check_val=check_list_str)
CutadaptOverlap = Field("cut_o", "Minimum adapter length (nt)",
                        int, check_val=check_nonneg_int)
CutadaptErrors = Field("cut_e", "Minimum adapter error rate",
                       float, check_val=check_probability)
CutadaptIndels = Field("cut_indels", "Allow indels in adapters", bool)
CutadaptNextSeq = Field("cut_nextseq", "Trim in NextSeq mode", bool)
CutadaptNoTrimmed = Field("cut_discard_trimmed", "Discard trimmed reads",
                          bool)
CutadaptNoUntrimmed = Field("cut_discard_untrimmed", "Discard untrimmed reads",
                            bool)
CutadaptMinLength = Field("cut_m", "Minimum length of read after trimming (nt)",
                          int, check_val=check_nonneg_int)
Bowtie2Local = Field("bt2_local", "Align in local mode", bool)
Bowtie2Discord = Field("bt2_discordant", "Keep discordant alignments", bool)
Bowtie2Dovetail = Field("bt2_dovetail",
                        "Consider dovetailed alignments concordant",
                        bool)
Bowtie2Contain = Field("bt2_contain",
                       "Consider nested alignments concordant",
                       bool)
Bowtie2Mixed = Field("bt2_mixed", "Align in mixed mode", bool)
Bowtie2Unal = Field("bt2_unal", "Output unaligned reads to SAM file", bool)
Bowtie2ScoreMin = Field("bt2_score_min", "Minimum alignment score", str)
Bowtie2MinLength = Field("bt2_i", "Minimum fragment length (nt)",
                         int, check_val=check_nonneg_int)
Bowtie2MaxLength = Field("bt2_x", "Maximum fragment length (nt)",
                         int, check_val=check_nonneg_int)
Bowtie2GBar = Field("bt2_gbar", "Gap buffer margin (nt)",
                    int, check_val=check_nonneg_int)
Bowtie2SeedLength = Field("bt2_l", "Seed length (nt)",
                          int, check_val=check_nonneg_int)
Bowtie2SeedInterval = Field("bt2_s", "Seed interval (nt)", str)
Bowtie2ExtTries = Field("bt2_d", "Maximum seed extension attempts",
                        int, check_val=check_nonneg_int)
Bowtie2Reseed = Field("bt2_r", "Maximum re-seeding attempts",
                      int, check_val=check_nonneg_int)
Bowtie2Dpad = Field("bt2_dpad", "Dynamic programming padding (nt)",
                    int, check_val=check_nonneg_int)
Bowtie2Orient = Field("bt2_orient", "Orientation of mates", str)
MinMapQual = Field("min_mapq", "Minimum mapping quality",
                   int, check_val=check_nonneg_int)
ReadsInit = Field("reads_init", "Number of reads initially",
                  int, check_val=check_nonneg_int)
ReadsTrim = Field("reads_trim", "Number of reads after trimming",
                  int, check_val=check_nonneg_int)
ReadsAlign = Field("reads_align", "Number of reads after alignment",
                   dict, iconv=iconv_dict_str_int,
                   check_val=check_dict_vals_nonneg_ints)
ReadsDedup = Field("reads_filter", "Number of reads after filtering",
                   dict, iconv=iconv_dict_str_int,
                   check_val=check_dict_vals_nonneg_ints)
ReadsRefs = Field("reads_refs", "Number of reads aligned by reference",
                  dict, iconv=iconv_dict_str_int,
                  check_val=check_dict_vals_nonneg_ints)

# Relation vector generation
NumReadsRel = Field("n_reads_rel", "Number of Reads", int,
                    check_val=check_nonneg_int)
NumBatchF = Field("n_batches", "Number of Batches", int,
                  check_val=check_nonneg_int)
ChecksumsF = Field("checksums", "MD5 Checksums of Batches", dict,
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
                   int,
                   check_val=check_nonneg_int)
ExclGUF = Field("exclude_gu", "Exclude G/U Bases", bool)
ExclUserPosF = Field("exclude_pos",
                     "Exclude User-Defined Positions",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int,
                     check_val=check_array_pos_ints)

# Bit vector filtering
MinInfoPosF = Field("min_ninfo_pos",
                    "Minimum Number of Informative Reads per Position",
                    int,
                    check_val=check_nonneg_int)
MinMutPosF = Field("min_fmut_pos",
                   "Minimum Fraction of Mutations per Position",
                   float,
                   oconv=get_oconv_float(),
                   check_val=check_probability)
MaxMutPosF = Field("max_fmut_pos",
                   "Maximum Fraction of Mutations per Position",
                   float,
                   oconv=get_oconv_float(),
                   check_val=check_probability)
MinMutGapF = Field("min_mut_gap",
                   "Minimum Gap Between Mutations (nt)",
                   int,
                   check_val=check_nonneg_int)
MinInfoReadF = Field("min_finfo_read",
                     "Minimum Fraction of Informative Positions per Read",
                     float,
                     oconv=get_oconv_float(),
                     check_val=check_probability)
MaxMutReadF = Field("max_fmut_read",
                    "Maximum Fraction of Mutations per Read",
                    float,
                    oconv=get_oconv_float(),
                    check_val=check_probability)
PosInitF = Field("pos_init",
                 "Positions Initially Given",
                 np.ndarray,
                 iconv=iconv_array_int,
                 oconv=oconv_array_int,
                 check_val=check_array_pos_ints)
PosCutPolyAF = Field("pos_polya",
                     "Positions Cut -- Poly(A) Sequence",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosCutGUF = Field("pos_gu",
                  "Positions Cut -- G/U Base",
                  np.ndarray,
                  iconv=iconv_array_int,
                  oconv=oconv_array_int,
                  check_val=check_array_pos_ints)
PosCutUserF = Field("pos_user",
                    "Positions Cut -- User-Specified",
                    np.ndarray,
                    iconv=iconv_array_int,
                    oconv=oconv_array_int,
                    check_val=check_array_pos_ints)
PosCutLoInfoF = Field("pos_min_ninfo",
                      "Positions Cut -- Too Few Informative Reads",
                      np.ndarray,
                      iconv=iconv_array_int,
                      oconv=oconv_array_int,
                      check_val=check_array_pos_ints)
PosCutLoMutF = Field("pos_min_fmut",
                     "Positions Cut -- Too Few Mutations",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosCutHiMutF = Field("pos_max_fmut",
                     "Positions Cut -- Too Many Mutations",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int,
                     check_val=check_array_pos_ints)
PosKeptF = Field("pos_kept",
                 "Positions Ultimately Kept",
                 np.ndarray,
                 iconv=iconv_array_int,
                 oconv=oconv_array_int,
                 check_val=check_array_pos_ints)
NumPosInitF = Field("n_pos_init",
                    "Number of Positions Initially Given",
                    int,
                    check_val=check_nonneg_int)
NumPosCutPolyAF = Field("n_pos_polya",
                        "Number of Positions Cut -- Poly(A) Sequence",
                        int,
                        check_val=check_nonneg_int)
NumPosCutGUF = Field("n_pos_gu",
                     "Number of Positions Cut -- G/U Base",
                     int,
                     check_val=check_nonneg_int)
NumPosCutUserF = Field("n_pos_user",
                       "Number of Positions Cut -- User-Specified",
                       int,
                       check_val=check_nonneg_int)
NumPosCutLoInfoF = Field("n_pos_min_ninfo",
                         "Number of Positions Cut -- Too Few Informative Reads",
                         int,
                         check_val=check_nonneg_int)
NumPosCutLoMutF = Field("n_pos_min_fmut",
                        "Number of Positions Cut -- Too Few Mutations",
                        int,
                        check_val=check_nonneg_int)
NumPosCutHiMutF = Field("n_pos_max_fmut",
                        "Number of Positions Cut -- Too Many Mutations",
                        int,
                        check_val=check_nonneg_int)
NumPosKeptF = Field("n_pos_kept",
                    "Number of Positions Ultimately Kept",
                    int,
                    check_val=check_nonneg_int)
NumReadsInitF = Field("n_reads_init",
                      "Number of Reads Initially Given",
                      int,
                      check_val=check_nonneg_int)
NumReadsLoInfoF = Field("n_reads_min_finfo",
                        "Number of Reads Cut -- Too Few Informative Positions",
                        int,
                        check_val=check_nonneg_int)
NumReadsHiMutF = Field("n_reads_max_fmut",
                       "Number of Reads Cut -- Too Many Mutations",
                       int,
                       check_val=check_nonneg_int)
NumReadsCloseMutF = Field("n_reads_min_gap",
                          "Number of Reads Cut -- Mutations Too Close Together",
                          int,
                          check_val=check_nonneg_int)
NumReadsKeptF = Field("n_reads_kept",
                      "Number of Reads Ultimately Kept",
                      int,
                      check_val=check_nonneg_int)
NumUniqReadKeptF = Field("n_uniq_reads",
                         "Number of Unique Bit Vectors",
                         int,
                         check_val=check_nonneg_int)

# EM clustering
MinIterClustF = Field("min_iter",
                      "Minimum EM Iterations per Cluster",
                      int,
                      check_val=check_nonneg_int)
MaxIterClustF = Field("max_iter",
                      "Maximum EM Iterations per Cluster",
                      int,
                      check_val=check_pos_int)
ClustConvThreshF = Field("conv_thresh",
                         "Convergence Threshold for Log Likelihood",
                         float,
                         oconv=get_oconv_float(),
                         check_val=check_pos_float)
MaxClustsF = Field("max_order",
                   "Maximum Number of Clusters",
                   int,
                   check_val=check_pos_int)
ClustNumRunsF = Field("num_runs",
                      "Number of Independent EM Runs",
                      int, check_val=check_pos_int)
NumClustsF = Field("best_order",
                   "Optimal Number of Clusters",
                   int,
                   check_val=check_pos_int)
ClustsBicF = Field("bic",
                   "Bayesian Information Criterion per Order",
                   dict,
                   iconv=iconv_int_keys,
                   oconv=get_oconv_dict_float(),
                   check_val=check_clusts_floats)
ClustsConvF = Field("converged",
                    "Iterations Until Convergence per Run",
                    dict,
                    iconv=iconv_int_keys,
                    check_val=check_clusts_list_nonneg_int)
ClustsLogLikesF = Field("log_likes",
                        "Log Likelihood per Run",
                        dict,
                        iconv=iconv_int_keys,
                        oconv=get_oconv_dict_list_float(),
                        check_val=check_clusts_list_floats)
ClustsLikeMeanF = Field("log_like_mean",
                        "Mean Log Likelihood per Order",
                        dict,
                        iconv=iconv_int_keys,
                        oconv=get_oconv_dict_float(),
                        check_val=check_clusts_floats)
ClustsLikeStdF = Field("log_like_std",
                       "Std. Dev. Log Likelihood per Order",
                       dict,
                       iconv=iconv_int_keys,
                       oconv=get_oconv_dict_float(),
                       check_val=check_clusts_floats)
ClustsVarInfoF = Field("var_info",
                       "Variation of Information per Order",
                       dict,
                       iconv=iconv_int_keys,
                       oconv=get_oconv_dict_float(),
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
        raise ValueError(f"Invalid report key: {repr(key)}")


def lookup_title(title: str) -> Field:
    if not title:
        raise ValueError("Got blank title for field")
    try:
        return field_titles()[title]
    except KeyError:
        raise ValueError(f"Invalid report field: {repr(title)}")


# Report classes

class Report(Output, ABC):
    """ Abstract base class for a report from a step. """

    @classmethod
    @abstractmethod
    def field_names(cls) -> tuple[str, ...]:
        """ Names of all fields of the report. """
        return tuple()

    @classmethod
    def from_dict(cls, odata: dict[str, Any]):
        """ Convert a dict of raw values (keyed by the titles of their
        fields) into a dict of encoded values (keyed by the keys of
        their fields), from which a new Report is instantiated. """
        if not isinstance(odata, dict):
            raise TypeError("Report classmethod from_data expected 'dict', "
                            f"but got {repr(type(odata).__name__)}")
        # Read every raw value, keyed by the title of its field.
        idata = dict()
        for title, value in odata.items():
            # Get the field corresponding to the title.
            field = lookup_title(title)
            # Cast the value to the input type and key it by the field.
            idata[field.key] = field.iconv(value)
        # Instantiate and return a new Report from the values.
        return cls(**idata)

    @classmethod
    def load(cls, file: Path):
        path_fields = path.parse(file, *cls.seg_types())
        with open(file) as f:
            report: Report | BatchReport = cls.from_dict(json.load(f))
        # Ensure that the path-related fields in the JSON data match the
        # actual path of the JSON file.
        for key, value in report.path_fields(path_fields[path.TOP]).items():
            if value != path_fields[key]:
                raise ValueError(f"Got different values for field {repr(key)} "
                                 f"from path ({repr(path_fields[key])}) and "
                                 f"contents ({repr(value)}) of report {file}")
        return report

    @classmethod
    def auto_fields(cls) -> dict[str, Any]:
        return {path.EXT: path.JSON_EXT}

    def __init__(self, **kwargs: Any | Callable[[Report], Any]):
        for name in self.field_names():
            value = kwargs.pop(name)
            if callable(value):
                # If the value of the keyword argument is callable, then
                # it must accept one argument -- self -- and return the
                # value of the attribute.
                value = value(self)
            self.__setattr__(name, value)
        if kwargs:
            raise ValueError(f"Invalid keywords for {self.__class__.__name__}: "
                             f"{list(kwargs)}")

    def get_field(self, field: Field):
        """ Return the value of a field of the report using the field
        instance directly, not its key. """
        return getattr(self, field.key)

    def to_dict(self):
        """ Return a dict of raw values of the fields, keyed by the
        titles of their fields. """
        odata = dict()
        for key in self.field_names():
            field = lookup_key(key)
            # Output only the fields with non-blank titles.
            value = self.get_field(field)
            if field.oconv is not None:
                # Convert the value to the proper output value.
                value = field.oconv(value)
            odata[field.title] = value
        return odata

    def save(self, top: Path, overwrite: bool = False):
        """ Save the report to a JSON file. """
        text = json.dumps(self.to_dict(), indent=4)
        save_path = self.get_path(top)
        with open(save_path, 'w' if overwrite else 'x') as f:
            f.write(text)
        logger.info(f"Wrote {self} to {save_path}")
        return save_path

    def __setattr__(self, key: str, value: Any):
        """ Validate the attribute name and value before setting it. """
        if key not in self.field_names():
            raise ValueError(
                f"Invalid field for {self.__class__.__name__}: {repr(key)}")
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
    @abstractmethod
    def field_names(cls):
        return super().field_names() + ("n_batches", "checksums")

    @classmethod
    @abstractmethod
    def _batch_types(cls) -> tuple[type[BaseBatch], ...]:
        """ Type(s) of batch(es) for the report. """

    @classmethod
    @cache
    def batch_types(cls) -> dict[str, type[BaseBatch]]:
        """ Type(s) of batch(es) for the report, keyed by name. """
        return {batch_type.btype(): batch_type
                for batch_type in cls._batch_types()}

    @classmethod
    def get_batch_type(cls, btype: str | None = None) -> type[BaseBatch]:
        """ Return a valid type of batch based on its name. """
        if btype is None:
            if (ntypes := len(cls.batch_types())) != 1:
                raise ValueError(f"btype is optional only if there is exactly "
                                 f"1 type of batch, but got {ntypes} types")
            return list(cls.batch_types().values())[0]
        return cls.batch_types()[btype]

    def __init__(self, *,
                 top: Path | None = None,
                 checksums: dict | Callable = calc_checksums,
                 **kwargs):
        super().__init__(**kwargs,
                         checksums=(partial(checksums, top=top)
                                    if callable(checksums)
                                    else checksums))

    @property
    def batch_nums(self):
        """ Batch numbers. """
        return range(self.get_field(NumBatchF))

    def get_batch_path(self, top: Path, batch: int, btype: str | None = None):
        """ Get the path to a batch of a specific number. """
        bt = self.get_batch_type(btype)
        return bt.build_path(batch=batch,
                             **self.path_fields(top, bt.auto_fields()))

    def digest_batch(self, top: Path, batch: int, btype: str | None = None):
        """ Compute the checksum of a specific batch. """
        return digest_file(self.get_batch_path(top, batch, btype))

    def compute_checksums(self, top: Path):
        """ Compute all checksums. """
        return {btype: [self.digest_batch(top, batch, btype)
                        for batch in self.batch_nums]
                for btype in self.batch_types()}

    def iter_batch_paths(self, top: Path, btype: str | None = None):
        """ Iterate through all batch paths. """
        for batch in self.batch_nums:
            yield self.get_batch_path(top, batch, btype)

    def find_invalid_batches(self, top: Path):
        """ Return all the batches of relation vectors that either do
        not exist or do not match their expected checksums. """
        missing: list[tuple[str, int]] = list()
        badsum: list[tuple[str, int]] = list()
        for btype, checksums in self.get_field(ChecksumsF).items():
            for batch, checksum in enumerate(checksums):
                label = btype, batch
                try:
                    if checksum != self.digest_batch(top, batch, btype):
                        # The file exists does not match the checksum.
                        badsum.append(label)
                except FileNotFoundError:
                    # The batch file does not exist.
                    missing.append(label)
        return missing, badsum

    def verify_batches(self, top: Path):
        missing, badsum = self.find_invalid_batches(top)
        if missing_files := [self.get_batch_path(top, batch, btype)
                             for btype, batch in missing]:
            raise FileNotFoundError(f"Missing batches: {missing_files}")
        if badsum_files := [self.get_batch_path(top, batch, btype)
                            for btype, batch in badsum]:
            raise ValueError(f"Invalid checksums: {badsum_files}")

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
