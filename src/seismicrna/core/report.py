from __future__ import annotations

import json
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from functools import cache
from inspect import getmembers
from logging import getLogger
from pathlib import Path
from typing import Any, Hashable, Callable

import numpy as np
from click import Option

from .arg import (opt_phred_enc,
                  opt_fastqc,
                  opt_cutadapt,
                  opt_cut_q1,
                  opt_cut_q2,
                  opt_cut_g1,
                  opt_cut_g2,
                  opt_cut_a1,
                  opt_cut_a2,
                  opt_cut_discard_trimmed,
                  opt_cut_discard_untrimmed,
                  opt_cut_o,
                  opt_cut_e,
                  opt_cut_indels,
                  opt_cut_nextseq,
                  opt_cut_m,
                  opt_bt2_d,
                  opt_bt2_r,
                  opt_bt2_dpad,
                  opt_bt2_orient,
                  opt_bt2_i,
                  opt_bt2_x,
                  opt_bt2_s,
                  opt_bt2_l,
                  opt_bt2_gbar,
                  opt_bt2_un,
                  opt_bt2_discordant,
                  opt_bt2_mixed,
                  opt_bt2_dovetail,
                  opt_bt2_contain,
                  opt_bt2_local,
                  opt_cram,
                  opt_min_reads,
                  opt_min_mapq,
                  opt_fold_temp,
                  opt_fold_md,
                  opt_fold_mfe,
                  opt_fold_max,
                  opt_fold_percent,
                  opt_quantile)
from .io import FileIO, ReadBatchIO, RefIO
from .rel import HalfRelPattern
from .version import __version__
from .write import need_write, write_mode

logger = getLogger(__name__)


# Field class

class Field(object):
    """ Field of a report. """

    __slots__ = "key", "title", "dtype", "default", "iconv", "oconv"

    def __init__(self,
                 key: str,
                 title: str,
                 dtype: type,
                 default: Any | None = None, *,
                 iconv: Callable[[Any], Any] | None = None,
                 oconv: Callable[[Any], Any] | None = None):
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
        default: Any
            Default value of the field
        iconv: Callable[[Any], Any] | None = None
            Convert an input value to the right type for instantiation.
            If omitted, just cast the value to the proper data type.
        oconv: Callable[[Any], Any] | None = None
            Convert an output value to the right type for exporting.
            If omitted, export the value with no change.
        """
        self.key = key
        self.title = title
        self.dtype = dtype
        self.default = default
        self.iconv = iconv if iconv is not None else self.dtype
        self.oconv = oconv

    def __str__(self):
        return f"{type(self).__name__} {repr(self.title)} ({self.key})"


class OptionField(Field):
    """ Field based on a command line option. """

    def __init__(self, option: Option, **kwargs):
        """
        Parameters
        ----------
        option: click.Option
            Option from which to make the field.
        **kwargs
            Additional keyword arguments passed to `Field`.
        """
        super().__init__(option.name,
                         option.help,
                         option.type,
                         option.default,
                         **kwargs)


# Field calculation functions

# Note that each function takes a single argument: a Report instance.
# So these functions could be implemented as Report instance methods.
# But this implementation could cause confusion because no one Report
# class can work with all these methods.

def calc_taken(report: Report):
    """ Calculate the time taken in minutes. """
    delta = report.get_field(TimeEndedF) - report.get_field(TimeBeganF)
    minutes = (delta.seconds + 1e-6 * delta.microseconds) / 60.
    if minutes < 0.:
        raise ValueError(f"Time taken must be positive, but got {minutes} min")
    return minutes


# Field definitions

DATETIME_FORMAT = "%Y-%m-%d at %H:%M:%S"
DECIMAL_PRECISION = 3  # general precision for decimals
PERC_VEC_PRECISION = 1
TIME_TAKEN_PRECISION = 2


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


# General fields

VersionF = Field("version", "Version of SEISMIC-RNA", str, __version__)
BranchesF = Field("branches", "Branches", list, list())
SampleF = Field("sample", "Name of Sample", str)
RefF = Field("ref", "Name of Reference", str)
SectF = Field("sect", "Name of Section", str)
End5F = Field("end5", "5' end of Section", int)
End3F = Field("end3", "3' end of Section", int)
MinReadsF = OptionField(opt_min_reads)
TimeBeganF = Field("began", "Time Began", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime)
TimeEndedF = Field("ended", "Time Ended", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime)
TimeTakenF = Field("taken", "Time Taken (minutes)", float, calc_taken,
                   oconv=get_oconv_float(TIME_TAKEN_PRECISION))

# Align fields

IsDemultF = Field("demultiplexed", "Use demultiplexed mode", bool)
IsPairedEndF = Field("paired_end", "Use paired-end mode", bool)
PhredEncF = OptionField(opt_phred_enc)
UseFastqcF = OptionField(opt_fastqc)
UseCutadaptF = OptionField(opt_cutadapt)
CutadaptQ1 = OptionField(opt_cut_q1)
CutadaptQ2 = OptionField(opt_cut_q2)
CutadaptG1 = OptionField(opt_cut_g1)
CutadaptA1 = OptionField(opt_cut_a1)
CutadaptG2 = OptionField(opt_cut_g2)
CutadaptA2 = OptionField(opt_cut_a2)
CutadaptOverlap = OptionField(opt_cut_o)
CutadaptErrors = OptionField(opt_cut_e)
CutadaptIndels = OptionField(opt_cut_indels)
CutadaptNextSeq = OptionField(opt_cut_nextseq)
CutadaptNoTrimmed = OptionField(opt_cut_discard_trimmed)
CutadaptNoUntrimmed = OptionField(opt_cut_discard_untrimmed)
CutadaptMinLength = OptionField(opt_cut_m)
Bowtie2Local = OptionField(opt_bt2_local)
Bowtie2Discord = OptionField(opt_bt2_discordant)
Bowtie2Dovetail = OptionField(opt_bt2_dovetail)
Bowtie2Contain = OptionField(opt_bt2_contain)
Bowtie2Mixed = OptionField(opt_bt2_mixed)
Bowtie2Un = OptionField(opt_bt2_un)
Bowtie2ScoreMin = Field("bt2_score_min",
                        "Minimum score for a valid alignment with Bowtie2",
                        str)
Bowtie2MinLength = OptionField(opt_bt2_i)
Bowtie2MaxLength = OptionField(opt_bt2_x)
Bowtie2GBar = OptionField(opt_bt2_gbar)
Bowtie2SeedLength = OptionField(opt_bt2_l)
Bowtie2SeedInterval = OptionField(opt_bt2_s)
Bowtie2ExtTries = OptionField(opt_bt2_d)
Bowtie2Reseed = OptionField(opt_bt2_r)
Bowtie2Dpad = OptionField(opt_bt2_dpad)
Bowtie2Orient = OptionField(opt_bt2_orient)
MinMapQualF = OptionField(opt_min_mapq)
CramOutF = OptionField(opt_cram)
ReadsInitF = Field("reads_init", "Number of reads initially", int)
ReadsTrimF = Field("reads_trim", "Number of reads after trimming", int)
ReadsAlignF = Field("reads_align", "Number of reads after alignment", dict,
                    iconv=iconv_dict_str_int)
ReadsDedupF = Field("reads_filter",
                    "Number of reads after filtering",
                    dict,
                    iconv=iconv_dict_str_int)
ReadsRefs = Field("reads_refs",
                  "Number of reads aligned by reference",
                  dict,
                  iconv=iconv_dict_str_int)

# Relate fields
NumReadsRelF = Field("n_reads_rel", "Number of Reads", int)
NumBatchF = Field("n_batches", "Number of Batches", int)
ChecksumsF = Field("checksums", "MD5 Checksums of Batches", dict)
RefseqChecksumF = Field("refseq_checksum",
                        "MD5 Checksum of Reference Sequence File",
                        str)
PooledSamplesF = Field("pooled_samples", "Pooled Samples", list)

# Mask fields
CountMutsF = Field("count_muts",
                   "Count the Following as Mutations",
                   HalfRelPattern,
                   iconv=HalfRelPattern.from_report_format,
                   oconv=HalfRelPattern.to_report_format)
CountRefsF = Field("count_refs",
                   "Count the Following as Matches",
                   HalfRelPattern,
                   iconv=HalfRelPattern.from_report_format,
                   oconv=HalfRelPattern.to_report_format)
ExclPolyAF = Field("exclude_polya",
                   "Exclude Poly(A) Sequences of at Least This Length (nt)",
                   int)
ExclGUF = Field("exclude_gu", "Exclude G/U Bases", bool)
ExclUserPosF = Field("exclude_pos",
                     "Exclude User-Defined Positions",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int)
MinInfoPosF = Field("min_ninfo_pos",
                    "Minimum Number of Informative Reads per Position",
                    int)
MinMutPosF = Field("min_fmut_pos",
                   "Minimum Fraction of Mutations per Position",
                   float,
                   oconv=get_oconv_float())
MaxMutPosF = Field("max_fmut_pos",
                   "Maximum Fraction of Mutations per Position",
                   float,
                   oconv=get_oconv_float())
MinMutGapF = Field("min_mut_gap",
                   "Minimum Gap Between Mutations (nt)",
                   int)
MinInfoReadF = Field("min_finfo_read",
                     "Minimum Fraction of Informative Positions per Read",
                     float,
                     oconv=get_oconv_float())
MaxMutReadF = Field("max_fmut_read",
                    "Maximum Fraction of Mutations per Read",
                    float,
                    oconv=get_oconv_float())
PosInitF = Field("pos_init",
                 "Positions Initially Given",
                 np.ndarray,
                 iconv=iconv_array_int,
                 oconv=oconv_array_int)
PosCutPolyAF = Field("pos_polya",
                     "Positions Cut -- Poly(A) Sequence",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int)
PosCutGUF = Field("pos_gu",
                  "Positions Cut -- G/U Base",
                  np.ndarray,
                  iconv=iconv_array_int,
                  oconv=oconv_array_int)
PosCutUserF = Field("pos_user",
                    "Positions Cut -- User-Specified",
                    np.ndarray,
                    iconv=iconv_array_int,
                    oconv=oconv_array_int)
PosCutLoInfoF = Field("pos_min_ninfo",
                      "Positions Cut -- Too Few Informative Reads",
                      np.ndarray,
                      iconv=iconv_array_int,
                      oconv=oconv_array_int)
PosCutLoMutF = Field("pos_min_fmut",
                     "Positions Cut -- Too Few Mutations",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int)
PosCutHiMutF = Field("pos_max_fmut",
                     "Positions Cut -- Too Many Mutations",
                     np.ndarray,
                     iconv=iconv_array_int,
                     oconv=oconv_array_int)
PosKeptF = Field("pos_kept",
                 "Positions Ultimately Kept",
                 np.ndarray,
                 iconv=iconv_array_int,
                 oconv=oconv_array_int)
NumPosInitF = Field("n_pos_init",
                    "Number of Positions Initially Given",
                    int)
NumPosCutPolyAF = Field("n_pos_polya",
                        "Number of Positions Cut -- Poly(A) Sequence",
                        int)
NumPosCutGUF = Field("n_pos_gu",
                     "Number of Positions Cut -- G/U Base",
                     int)
NumPosCutUserF = Field("n_pos_user",
                       "Number of Positions Cut -- User-Specified",
                       int)
NumPosCutLoInfoF = Field("n_pos_min_ninfo",
                         "Number of Positions Cut -- Too Few Informative Reads",
                         int)
NumPosCutLoMutF = Field("n_pos_min_fmut",
                        "Number of Positions Cut -- Too Few Mutations",
                        int)
NumPosCutHiMutF = Field("n_pos_max_fmut",
                        "Number of Positions Cut -- Too Many Mutations",
                        int)
NumPosKeptF = Field("n_pos_kept",
                    "Number of Positions Ultimately Kept",
                    int)
NumReadsInitF = Field("n_reads_init",
                      "Number of Reads Initially Given",
                      int)
NumReadsLoInfoF = Field("n_reads_min_finfo",
                        "Number of Reads Cut -- Too Few Informative Positions",
                        int)
NumReadsHiMutF = Field("n_reads_max_fmut",
                       "Number of Reads Cut -- Too Many Mutations",
                       int)
NumReadsCloseMutF = Field("n_reads_min_gap",
                          "Number of Reads Cut -- Mutations Too Close Together",
                          int)
NumReadsKeptF = Field("n_reads_kept",
                      "Number of Reads Ultimately Kept",
                      int)
NumUniqReadKeptF = Field("n_uniq_reads",
                         "Number of Unique Bit Vectors",
                         int)

# Cluster fields

MinIterClustF = Field("min_iter",
                      "Minimum EM Iterations per Cluster",
                      int)
MaxIterClustF = Field("max_iter",
                      "Maximum EM Iterations per Cluster",
                      int)
ClustConvThreshF = Field("conv_thresh",
                         "Convergence Threshold for Log Likelihood",
                         float,
                         oconv=get_oconv_float())
MaxClustsF = Field("max_order", "Maximum Number of Clusters", int)
ClustNumRunsF = Field("num_runs", "Number of Independent EM Runs", int)
NumClustsF = Field("best_order", "Optimal Number of Clusters", int)
ClustsBicF = Field("bic", "Bayesian Information Criterion per Order", dict,
                   iconv=iconv_int_keys,
                   oconv=get_oconv_dict_float())
ClustsConvF = Field("converged", "Iterations Until Convergence per Run", dict,
                    iconv=iconv_int_keys)
ClustsLogLikesF = Field("log_likes", "Log Likelihood per Run", dict,
                        iconv=iconv_int_keys,
                        oconv=get_oconv_dict_list_float())
ClustsLikeMeanF = Field("log_like_mean", "Mean Log Likelihood per Order", dict,
                        iconv=iconv_int_keys,
                        oconv=get_oconv_dict_float())
ClustsLikeStdF = Field("log_like_std",
                       "Std. Dev. Log Likelihood per Order",
                       dict,
                       iconv=iconv_int_keys,
                       oconv=get_oconv_dict_float())
ClustsVarInfoF = Field("var_info",
                       "Variation of Information per Order",
                       dict,
                       iconv=iconv_int_keys,
                       oconv=get_oconv_dict_float())

# Fold fields

ProfileF = Field("profile", "Name of Profile", str)
Quantile = OptionField(opt_quantile)
FoldTempF = OptionField(opt_fold_temp)
FoldMaxDistF = OptionField(opt_fold_md)
FoldMinFreeEnergyF = OptionField(opt_fold_mfe)
FoldMaxStructsF = OptionField(opt_fold_max)
FoldPercent = OptionField(opt_fold_percent)


# Field managing functions

@cache
def fields():
    return [member for _, member in getmembers(sys.modules[__name__])
            if isinstance(member, Field)]


@cache
def field_keys() -> dict[str, Field]:
    return {field.key: field for field in fields()}


@cache
def field_titles() -> dict[str, Field]:
    return {field.title: field for field in fields() if field.title}


def lookup_key(key: str):
    """ Get a field by its key. """
    try:
        return field_keys()[key]
    except KeyError:
        raise ValueError(f"Invalid report key: {repr(key)}")


def lookup_title(title: str):
    """ Get a field by its title. """
    if not title:
        raise ValueError("Got blank title for field")
    try:
        return field_titles()[title]
    except KeyError:
        raise ValueError(f"Invalid report field: {repr(title)}")


def key_to_title(key: str):
    """ Map a field's key to its title. """
    return lookup_key(key).title


def default_key(key: str):
    """ Get the default value of a field by its key. """
    if (default := lookup_key(key).default) is None:
        raise ValueError(f"Field {repr(key_to_title(key))} has no default")
    return default


# Report classes

class Report(FileIO, ABC):
    """ Abstract base class for a report from a step. """

    @classmethod
    @abstractmethod
    def fields(cls):
        """ All fields of the report. """
        return [BranchesF, TimeBeganF, TimeEndedF, TimeTakenF, VersionF]

    @classmethod
    @cache
    def field_keys(cls):
        """ Keys of all fields of the report. """
        return [field.key for field in cls.fields()]

    @classmethod
    def from_dict(cls, odata: dict[str, Any]):
        """ Convert a dict of raw values (keyed by the titles of their
        fields) into a dict of encoded values (keyed by the keys of
        their fields), from which a new Report is instantiated. """
        if not isinstance(odata, dict):
            raise TypeError(f"Expected dict, but got {type(odata).__name__}")
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
    def load(cls, file: Path) -> Report:
        with open(file) as f:
            report = cls.from_dict(json.load(f))
        # Ensure that the path-related fields in the JSON data match the
        # actual path of the JSON file.
        top, path_fields = cls.parse_path(file)
        for key, value in report.path_field_values().items():
            if value != path_fields.get(key):
                raise ValueError(f"Got different values for field {repr(key)} "
                                 f"in path ({repr(path_fields.get(key))}) and "
                                 f"contents ({repr(value)}) of report {file}")
        return report

    @classmethod
    def _auto_default_fields(cls):
        return [BranchesF, TimeTakenF, VersionF]

    @classmethod
    def _auto_init_kwargs(cls, **kwargs):
        """ Automatic keyword arguments for __init__. """
        return {field.key: field.default
                for field in cls._auto_default_fields()} | kwargs

    def __init__(self, **kwargs: Any | Callable[[Report], Any]):
        kwargs = self._auto_init_kwargs(**kwargs)
        defaulted = dict()
        for key in self.field_keys():
            # Try to get the value of the field from the report.
            try:
                value = kwargs.pop(key)
            except KeyError:
                # If the report file is missing that field (e.g. because
                # it came from a different version of SEISMIC-RNA), then
                # for cross-version compatibility, use the default value
                # of the field.
                value = default_key(key)
                defaulted[key_to_title(key)] = value
            if callable(value):
                # If the value of the keyword argument is callable, then
                # it must accept one argument -- self -- and return the
                # value of the attribute.
                value = value(self)
            setattr(self, key, value)
        if kwargs:
            # If the report file has extra fields (e.g. because it came
            # from a different version of SEISMIC-RNA), then just log a
            # warning and ignore the extra fields (to make different
            # versions compatible).
            logger.warning(f"Got extra fields for {type(self).__name__}: "
                           f"{list(kwargs)}")
        if defaulted:
            # If the report file was missing keyword arguments that have
            # default values, AND if parsing the report file succeeded,
            # then warn about the default values.
            logger.warning(f"Missing fields for {type(self).__name__} "
                           f"and using defaults: {defaulted}")

    def get_field(self, field: Field, missing_ok: bool = False):
        """ Return the value of a field of the report using the field
        instance directly, not its key. """
        try:
            return getattr(self, field.key)
        except AttributeError:
            if missing_ok:
                return None
            raise

    def to_dict(self):
        """ Return a dict of raw values of the fields, keyed by the
        titles of their fields. """
        odata = dict()
        for key in self.field_keys():
            field = lookup_key(key)
            # Output only the fields with non-blank titles.
            value = self.get_field(field)
            if field.oconv is not None:
                # Convert the value to the proper output value.
                value = field.oconv(value)
            odata[field.title] = value
        return odata

    def save(self, top: Path, force: bool = False):
        """ Save the report to a JSON file. """
        text = json.dumps(self.to_dict(), indent=4)
        save_path = self.get_path(top)
        if need_write(save_path, force):
            with open(save_path, write_mode(force)) as f:
                f.write(text)
            logger.info(f"Wrote {self} to {save_path}")
        return save_path

    def __setattr__(self, key: str, value: Any):
        """ Validate the attribute name and value before setting it. """
        if key not in self.field_keys():
            raise ValueError(
                f"Invalid field for {type(self).__name__}: {repr(key)}"
            )
        super().__setattr__(key, value)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.to_dict() == other.to_dict()


class RefseqReport(Report, RefIO, ABC):
    """ Report associated with a reference sequence file. """

    @classmethod
    @abstractmethod
    def fields(cls):
        return [RefseqChecksumF] + super().fields()


class BatchedReport(Report, ABC):
    """ Report with a number of data batches (one file per batch). """

    @classmethod
    @abstractmethod
    def fields(cls):
        return [NumBatchF, ChecksumsF] + super().fields()

    @classmethod
    @abstractmethod
    def _batch_types(cls) -> tuple[type[ReadBatchIO], ...]:
        """ Type(s) of batch(es) for the report. """

    @classmethod
    @cache
    def batch_types(cls) -> dict[str, type[ReadBatchIO]]:
        """ Type(s) of batch(es) for the report, keyed by name. """
        return {batch_type.btype(): batch_type
                for batch_type in cls._batch_types()}

    @classmethod
    def get_batch_type(cls, btype: str | None = None) -> type[ReadBatchIO]:
        """ Return a valid type of batch based on its name. """
        if btype is None:
            batch_types = list(cls.batch_types().values())
            if (ntypes := len(batch_types)) != 1:
                raise ValueError(f"btype is optional only if there is exactly "
                                 f"one type of batch, but got {ntypes} types")
            return batch_types[0]
        return cls.batch_types()[btype]


class BatchedRefseqReport(BatchedReport, RefseqReport, ABC):
    """ Convenience class used as a base for several Report classes. """

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
