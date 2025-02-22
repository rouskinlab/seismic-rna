from __future__ import annotations

import json
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from functools import cache
from inspect import getmembers
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Hashable, Iterable

import numpy as np
from click import Option

from .arg import (opt_phred_enc,
                  opt_fastp,
                  opt_fastp_5,
                  opt_fastp_3,
                  opt_fastp_w,
                  opt_fastp_m,
                  opt_fastp_poly_g_min_len,
                  opt_fastp_poly_x,
                  opt_fastp_poly_x_min_len,
                  opt_fastp_adapter_trimming,
                  opt_fastp_adapter_1,
                  opt_fastp_adapter_2,
                  opt_fastp_adapter_fasta,
                  opt_fastp_detect_adapter_for_pe,
                  opt_fastp_poly_g,
                  opt_fastp_min_length,
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
                  opt_sep_strands,
                  opt_f1r2_fwd,
                  opt_rev_label,
                  opt_min_reads,
                  opt_min_mapq,
                  opt_insert3,
                  opt_ambindel,
                  opt_overhangs,
                  opt_clip_end5,
                  opt_clip_end3,
                  opt_fold_temp,
                  opt_fold_md,
                  opt_fold_mfe,
                  opt_fold_max,
                  opt_fold_percent,
                  opt_quantile,
                  opt_quick_unbias,
                  opt_quick_unbias_thresh,
                  opt_min_ncov_read,
                  opt_min_finfo_read,
                  opt_min_mut_gap,
                  opt_min_ninfo_pos,
                  opt_max_fmut_pos,
                  opt_max_mask_iter,
                  opt_em_runs,
                  opt_em_thresh,
                  opt_min_em_iter,
                  opt_max_em_iter,
                  opt_max_fmut_read,
                  opt_min_clusters,
                  opt_max_clusters,
                  opt_jackpot,
                  opt_jackpot_conf_level,
                  opt_max_jackpot_quotient,
                  opt_max_pearson_run,
                  opt_min_marcd_run,
                  opt_max_loglike_vs_best,
                  opt_min_pearson_vs_best,
                  opt_max_marcd_vs_best,
                  opt_try_all_ks,
                  opt_write_all_ks,
                  opt_mask_gu,
                  opt_mask_polya,
                  opt_mask_discontig,
                  opt_min_phred)
from .error import InconsistentValueError
from .io import SampleFileIO, ReadBatchIO, RefFileIO, RegFileIO
from .logs import logger
from .path import flatten_branches
from .rel import HalfRelPattern
from .version import __version__
from .write import need_write, write_mode


# Field class

class ReportField(object):
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


class OptionReportField(ReportField):
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


def calc_dt_minutes(began: datetime, ended: datetime):
    """ Calculate the time taken in minutes. """
    delta = ended - began
    minutes = (delta.seconds + 1e-6 * delta.microseconds) / 60.
    if minutes < 0.:
        raise ValueError(f"Time taken must be positive, but got {minutes} min")
    return minutes


def calc_taken(report: Report):
    """ Calculate the time taken in minutes. """
    return calc_dt_minutes(report.get_field(TimeBeganF),
                           report.get_field(TimeEndedF))


# Field definitions

DATETIME_FORMAT = "%Y-%m-%d at %H:%M:%S"
DECIMAL_PRECISION = 3  # general precision for decimals
PERC_VEC_PRECISION = 1
TIME_TAKEN_PRECISION = 2


def iconv_int_keys(mapping: dict[Any, Any]):
    return {int(key): value for key, value in mapping.items()}


def iconv_array_int(nums: list[int]):
    return np.asarray(nums, dtype=int)


def iconv_dict_str_int(mapping: dict[Any, Any]) -> dict[str, int]:
    return {str(key): int(value) for key, value in mapping.items()}


def iconv_dict_str_dict_int_dict_int_int(
        mapping: dict[Any, dict[Any, dict[Any, Any]]]
) -> dict[str, dict[int, dict[int, int]]]:
    return {str(key1): {int(key2): {int(key3): int(val3)
                                    for key3, val3 in val2.items()}
                        for key2, val2 in val1.items()}
            for key1, val1 in mapping.items()}


@cache
def get_oconv_float(precision: int = DECIMAL_PRECISION):
    def oconv_float(num):
        return float(round(num, precision))

    return oconv_float


@cache
def get_oconv_list(dtype: type, precision: int = DECIMAL_PRECISION):
    if dtype is float:
        oconv_func = get_oconv_float(precision)
    else:
        oconv_func = dtype

    def oconv_list(nums: Iterable):
        return list(map(oconv_func, nums))

    return oconv_list


@cache
def get_oconv_dict(dtype: type, precision: int = DECIMAL_PRECISION):
    if dtype is float:
        oconv_func = get_oconv_float(precision)
    else:
        oconv_func = dtype

    def oconv_dict(dnum: dict):
        return {d: oconv_func(num) for d, num in dnum.items()}

    return oconv_dict


@cache
def get_oconv_dict_list(dtype: type, precision: int = DECIMAL_PRECISION):
    oconv_list = get_oconv_list(dtype, precision=precision)

    def oconv_dict_list(dnums: dict[Hashable, Iterable]):
        return {d: oconv_list(nums) for d, nums in dnums.items()}

    return oconv_dict_list


def iconv_datetime(text: str):
    return datetime.strptime(text, DATETIME_FORMAT)


def oconv_datetime(dtime: datetime):
    return dtime.strftime(DATETIME_FORMAT)


# General fields
VersionF = ReportField("version", "Version of SEISMIC-RNA", str, __version__)
BranchesF = ReportField("branches", "Branches", dict)
SampleF = ReportField("sample", "Sample", str)
RefF = ReportField("ref", "Reference", str)
RegF = ReportField("reg", "Region", str)
End5F = ReportField("end5", "Region 5' end", int)
End3F = ReportField("end3", "Region 3' end", int)
MinReadsF = OptionReportField(opt_min_reads)
TimeBeganF = ReportField("began",
                         "Time began",
                         datetime,
                         iconv=iconv_datetime,
                         oconv=oconv_datetime)
TimeEndedF = ReportField("ended",
                         "Time ended",
                         datetime,
                         iconv=iconv_datetime,
                         oconv=oconv_datetime)
TimeTakenF = ReportField("taken",
                         "Time taken (minutes)",
                         float,
                         calc_taken,
                         oconv=get_oconv_float(TIME_TAKEN_PRECISION))

# Align fields
IsDemultF = ReportField("demultiplexed", "Use demultiplexed mode", bool)
IsPairedEndF = ReportField("paired_end", "Use paired-end mode", bool)
PhredEncF = OptionReportField(opt_phred_enc)
UseFastpF = OptionReportField(opt_fastp)
Fastp5F = OptionReportField(opt_fastp_5)
Fastp3F = OptionReportField(opt_fastp_3)
FastpWF = OptionReportField(opt_fastp_w)
FastpMF = OptionReportField(opt_fastp_m)
FastpPolyGF = OptionReportField(opt_fastp_poly_g)
FastpPolyGMinLenF = OptionReportField(opt_fastp_poly_g_min_len)
FastpPolyXF = OptionReportField(opt_fastp_poly_x)
FastpPolyXMinLenF = OptionReportField(opt_fastp_poly_x_min_len)
FastpAdapterTrimmingF = OptionReportField(opt_fastp_adapter_trimming)
FastpAdapter1F = OptionReportField(opt_fastp_adapter_1)
FastpAdapter2F = OptionReportField(opt_fastp_adapter_2)
FastpAdapterFastaF = OptionReportField(opt_fastp_adapter_fasta)
FastpDetectAdapterForPEF = OptionReportField(opt_fastp_detect_adapter_for_pe)
FastpMinLengthF = OptionReportField(opt_fastp_min_length)
Bowtie2Local = OptionReportField(opt_bt2_local)
Bowtie2Discord = OptionReportField(opt_bt2_discordant)
Bowtie2Dovetail = OptionReportField(opt_bt2_dovetail)
Bowtie2Contain = OptionReportField(opt_bt2_contain)
Bowtie2Mixed = OptionReportField(opt_bt2_mixed)
Bowtie2Un = OptionReportField(opt_bt2_un)
Bowtie2ScoreMin = ReportField("bt2_score_min",
                              "Discard alignments that score below this threshold",
                              str)
Bowtie2MinLengthF = OptionReportField(opt_bt2_i)
Bowtie2MaxLengthF = OptionReportField(opt_bt2_x)
Bowtie2GBarF = OptionReportField(opt_bt2_gbar)
Bowtie2SeedLength = OptionReportField(opt_bt2_l)
Bowtie2SeedInterval = OptionReportField(opt_bt2_s)
Bowtie2ExtTries = OptionReportField(opt_bt2_d)
Bowtie2Reseed = OptionReportField(opt_bt2_r)
Bowtie2Dpad = OptionReportField(opt_bt2_dpad)
Bowtie2Orient = OptionReportField(opt_bt2_orient)
MinMapQualF = OptionReportField(opt_min_mapq)
SepStrandsF = OptionReportField(opt_sep_strands)
F1R2FwdF = OptionReportField(opt_f1r2_fwd)
RevLabelF = OptionReportField(opt_rev_label)
AlignReadsInitF = ReportField("align_reads_init", "Number of reads in the FASTQ file(s)", int)
ReadsTrimF = ReportField("reads_trim", "Number of reads after trimming", int)
ReadsAlignF = ReportField("reads_align",
                          "Number of reads after alignment",
                          dict,
                          iconv=iconv_dict_str_int)
ReadsDedupF = ReportField("reads_filter",
                          "Number of reads after filtering",
                          dict,
                          iconv=iconv_dict_str_int)
ReadsRefsF = ReportField("reads_refs",
                         "Number of reads aligned to each reference",
                         dict,
                         iconv=iconv_dict_str_int)
RefFastaChecksumF = ReportField("ref_fasta_checksum",
                                "Checksum of the reference fasta (SHA-512)",
                                str)
FastqChecksumsF = ReportField("fastq_checksums",
                              "Checksum(s) of the input fastq(s) (SHA-512)",
                              dict)
XamChecksumF = ReportField("xam_checksum",
                           "Checksum of the input xam (SHA-512)",
                           str)

# Relate fields
NumReadsXamF = ReportField("n_reads_xam",
                           "Number of reads in SAM/BAM/CRAM file",
                           int)
NumReadsRelF = ReportField("n_reads_rel",
                           "Number of reads processed by relate",
                           int)
NumBatchesF = ReportField("n_batches", "Number of batches", int)
ChecksumsF = ReportField("checksums", "Checksums of batches (SHA-512)", dict)
RefseqChecksumF = ReportField("refseq_checksum",
                              "Checksum of reference sequence (SHA-512)",
                              str)
Insert3F = OptionReportField(opt_insert3)
AmbindelF = OptionReportField(opt_ambindel)
OverhangsF = OptionReportField(opt_overhangs)
MinPhredF = OptionReportField(opt_min_phred)
ClipEnd5F = OptionReportField(opt_clip_end5)
ClipEnd3F = OptionReportField(opt_clip_end3)

# Pool fields
PooledSamplesF = ReportField("pooled_samples", "Pooled samples", list)

# Mask fields
mask_iter_no_convergence = 0
CountMutsF = ReportField("count_muts",
                         "Count as mutations",
                         HalfRelPattern,
                         iconv=HalfRelPattern.from_report_format,
                         oconv=HalfRelPattern.to_report_format)
CountRefsF = ReportField("count_refs",
                         "Count as matches",
                         HalfRelPattern,
                         iconv=HalfRelPattern.from_report_format,
                         oconv=HalfRelPattern.to_report_format)
ExclPolyAF = OptionReportField(opt_mask_polya)
ExclGUF = OptionReportField(opt_mask_gu)
ExclListPosF = ReportField("mask_pos",
                           "Mask additional positions from a list",
                           np.ndarray,
                           iconv=iconv_array_int,
                           oconv=get_oconv_list(int))
MinNInfoPosF = OptionReportField(opt_min_ninfo_pos)
MaxFMutPosF = OptionReportField(opt_max_fmut_pos)
MinNCovReadF = OptionReportField(opt_min_ncov_read)
DiscontigF = OptionReportField(opt_mask_discontig)
MinMutGapF = OptionReportField(opt_min_mut_gap)
QuickUnbiasF = OptionReportField(opt_quick_unbias)
QuickUnbiasThreshF = OptionReportField(opt_quick_unbias_thresh)
MinFInfoReadF = OptionReportField(opt_min_finfo_read)
MaxFMutReadF = OptionReportField(opt_max_fmut_read)
MaxMaskIterF = OptionReportField(opt_max_mask_iter)
PosCutPolyAF = ReportField("pos_polya",
                           "Positions in stretches of consecutive A bases",
                           np.ndarray,
                           iconv=iconv_array_int,
                           oconv=get_oconv_list(int))
PosCutGUF = ReportField("pos_gu",
                        "Positions with G or U bases",
                        np.ndarray,
                        iconv=iconv_array_int,
                        oconv=get_oconv_list(int))
PosCutListF = ReportField("pos_list",
                          "Positions masked from a list",
                          np.ndarray,
                          iconv=iconv_array_int,
                          oconv=get_oconv_list(int))
PosCutLoInfoF = ReportField("pos_min_ninfo",
                            "Positions with too few informative base calls",
                            np.ndarray,
                            iconv=iconv_array_int,
                            oconv=get_oconv_list(int))
PosCutHiMutF = ReportField("pos_max_fmut",
                           "Positions with too many mutations",
                           np.ndarray,
                           iconv=iconv_array_int,
                           oconv=get_oconv_list(int))
PosKeptF = ReportField("pos_kept",
                       "Positions kept after masking",
                       np.ndarray,
                       iconv=iconv_array_int,
                       oconv=get_oconv_list(int))
NumPosInitF = ReportField("n_pos_init",
                          "Total number of positions in the region",
                          int)
NumPosCutPolyAF = ReportField("n_pos_polya",
                              "Number of positions in stretches of consecutive "
                              "A bases",
                              int)
NumPosCutGUF = ReportField("n_pos_gu",
                           "Number of positions with G or U bases",
                           int)
NumPosCutListF = ReportField("n_pos_list",
                             "Number of positions masked from a list",
                             int)
NumPosCutLoInfoF = ReportField("n_pos_min_ninfo",
                               "Number of positions with too few informative "
                               "base calls",
                               int)
NumPosCutHiMutF = ReportField("n_pos_max_fmut",
                              "Number of positions with too many mutations",
                              int)
NumPosKeptF = ReportField("n_pos_kept",
                          "Number of positions kept after masking",
                          int)
NumReadsInitF = ReportField("n_reads_init",
                            "Total number of reads before masking",
                            int)
NumReadCutListF = ReportField("n_reads_list",
                              "Number of reads masked from a list",
                              int)
NumReadsLoNCovF = ReportField(
    "n_reads_min_ncov",
    "Number of reads with too few bases covering the region",
    int
)
NumDiscontigF = ReportField("n_reads_discontig",
                            "Number of reads with discontiguous mates",
                            int)
NumReadsLoInfoF = ReportField("n_reads_min_finfo",
                              "Number of reads with too few informative "
                              "base calls",
                              int)
NumReadsHiMutF = ReportField("n_reads_max_fmut",
                             "Number of reads with too many mutations",
                             int)
NumReadsCloseMutF = ReportField("n_reads_min_gap",
                                "Number of reads with two mutations too close",
                                int)
NumReadsKeptF = ReportField("n_reads_kept",
                            "Number of reads kept after masking",
                            int)
NumMaskIterF = ReportField("n_mask_iter",
                           "Number of iterations until convergence "
                           f"({mask_iter_no_convergence} if not converged)",
                           int)

# Cluster fields

NumUniqReadKeptF = ReportField("n_uniq_reads",
                               "Number of unique reads",
                               int)
MinIterClustF = OptionReportField(opt_min_em_iter)
MaxIterClustF = OptionReportField(opt_max_em_iter)
ClustConvThreshF = OptionReportField(opt_em_thresh)
MinClustsF = OptionReportField(opt_min_clusters)
MaxClustsF = OptionReportField(opt_max_clusters)
JackpotF = OptionReportField(opt_jackpot)
JackpotConfLevelF = OptionReportField(opt_jackpot_conf_level)
MaxJackpotQuotientF = OptionReportField(opt_max_jackpot_quotient)
MaxPearsonRunF = OptionReportField(opt_max_pearson_run)
MinMARCDRunF = OptionReportField(opt_min_marcd_run)
MaxLogLikeVsBestF = OptionReportField(opt_max_loglike_vs_best)
MinPearsonVsBestF = OptionReportField(opt_min_pearson_vs_best)
MaxMARCDVsBestF = OptionReportField(opt_max_marcd_vs_best)
TryAllKsF = OptionReportField(opt_try_all_ks)
WriteAllKsF = OptionReportField(opt_write_all_ks)
ClustNumRunsF = OptionReportField(opt_em_runs)
EMKPassingF = ReportField("em_k_passing",
                          "Whether each number of clusters (K) passed filters",
                          dict,
                          iconv=iconv_int_keys,
                          oconv=get_oconv_dict(bool))
KsWrittenF = ReportField("ks_written",
                         "Numbers of clusters written to batches",
                         list)
BestKF = ReportField("best_k", "Best number of clusters", int)

# Join fields

JoinedRegionsF = ReportField("joined_regions", "Joined regions", list)
JoinedClustersF = ReportField("joined_clusters",
                              "Joined clusters",
                              dict,
                              iconv=iconv_dict_str_dict_int_dict_int_int)

# Fold fields

ProfileF = ReportField("profile", "Profile", str)
Quantile = OptionReportField(opt_quantile)
FoldTempF = OptionReportField(opt_fold_temp)
FoldMaxDistF = OptionReportField(opt_fold_md)
FoldMinFreeEnergyF = OptionReportField(opt_fold_mfe)
FoldMaxStructsF = OptionReportField(opt_fold_max)
FoldPercent = OptionReportField(opt_fold_percent)


# Field exceptions


class ReportFieldError(RuntimeError):
    """ Any error involving a field of a report. """


class ReportFieldTypeError(ReportFieldError, TypeError):
    pass


class ReportFieldValueError(ReportFieldError, ValueError):
    pass


class ReportFieldKeyError(ReportFieldError, KeyError):
    pass


class ReportFieldAttributeError(ReportFieldError, AttributeError):
    pass


class InvalidReportFieldKeyError(ReportFieldKeyError):
    """ The key does not belog to an actual report field. """


class InvalidReportFieldTitleError(ReportFieldKeyError):
    """ The title does not belog to an actual report field. """


class MissingFieldWithNoDefaultError(ReportFieldValueError):
    """ The default value is requested of a field with no default. """


class ReportDoesNotHaveFieldError(ReportFieldAttributeError):
    """ A report does not contain this type of field. """


# Field managing functions

@cache
def fields():
    return [member for _, member in getmembers(sys.modules[__name__])
            if isinstance(member, ReportField)]


@cache
def field_keys() -> dict[str, ReportField]:
    keys = dict()
    for field in fields():
        if field.key:
            assert field.key not in keys
            keys[field.key] = field
    return keys


@cache
def field_titles() -> dict[str, ReportField]:
    titles = dict()
    for field in fields():
        if field.title:
            assert field.title not in titles
            titles[field.title] = field
    return titles


def lookup_key(key: str):
    """ Get a field by its key. """
    try:
        return field_keys()[key]
    except KeyError:
        raise InvalidReportFieldKeyError(key) from None


def lookup_title(title: str):
    """ Get a field by its title. """
    if not title:
        raise ValueError("Got blank title for field")
    try:
        return field_titles()[title]
    except KeyError:
        raise InvalidReportFieldTitleError(title) from None


def key_to_title(key: str):
    """ Map a field's key to its title. """
    return lookup_key(key).title


def default_key(key: str):
    """ Get the default value of a field by its key. """
    if (default := lookup_key(key).default) is None:
        raise MissingFieldWithNoDefaultError(key_to_title(key))
    return default


# Report classes

class Report(SampleFileIO, ABC):
    """ Abstract base class for a report from a step. """

    @classmethod
    def get_ident_report_fields(cls) -> list[ReportField]:
        """ Identification fields of the report. """
        return [SampleF, BranchesF]

    @classmethod
    def get_param_report_fields(cls) -> list[ReportField]:
        """ Parameter fields of the report. """
        return []

    @classmethod
    def get_result_report_fields(cls) -> list[ReportField]:
        """ Result fields of the report. """
        return []

    @classmethod
    def get_checksum_report_fields(cls) -> list[ReportField]:
        """ Checksum fields of the report. """
        return []

    @classmethod
    def get_meta_report_fields(cls) -> list[ReportField]:
        """ Metadata fields of the report. """
        return [TimeBeganF,
                TimeEndedF,
                TimeTakenF,
                VersionF]

    @classmethod
    @cache
    def get_report_fields(cls):
        """ All fields of the report. """
        return list(chain(cls.get_ident_report_fields(),
                          cls.get_param_report_fields(),
                          cls.get_result_report_fields(),
                          cls.get_checksum_report_fields(),
                          cls.get_meta_report_fields()))

    @classmethod
    @cache
    def get_field_keys(cls):
        """ Keys of all fields of the report. """
        return [field.key for field in cls.get_report_fields()]

    @classmethod
    @cache
    def get_field_keys_set(cls):
        """ Same as get_field_keys but caches and returns a set for fast
        membership checking. """
        return set(cls.get_field_keys())

    @classmethod
    def from_dict(cls, odata: dict[str, Any]):
        """ Convert a dict of raw values (keyed by the titles of their
        fields) into a dict of encoded values (keyed by the keys of
        their fields), from which a new Report is instantiated. """
        if not isinstance(odata, dict):
            raise TypeError(odata)
        # Read every raw value, keyed by the title of its field.
        idata = dict()
        for title, value in odata.items():
            # Get the field corresponding to the title.
            try:
                field = lookup_title(title)
            except InvalidReportFieldTitleError as error:
                logger.warning(error)
            else:
                # Cast the value to the input type; key it by the field.
                idata[field.key] = field.iconv(value)
        # Instantiate and return a new Report from the values.
        return cls(**idata)

    @classmethod
    def load(cls, file: str | Path) -> Report:
        logger.routine(f"Began loading {cls.__name__} from {file}")
        with open(file) as f:
            report = cls.from_dict(json.load(f))
        # Ensure that the path-related fields in the JSON data match the
        # actual path of the JSON file.
        top, path_fields = cls.parse_path(file)
        for key, value in report.get_path_field_values().items():
            if key == BranchesF.key:
                # The branches in the report must be flattened to match
                # the value from parsin the path.
                value = flatten_branches(value)
            if value != path_fields.get(key):
                raise InconsistentValueError(
                    f"Got different values for {repr(key)} in path "
                    f"({repr(path_fields.get(key))}) and contents "
                    f"({repr(value)}) of {file}"
                )
        logger.routine(f"Ended loading {cls.__name__} from {file}")
        return report

    @classmethod
    def _get_auto_default_fields(cls):
        return [TimeTakenF, VersionF]

    @classmethod
    def _get_auto_init_kwargs(cls):
        """ Automatic keyword arguments for __init__. """
        return {field.key: field.default
                for field in cls._get_auto_default_fields()}

    def __init__(self, **kwargs: Any | Callable[[Report], Any]):
        kwargs = self._get_auto_init_kwargs() | kwargs
        defaulted = dict()
        for key in self.get_field_keys():
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
            logger.warning(f"{type(self).__name__} got extra fields, "
                           f"which it discarded: {list(kwargs)}")
        if defaulted:
            # If the report file was missing keyword arguments that have
            # default values, AND if parsing the report file succeeded,
            # then warn about the default values.
            logger.warning(f"{type(self).__name__} is missing fields "
                           f"and using defaults: {defaulted}")

    def get_field(self, field: ReportField, missing_ok: bool = False):
        """ Return the value of a field of the report using the field
        instance directly, not its key. """
        try:
            return getattr(self, field.key)
        except AttributeError:
            if missing_ok:
                return None
            raise ReportDoesNotHaveFieldError(
                f"{type(self).__name__}.{field.key}"
            ) from None

    def to_dict(self):
        """ Return a dict of raw values of the fields, keyed by the
        titles of their fields. """
        odata = dict()
        for key in self.get_field_keys():
            field = lookup_key(key)
            # Output only the fields with non-blank titles.
            value = self.get_field(field)
            if field.oconv is not None:
                # Convert the value to the proper output value.
                value = field.oconv(value)
            odata[field.title] = value
        return odata

    def save(self, top: Path, force: bool = False):
        text = json.dumps(self.to_dict(), indent=4)
        save_path = self.get_path(top)
        if need_write(save_path, force):
            with open(save_path, write_mode(force)) as f:
                f.write(text)
            logger.action(f"Wrote {self} to {save_path}")
        return save_path

    def __setattr__(self, key: str, value: Any):
        """ Validate the attribute name and value before setting it. """
        # Cache a set of the field keys for fast membership checking.
        if key not in self.get_field_keys_set():
            raise ReportDoesNotHaveFieldError(f"{type(self).__name__}.{key}")
        super().__setattr__(key, value)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.to_dict() == other.to_dict()


class RefReport(Report, RefFileIO, ABC):

    @classmethod
    def get_ident_report_fields(cls):
        return [*super().get_ident_report_fields(), RefF]


class RegReport(RefReport, RegFileIO, ABC):

    @classmethod
    def get_ident_report_fields(cls):
        return [*super().get_ident_report_fields(), RegF]


class BatchedReport(Report, ABC):
    """ Report with a number of data batches (one file per batch). """

    @classmethod
    def get_checksum_report_fields(cls):
        return [NumBatchesF,
                ChecksumsF,
                *super().get_checksum_report_fields()]

    @classmethod
    @abstractmethod
    def _get_batch_types(cls) -> list[type[ReadBatchIO]]:
        """ Type(s) of batch(es) for the report. """

    @classmethod
    @cache
    def get_batch_types(cls) -> dict[str, type[ReadBatchIO]]:
        """ Type(s) of batch(es) for the report, keyed by name. """
        return {batch_type.btype(): batch_type
                for batch_type in cls._get_batch_types()}

    @classmethod
    def get_batch_type(cls, btype: str | None = None) -> type[ReadBatchIO]:
        """ Return a valid type of batch based on its name. """
        if btype is None:
            batch_types = list(cls.get_batch_types().values())
            if (ntypes := len(batch_types)) != 1:
                raise ValueError("btype is optional only if there is exactly "
                                 f"one type of batch, but got {ntypes} types")
            return batch_types[0]
        return cls.get_batch_types()[btype]
