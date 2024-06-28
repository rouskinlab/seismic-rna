from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path
from subprocess import CompletedProcess

from ..core import path
from ..core.extern import (GUNZIP_CMD,
                           WORD_COUNT_CMD,
                           ShellCommand,
                           args_to_cmd,
                           cmds_to_pipe)

logger = getLogger(__name__)

FQ_LINES_PER_READ = 4
PHRED_ENCS = {33, 64}


def fastq_gz(fastq_file: Path):
    """ Return whether a FASTQ file is compressed with gzip. """
    ext = "".join(fastq_file.suffixes)
    if ext not in path.FastqExt.options:
        raise ValueError(f"Invalid FASTQ extension: {ext}")
    return fastq_file.suffix == ".gz"


def get_args_count_fastq_reads(fastq_file: Path):
    """ Count the reads in a FASTQ file. """
    if fastq_gz(fastq_file):
        return cmds_to_pipe([
            args_to_cmd([GUNZIP_CMD, "--stdout", fastq_file]),
            args_to_cmd([WORD_COUNT_CMD, "-l"])
        ])
    return args_to_cmd([WORD_COUNT_CMD, "-l", fastq_file])


def parse_stdout_count_fastq_reads(process: CompletedProcess):
    """ Parse the output of word count to find the number of reads. """
    n_lines = int(process.stdout.strip().split()[0])
    n_reads, n_extra = divmod(n_lines, FQ_LINES_PER_READ)
    if n_extra:
        raise ValueError(f"Got {n_lines} lines, but expected a multiple of "
                         f"{FQ_LINES_PER_READ}")
    return n_reads


def count_fastq_reads(fastq_file: Path):
    """ Count the reads in a FASTQ file. """
    step = ShellCommand("counting reads in",
                        get_args_count_fastq_reads,
                        parse_stdout_count_fastq_reads,
                        opath=False)
    return step(fastq_file)


def format_phred_arg(phred_enc: int):
    """ Format a Phred score argument for Bowtie2. """
    if phred_enc not in PHRED_ENCS:
        logger.warning(f"Expected phred_enc to be one of {PHRED_ENCS}, "
                       f"but got {phred_enc}, which may cause problems")
    return f"--phred{phred_enc}"


class FastqUnit(object):
    """
    Unified interface for the following sets of sequencing reads:

    - One FASTQ file of single-end reads from one sample
    - One FASTQ file of interleaved, paired-end reads from one sample
    - Two FASTQ files of mate 1 and 2 paired-end reads from one sample
    - One FASTQ file of single-end reads originating from one reference
      sequence in one sample
    - One FASTQ file of interleaved, paired-end reads originating from
      one reference sequence in one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads originating
      from one reference sequence in one sample
    """

    MAX_PHRED_ENC = 2 ** 7 - 1

    KEY_SINGLE = "fastqz"
    KEY_INTER = "fastqy"
    KEY_MATED = "fastqx"
    KEY_DSINGLE = "dmfastqz"
    KEY_DINTER = "dmfastqy"
    KEY_DMATED = "dmfastqx"
    KEY_MATE1 = "fastq1"
    KEY_MATE2 = "fastq2"

    BOWTIE2_FLAGS = {KEY_SINGLE: "-U",
                     KEY_INTER: "--interleaved",
                     KEY_MATE1: "-1",
                     KEY_MATE2: "-2"}

    def __init__(self, *,
                 fastqz: Path | None = None,
                 fastqy: Path | None = None,
                 fastq1: Path | None = None,
                 fastq2: Path | None = None,
                 phred_enc: int,
                 one_ref: bool):
        if fastqz:
            if fastqy or fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.paths: dict[str, Path] = {self.KEY_SINGLE: fastqz}
            self.paired = False
            self.interleaved = False
        elif fastqy:
            if fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.paths: dict[str, Path] = {self.KEY_INTER: fastqy}
            self.paired = True
            self.interleaved = True
        elif fastq1:
            if not fastq2:
                raise TypeError("Got fastq1 but not fastq2")
            self.paths: dict[str, Path] = {self.KEY_MATE1: fastq1,
                                           self.KEY_MATE2: fastq2}
            self.paired = True
            self.interleaved = False
        elif fastq2:
            raise TypeError("Got fastq2 but not fastq1")
        if phred_enc < 0 or phred_enc > self.MAX_PHRED_ENC:
            raise ValueError(f"Invalid Phred encoding: {phred_enc}")
        self.phred_enc = phred_enc
        self.one_ref = one_ref
        self.sample, self.ref, self.exts = self.get_sample_ref_exts()
        logger.debug(f"Instantiated a {type(self).__name__} with "
                     + ", ".join(f"{k} = {v}" for k, v in self.paths.items())
                     + f", phred_enc = {phred_enc}, one_ref = {one_ref}")

    @cached_property
    def phred_arg(self):
        return format_phred_arg(self.phred_enc)

    @cached_property
    def kind(self):
        cls = type(self).__name__
        if self.paired:
            if self.interleaved:
                return f"{cls} of paired-end reads interleaved in one file"
            return f"{cls} of paired-end reads in separate files"
        return f"{cls} of single-end reads in one file"

    @cached_property
    def parent(self):
        """ Return the parent directory of the FASTQ file(s). """
        parents = [inp.parent for inp in self.paths.values()]
        if not parents:
            raise TypeError("Not parent directory")
        if any(parent != parents[0] for parent in parents[1:]):
            raise ValueError("More than one parent directory")
        return parents[0]

    @cached_property
    def seg_types(self) -> dict[str, tuple[path.Segment, ...]]:
        if self.one_ref:
            seg_types = {self.KEY_SINGLE: path.DMFASTQ_SEGS,
                         self.KEY_INTER: path.DMFASTQ_SEGS,
                         self.KEY_MATE1: path.DMFASTQ1_SEGS,
                         self.KEY_MATE2: path.DMFASTQ2_SEGS}
        else:
            seg_types = {self.KEY_SINGLE: path.FASTQ_SEGS,
                         self.KEY_INTER: path.FASTQ_SEGS,
                         self.KEY_MATE1: path.FASTQ1_SEGS,
                         self.KEY_MATE2: path.FASTQ2_SEGS}
        return {key: seg_types[key] for key in self.paths}

    @cached_property
    def n_reads(self) -> int:
        """ Number of reads in the FASTQ file(s). """
        n_reads = list({count_fastq_reads(fq) for fq in self.paths.values()})
        if len(n_reads) != 1:
            raise ValueError(
                f"Expected one unique number of reads, but got {len(n_reads)}")
        return n_reads[0]

    def get_sample_ref_exts(self):
        """ Return the sample and reference of the FASTQ file(s). """
        samples: set[str] = set()
        refs: set[str | None] = set()
        exts: dict[str, str] = dict()
        for key, fq in self.paths.items():
            fq_fields = path.parse(fq, *self.seg_types[key])
            samples.add(fq_fields[path.SAMP])
            refs.add(fq_fields.get(path.REF))
            exts[key] = fq_fields[path.EXT]
        if len(samples) > 1:
            raise ValueError(f"Sample names of {self} disagree: "
                             + " ≠ ".join(samples))
        if len(refs) > 1:
            raise ValueError(f"Ref names of {self} disagree: "
                             + " ≠ ".join(map(str, refs)))
        return list(samples)[0], list(refs)[0], exts

    def fields(self, key: str):
        fields = {path.SAMP: self.sample}
        if self.ref is not None:
            fields[path.REF] = self.ref
        fields[path.EXT] = self.exts[key]
        return fields

    @property
    def cutadapt_input_args(self):
        """ Return input file arguments for Cutadapt. """
        return tuple(self.paths.values())

    @property
    def bowtie2_inputs(self):
        """ Return input file arguments for Bowtie2. """
        return tuple(chain(*[(self.BOWTIE2_FLAGS[key], fq)
                             for key, fq in self.paths.items()]))

    def to_new(self, *new_segments: path.Segment, **new_fields):
        """ Return a new FASTQ unit with updated path fields. """
        new_paths = dict()
        for key, self_path in self.paths.items():
            combined_segments = new_segments + self.seg_types[key]
            combined_fields = self.fields(key) | new_fields
            new_paths[key] = path.build(*combined_segments, **combined_fields)
        return self.__class__(**new_paths,
                              phred_enc=self.phred_enc,
                              one_ref=self.one_ref)

    @classmethod
    def _from_files(cls, /, *,
                    phred_enc: int,
                    one_ref: bool,
                    fqs: list[Path],
                    key: str):
        if key not in (cls.KEY_SINGLE, cls.KEY_INTER):
            raise ValueError(f"Invalid key: '{key}'")
        segs = path.DMFASTQ_SEGS if one_ref else path.FASTQ_SEGS
        for fq in path.find_files_chain(fqs, segs):
            try:
                yield cls(phred_enc=phred_enc, one_ref=one_ref, **{key: fq})
            except Exception as error:
                logger.error(f"Failed to load FASTQ file {fq}: {error}")

    @classmethod
    def _from_mates(cls, /, *,
                    phred_enc: int,
                    one_ref: bool,
                    fqs: list[Path]):
        # Determine the key and segments based on whether the FASTQs are
        # demultiplexed
        if one_ref:
            seg1s = path.DMFASTQ1_SEGS
            seg2s = path.DMFASTQ2_SEGS
        else:
            seg1s = path.FASTQ1_SEGS
            seg2s = path.FASTQ2_SEGS
        # List all FASTQ mate 1 and mate 2 files.
        fq1s = list(path.find_files_chain(fqs, seg1s))
        fq2s = list(path.find_files_chain(fqs, seg2s))

        # Determine the sample and/or reference name of each file.
        def by_tag(fqs_: list[Path], segs: list[path.Segment]):
            tags: dict[tuple[str, str | None], Path] = dict()
            for fq in fqs_:
                fields = path.parse(fq, *segs)
                tag_ = fields[path.SAMP], fields.get(path.REF)
                if tag_ in tags:
                    logger.warning(f"Duplicate sample and reference: {tag_}")
                else:
                    tags[tag_] = fq
            return tags

        tag1s = by_tag(fq1s, seg1s)
        tag2s = by_tag(fq2s, seg2s)
        # Check for any mates with only one file.
        set1s, set2s = set(tag1s), set(tag2s)
        if miss1 := set2s - set1s:
            logger.error(f"Missing FASTQ mate 1 files: {miss1}")
        if miss2 := set1s - set2s:
            logger.error(f"Missing FASTQ mate 2 files: {miss2}")
        # Yield a FASTQ unit for each pair of mated files.
        for tag in set1s & set2s:
            fq_args = {cls.KEY_MATE1: tag1s[tag], cls.KEY_MATE2: tag2s[tag]}
            try:
                yield cls(phred_enc=phred_enc, one_ref=one_ref, **fq_args)
            except Exception as error:
                logger.error(f"Failed to load FASTQ pair {fq_args}: {error}")

    @classmethod
    def from_paths(cls, /, *, phred_enc: int, **fastq_args: list[Path]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.

        Parameters
        ----------
        phred_enc: int
            ASCII offset for encoding Phred scores
        fastq_args: list[Path]
            FASTQ files, given as lists of paths:
            - fastqz: FASTQ files of single-end reads
            - fastqy: FASTQ files of interleaved paired-end reads
            - fastqx: mated FASTQ files of paired-end reads
            - dmfastqz: demultiplexed FASTQ files of single-end reads
            - dmfastqy: demultiplexed FASTQ files of interleaved paired-end reads
            - dmfastqx: demultiplexed mated FASTQ files of paired-end reads

        Yield
        -----
        FastqUnit
            FastqUnit representing the FASTQ or pair of FASTQ files.
            The order is determined primarily by the order of keyword
            arguments; within each keyword argument, by the order of
            file or directory paths; and for directories, by the order
            in which `os.path.listdir` returns file paths.
        """
        # List all FASTQ files.
        # single-end
        yield from cls._from_files(phred_enc=phred_enc,
                                   one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_SINGLE, ()),
                                   key=cls.KEY_SINGLE)
        # interleaved paired-end
        yield from cls._from_files(phred_enc=phred_enc,
                                   one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_INTER, ()),
                                   key=cls.KEY_INTER)
        # mated paired-end
        yield from cls._from_mates(phred_enc=phred_enc,
                                   one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_MATED, ()))
        # demultiplexed single-end
        yield from cls._from_files(phred_enc=phred_enc,
                                   one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DSINGLE, ()),
                                   key=cls.KEY_SINGLE)
        # demultiplexed interleaved paired-end
        yield from cls._from_files(phred_enc=phred_enc,
                                   one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DINTER, ()),
                                   key=cls.KEY_INTER)
        # demultiplexed mated paired-end
        yield from cls._from_mates(phred_enc=phred_enc,
                                   one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DMATED, ()))

    def __str__(self):
        return " ".join(
            [self.kind, " and ".join(map(str, self.paths.values()))]
        )

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
