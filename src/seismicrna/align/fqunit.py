from functools import cached_property
from itertools import chain
from pathlib import Path
from subprocess import CompletedProcess
from typing import Iterable

from ..core import path
from ..core.extern import (GUNZIP_CMD,
                           WORD_COUNT_CMD,
                           ShellCommand,
                           args_to_cmd,
                           cmds_to_pipe)
from ..core.io import calc_sha512_path
from ..core.logs import logger
from ..core.ngs import DuplicateSampleReferenceError

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


class MissingFastqMate(FileNotFoundError):
    """ Missing a file in a pair of paired-end FASTQ files. """


class MissingFastqMate1(MissingFastqMate):
    """ Missing mate 1 in a pair of paired-end FASTQ files. """


class MissingFastqMate2(MissingFastqMate):
    """ Missing mate 2 in a pair of paired-end FASTQ files. """


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
    def seg_types(self) -> dict[str, tuple[path.PathSegment, ...]]:
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
                f"Expected one unique number of reads, but got {len(n_reads)}"
            )
        return n_reads[0]

    def get_sample_ref_exts(self):
        """ Return the sample and reference of the FASTQ file(s). """
        samples: set[str] = set()
        refs: set[str | None] = set()
        exts: dict[str, str] = dict()
        for key, fq in self.paths.items():
            fq_fields = path.parse(fq, self.seg_types[key])
            samples.add(fq_fields[path.SAMPLE])
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
        fields = {path.SAMPLE: self.sample}
        if self.ref is not None:
            fields[path.REF] = self.ref
        fields[path.EXT] = self.exts[key]
        return fields

    @property
    def bowtie2_inputs(self):
        """ Return input file arguments for Bowtie2. """
        return tuple(chain(*[(self.BOWTIE2_FLAGS[key], fq)
                             for key, fq in self.paths.items()]))

    def to_new(self, *new_segments: path.PathSegment, **new_fields):
        """ Return a new FASTQ unit with updated path fields. """
        new_paths = dict()
        for key, self_path in self.paths.items():
            combined_segments = new_segments + self.seg_types[key]
            combined_fields = self.fields(key) | new_fields
            new_paths[key] = path.build(combined_segments, combined_fields)
        return self.__class__(**new_paths,
                              phred_enc=self.phred_enc,
                              one_ref=self.one_ref)

    @cached_property
    def checksums(self):
        return {name: calc_sha512_path(path) for name, path in self.paths.items()}

    @classmethod
    def _from_files(cls, /, *,
                    phred_enc: int,
                    one_ref: bool,
                    fqs: Iterable[str | Path],
                    key: str):
        if key != cls.KEY_SINGLE and key != cls.KEY_INTER:
            raise ValueError(f"Invalid key: {repr(key)}")
        if one_ref:
            segs = path.DMFASTQ_SEGS
        else:
            segs = path.FASTQ_SEGS
        sample_refs = set()
        for fq in path.find_files_chain(fqs, segs):
            logger.detail(f"Generating {cls} from {fq}")
            try:
                fq_unit = cls(phred_enc=phred_enc, one_ref=one_ref, **{key: fq})
            except Exception as error:
                logger.error(error)
            else:
                sample_ref = fq_unit.sample, fq_unit.ref
                if sample_ref in sample_refs:
                    raise DuplicateSampleReferenceError(sample_ref)
                logger.detail(f"Generated {cls.__name__} for {sample_ref}")
                yield fq_unit

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

        # Determine the sample and reference name of each file.
        def find_sample_ref(fqs_: list[Path], segs: list[path.PathSegment]):
            sample_refs: dict[tuple[str, str | None], Path] = dict()
            for fq in fqs_:
                fields = path.parse(fq, segs)
                sample_ref_ = fields[path.SAMPLE], fields.get(path.REF)
                if sample_ref_ in sample_refs:
                    raise DuplicateSampleReferenceError(sample_ref_)
                sample_refs[sample_ref_] = fq
            return sample_refs

        sample_ref_1s = find_sample_ref(fq1s, seg1s)
        sample_ref_2s = find_sample_ref(fq2s, seg2s)
        # Check for any mates with only one file.
        set1s, set2s = set(sample_ref_1s), set(sample_ref_2s)
        if missing1 := set2s - set1s:
            raise MissingFastqMate1(missing1)
        if missing2 := set1s - set2s:
            raise MissingFastqMate2(missing2)
        # Yield a FASTQ unit for each pair of mated files.
        for sample_ref in set1s & set2s:
            fq_args = {cls.KEY_MATE1: sample_ref_1s[sample_ref],
                       cls.KEY_MATE2: sample_ref_2s[sample_ref]}
            try:
                fq_unit = cls(phred_enc=phred_enc, one_ref=one_ref, **fq_args)
            except Exception as error:
                logger.error(error)
            else:
                logger.detail(f"Generated {cls.__name__} for {sample_ref}")
                yield fq_unit

    @classmethod
    def from_paths(cls, /, *, phred_enc: int, **fastq_args: Iterable[str | Path]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.

        Parameters
        ----------
        phred_enc: int
            ASCII offset for encoding Phred scores
        fastq_args: list[Path]
            FASTQ files, given as iterables of paths:
            - fastqz: single-end
            - fastqy: interleaved paired-end
            - fastqx: mated paired-end
            - dmfastqz: demultiplexed single-end
            - dmfastqy: demultiplexed interleaved paired-end
            - dmfastqx: demultiplexed mated paired-end

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
        logger.routine(f"Began generating {cls.__name__} instances")
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
        logger.routine(f"Ended generating {cls.__name__} instances")

    def __str__(self):
        return " ".join(
            [self.kind, " and ".join(map(str, self.paths.values()))]
        )

    def __repr__(self):
        return str(self)
