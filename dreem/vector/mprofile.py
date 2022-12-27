from collections import defaultdict
import itertools
import os
import pathlib
import re
from tqdm import tqdm
from datetime import datetime
from hashlib import file_digest
from multiprocessing import Pool
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd

from dreem.util.dflt import NUM_PROCESSES
from dreem.util.fa import FastaParser
from dreem.util.seq import DNA
from dreem.vector.samview import SamViewer
from dreem.vector.vector import SamRecord


DEFAULT_BATCH_SIZE = 100_000_000  # 100 megabytes


def ref_from_bam_file(bam_file: str):
    return pathlib.PosixPath(bam_file).stem.encode()


def sample_from_bam_file(bam_file: str):
    return os.path.basename(os.path.dirname(os.path.abspath(bam_file)))


def samples_from_bam_files(bam_files: List[str]):
    samples = defaultdict(list)
    for bam_file in bam_files:
        samples[sample_from_bam_file(bam_file)].append(bam_file)
    return dict(samples)


class Region(object):
    __slots__ = ["ref_name", "first", "last", "ref_seq"]

    def __init__(self, ref_name: bytes, first: int, last: int, ref_seq: DNA):
        if first <= 0:
            raise ValueError(f"first ({first}) must be >= 1")
        self.first = first
        if last > len(ref_seq):
            raise ValueError(f"last ({last}) must be <= the length of "
                             f"ref_seq ({len(ref_seq)}) {ref_seq}.")
        if -len(ref_seq) <= last < 0:
            # This option allows using non-positive end coordinates to mean
            # distance from the 3' end (similar to Python's indexing), with
            # -1 meaning up to and including the last coordinate, -2 meaning
            # up to an including the second-to-last coordinate, and so on
            # until -len(ref_seq), which means the first coordinate.
            last += len(ref_seq)
        if last < first:
            raise ValueError(f"last ({last}) must be >= first ({first})")
        self.last = last
        self.ref_seq = ref_seq
        self.ref_name = ref_name
    
    @property
    def spanning(self):
        return self.first == 1 and self.last == len(self.ref_seq)
    
    @property
    def region_seq(self):
        return self.ref_seq[self.first - 1: self.last]

    @property
    def length(self):
        return self.last - self.first + 1

    @property
    def positions(self):
        return np.arange(self.first, self.last + 1)

    @property
    def ref_coords(self):
        return self.ref_name, self.first, self.last
    
    @property
    def columns(self):
        return [f"{chr(base)}{pos}" for base, pos
                in zip(self.ref_seq, self.positions)]
    
    def __eq__(self, other: object) -> bool:
        if isinstance(other, Region):
            return self.ref_coords == other.ref_coords
        else:
            return NotImplemented        


class PrimerRegion(Region):
    primer_gap = 0

    def __init__(self, ref_name: str, ref_seq: DNA,
                 first: Optional[int] = None, last: Optional[int] = None,
                 fwd: Optional[DNA] = None, rev: Optional[DNA] = None):
        if first is None:
            first = (1 if fwd is None else self.locate_primer(
                ref_seq, fwd)[1] + self.primer_gap + 1)
        if last is None:
            last = (len(ref_seq) if rev is None else self.locate_primer(
                ref_seq, rev.rc)[0] - self.primer_gap - 1)
        super().__init__(ref_name, first, last, ref_seq)

    @staticmethod
    def locate_primer(target: DNA, subseq: DNA):
        matches = list(re.finditer(subseq, target))
        if not matches:
            raise ValueError(f"Primer '{subseq}' is not in target '{target}'")
        if len(matches) > 1:
            raise ValueError(f"Primer '{subseq}' occurs {len(matches)} times "
                             f"in target '{target}'")
        return matches[0].start() + 1, matches[0].end()


class MutationalProfile(Region):
    __slots__ = ["sample_name"]

    def __init__(self, sample_name: str, ref_name: str,
                 first: int, last: int, ref_seq: DNA):
        super().__init__(ref_name, first, last, ref_seq)
        self.sample_name = sample_name

    @property
    def short_path(self):
        return os.path.join(self.sample_name,
                            self.ref_name.decode(),
                            f"{self.first}-{self.last}")


class VectorIO(MutationalProfile):
    __slots__ = ["out_dir", "num_batches", "num_vectors", "checksums"]

    digest_algo = "md5"

    def __init__(self, out_dir: str, sample_name: str, ref_name: str,
                 first: int, last: int, ref_seq: DNA, num_batches: int = 0,
                 num_vectors: int = 0, checksums: Optional[List[str]] = None):
        super().__init__(sample_name, ref_name, first, last, ref_seq)
        self.out_dir = os.path.abspath(out_dir)
        self.num_batches = num_batches
        self.num_vectors = num_vectors
        self.checksums = list() if checksums is None else checksums

    @property
    def full_path(self):
        return os.path.join(self.out_dir, self.short_path)

    @property
    def report_file(self):
        return f"{self.full_path}_report.txt"

    @property
    def mv_dir(self):
        return self.full_path

    def get_mv_filename(self, batch_num: int):
        return os.path.join(self.mv_dir, f"{batch_num}.orc")
    
    @property
    def batch_nums(self):
        return range(self.num_batches)

    @property
    def mv_files(self):
        return list(map(self.get_mv_filename, self.batch_nums))

    @classmethod
    def digest_file(cls, path):
        with open(path, "rb") as f:
            digest = file_digest(f, cls.digest_algo).hexdigest()
        return digest


class Report(VectorIO):
    __slots__ = ["speed", "duration", "began", "ended"]

    fields = {"Sample Name": str, "Ref Name": bytes, "First": int, "Last": int,
              "Ref Seq": DNA, "Num Batches": int, "Num Vectors": int,
              "Speed": float, "Duration": float, "Began": datetime,
              "Ended": datetime, "Checksums": list}

    units = {"Speed": "vec/s", "Duration": "s",
             "Checksums": VectorIO.digest_algo}

    datetime_fmt = "%H:%M:%S on %Y-%m-%d"

    def __init__(self, out_dir: str, sample_name: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA, num_batches: int,
                 num_vectors: int, began: datetime, ended: datetime,
                 checksums: List[str], speed: float = 0.0,
                 duration: float = 0.0):
        super().__init__(out_dir, sample_name, ref_name, first, last,
                         ref_seq, num_batches, num_vectors, checksums)
        self.began = began
        self.ended = ended
        if duration > 0.0:
            self.duration = duration
        else:
            dt = self.ended - self.began
            self.duration = dt.seconds + dt.microseconds / 1E6
        if speed > 0.0:
            self.speed = speed
        else:
            try:
                self.speed = self.num_vectors / self.duration
            except ZeroDivisionError:
                self.speed = float("nan")

    @classmethod
    def append_unit(cls, label: str):
        return f"{label} ({unit})" if (unit := cls.units.get(label)) else label

    @classmethod
    def truncate_unit(cls, label: str):
        # Note: this class method only works if every label that has a unit
        # (namely "Speed" and "Duration") contains no whitespace.
        return word if cls.units.get(word := label.split(" ")[0]) else label

    @classmethod
    def attr_to_label(cls, attr: str):
        return cls.append_unit(" ".join(map(str.capitalize, attr.split("_"))))

    @classmethod
    def label_to_attr(cls, label: str):
        return cls.truncate_unit(label).replace(" ", "_").lower()

    @classmethod
    def format_val(cls, val) -> str:
        dtype = type(val)
        if dtype is str or dtype is int or dtype is DNA:
            return val
        if dtype is bytes:
            return val.decode()
        if dtype is float:
            return round(val, 2)
        if dtype is datetime:
            return val.strftime(cls.datetime_fmt)
        if dtype is list:
            return ", ".join(val)
        raise ValueError(dtype)

    @classmethod
    def parse_valstr(cls, label: str, valstr: str):
        dtype = cls.fields[cls.truncate_unit(label)]
        if dtype is str or dtype is int or dtype is float:
            return dtype(valstr)
        if dtype is bytes:
            return valstr.encode()
        if dtype is DNA:
            return DNA(valstr.encode())
        if dtype is datetime:
            return datetime.strptime(valstr, cls.datetime_fmt)
        if dtype is list:
            return valstr.split(", ")
        raise ValueError(dtype)

    def save(self):
        width = max(map(len, map(self.append_unit, self.fields.keys())))
        pattern = "{label: <" + str(width) + "}\t{val}\n"
        with open(self.report_file, "w") as f:
            for label in self.fields.keys():
                attr = self.label_to_attr(label)
                val = self.format_val(self.__getattribute__(attr))
                f.write(pattern.format(label=self.append_unit(label), val=val))
    
    '''
    @classmethod
    def load(cls, report_file):
        vals = {"out_dir": cls.proj_dir_from_file(report_file)}
        with open(report_file) as f:
            for line in f:
                label, valstr = map(str.rstrip, line.split("\t"))
                attr = cls.label_to_attr(label)
                vals[attr] = cls.parse_valstr(label, valstr)
        return cls(**vals)
    '''


class VectorWriter(VectorIO):
    """
    Computes mutation vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """
    __slots__ = ["_bam_file", "_parallel_reads", "_seqbytes"]

    def __init__(self, out_dir: str, bam_file: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA, parallel_reads: bool):
        sample = sample_from_bam_file(bam_file)
        super().__init__(out_dir, sample, ref_name, first, last, ref_seq)
        self._bam_file = bam_file
        self._parallel_reads = parallel_reads
        self._seqbytes = bytes(self.region_seq)

    def _comp_vector(self, rec: SamRecord):
        """
        """
        assert rec.ref_name == self.ref_name
        muts = rec.vectorize(self._seqbytes, self.first, self.last)
        assert muts != bytes(len(muts))
        return muts
    
    def _write_vector_batch(self, muts: np.ndarray, batch_num: int):
        df = pd.DataFrame(data=muts, columns=self.columns, copy=False)
        mv_file = self.get_mv_filename(batch_num)
        df.to_orc(mv_file, engine="pyarrow")
        return mv_file

    def _gen_vector_batch(self, sam_viewer: SamViewer, batch_num: int,
                          start: int, stop: int):
        with sam_viewer as sv:
            mut_bytes = b"".join(map(self._comp_vector,
                                     sv.get_records(start, stop)))
        muts = np.frombuffer(mut_bytes, dtype=np.byte)
        n_records, rem = divmod(len(muts), self.length)
        assert rem == 0
        muts.resize((n_records, self.length))
        mv_file = self._write_vector_batch(muts, batch_num)
        checksum = self.digest_file(mv_file)
        return n_records, checksum

    def _gen_vectors(self):
        os.makedirs(self.mv_dir, exist_ok=False)
        with SamViewer(self._bam_file, self.ref_name, self.first, self.last,
                       self.spanning) as sv:
            batch_size = max(1, DEFAULT_BATCH_SIZE // self.length)
            print(self.length, DEFAULT_BATCH_SIZE, batch_size)
            indexes = list(sv.get_batch_indexes(batch_size))
            starts = indexes[:-1]
            stops = indexes[1:]
            self.num_batches = len(starts)
            assert self.num_batches == len(stops)
            print(starts, stops, self.num_batches)
            svs = [SamViewer(sv.working_path, self.ref_name, self.first,
                             self.last, self.spanning, make=False, remove=False)
                   for _ in self.batch_nums]
            args = list(zip(svs, self.batch_nums, starts, stops))
            if self._parallel_reads:
                n_procs = max(1, min(NUM_PROCESSES, self.num_batches))
                with Pool(n_procs, maxtasksperchild=1) as pool:
                    results = pool.starmap(self._gen_vector_batch, args,
                                           chunksize=1)
            else:
                results = list(itertools.starmap(self._gen_vector_batch, args))
            assert len(results) == self.num_batches
            nums_vectors, self.checksums = map(list, zip(*results))
            print(nums_vectors)
            self.num_vectors = sum(nums_vectors)
            print(self.num_vectors)

    def _write_report(self, t_start, t_end):
        Report(self.out_dir, self.sample_name, self.ref_name, self.first,
               self.last, self.ref_seq, self.num_batches, self.num_vectors,
               t_start, t_end, self.checksums).save()
    
    def gen_vectors(self):
        print(f"Mutational Profile {self.short_path}: computing vectors")
        t_start = datetime.now()
        self._gen_vectors()
        t_end = datetime.now()
        print(f"Mutational Profile {self.short_path}: writing report")
        self._write_report(t_start, t_end)
        print(f"Mutational Profile {self.short_path}: finished")


class VectorWriterSpawner(object):
    __slots__ = ["out_dir", "bam_files", "ref_file", "coords", "primers",
                 "fill", "parallel_profiles", "parallel_reads"]

    def __init__(self,
                 out_dir: str,
                 fasta: str,
                 bam_files: List[str],
                 coords: List[Tuple[bytes, int, int]],
                 primers: List[Tuple[bytes, DNA, DNA]],
                 fill: bool,
                 parallel: str):
        self.out_dir = out_dir
        self.bam_files = bam_files
        self.ref_file = fasta
        self.coords = coords
        self.primers = primers
        self.fill = fill
        if parallel == "auto":
            parallel = ("reads" if self.num_samples == self.num_regions == 1
                        else "profiles")
        if parallel == "profiles":
            self.parallel_profiles = True
            self.parallel_reads = False
        elif parallel == "reads":
            self.parallel_profiles = False
            self.parallel_reads = True
        elif parallel == "off":
            self.parallel_profiles = False
            self.parallel_reads = False
        else:
            raise ValueError(f"Invalid value for parallel: '{parallel}'")
    
    @property
    def num_samples(self):
        return len(samples_from_bam_files(self.bam_files))
    
    @property
    def num_regions(self):
        return len(self.coords) + len(self.primers)
    
    @property
    def ref_seqs(self):
        seqs = dict(FastaParser(self.ref_file).parse())
        if not seqs:
            raise ValueError(f"'{self.ref_file}' contained no sequences")
        return seqs
    
    @property
    def regions(self):
        regions: Dict[bytes, List[PrimerRegion]] = defaultdict(list)
        ref_seqs = self.ref_seqs

        def add_region(region: PrimerRegion):
            if any(region == other for other in regions[region.ref_name]):
                raise ValueError(f"Duplicate region: {region.ref_coords}")
            regions[region.ref_name].append(region)
        
        def ref_to_seq(ref_seqs: Dict[bytes, DNA], ref: bytes):
            try:
                return ref_seqs[ref]
            except KeyError:
                raise ValueError(f"No reference named '{ref.decode()}'")

        for ref, first, last in self.coords:
            add_region(PrimerRegion(ref, ref_to_seq(ref_seqs, ref),
                                    first=first, last=last))
        for ref, fwd, rev in self.primers:
            add_region(PrimerRegion(ref, ref_to_seq(ref_seqs, ref),
                                    fwd=fwd, rev=rev))
        if self.fill:
            for ref, seq in ref_seqs.items():
                if ref not in regions:
                    add_region(PrimerRegion(ref, seq))
        return regions

    @property
    def writers(self):
        regions = self.regions
        no_writers = True
        for bam_file in self.bam_files:
            ref_name = ref_from_bam_file(bam_file)
            for region in regions[ref_name]:
                assert region.ref_name == ref_name
                no_writers = False
                yield VectorWriter(self.out_dir, bam_file, ref_name,
                                   region.first, region.last, region.ref_seq,
                                   self.parallel_reads)
        if no_writers:
            raise ValueError("No samples and/or regions were given.")

    def gen_mut_profiles(self, processes: int = 0):
        if self.parallel_profiles:
            writers = list(self.writers)
            with Pool(processes if processes
                      else min(NUM_PROCESSES, len(writers)),
                      maxtasksperchild=1) as pool:
                pool.map(VectorWriter.gen_vectors, writers,
                         chunksize=1)
        else:
            for writer in self.writers:
                writer.gen_vectors()


'''
class VectorReader(VectorIO):
    @property
    def shape(self):
        return self.num_vectors, self.length

    def run_checksums(self):
        if len(mv_files := self.mv_files) != len(self.checksums):
            raise ValueError(f"Got {len(mv_files)} files but "
                             f"{len(self.checksums)} checksums")
        for mv_file, checksum in zip(mv_files, self.checksums):
            digest = self.digest_file(mv_file)
            if digest != checksum:
                raise ValueError(f"Hex digest of {mv_file} ({digest}) "
                                 f"did not match checksum ({checksum})")

    @property
    def vectors(self):
        self.run_checksums()
        # FIXME: I suspect there are more efficient ways to load these files
        mvs = np.full(self.shape, fill_value=BLANK, dtype=bytes)
        row = 0
        for mv_file in self.mv_files:
            data = orc.read_table(mv_file).to_pandas().values
            n_vectors, n_cols = data.shape
            if n_cols != self.length:
                raise ValueError(f"Expected DataFrame with {self.length}"
                                 f" columns but got {n_cols} columns.")
            mvs[row: (row := row + n_vectors)] = data
        if row != self.num_vectors:
            raise ValueError(f"Expected DataFrame with {self.num_vectors}"
                             f" rows but got {row} rows.")
        vectors = pd.DataFrame(data=mvs, columns=self.positions)
        return VectorSet(self.sample_name, self.ref_name, self.first,
                         self.last, self.ref_seq, vectors)

    @classmethod
    def load(cls, report_file):
        rep = Report.load(report_file)
        return cls(rep.out_dir, rep.sample_name, rep.ref_name, rep.first,
                   rep.last, rep.ref_seq, rep.num_batches, rep.num_vectors,
                   rep.checksums)


class VectorSet(MutationalProfile):
    __slots__ = ["_vectors", "_cover_count", "_mm_count", "_del_count"]

    color_dict = {"A": "red", "C": "blue", "G": "orange", "T": "green"}

    def __init__(self, sample_name: str, ref_name: str, first: int, last: int,
                 ref_seq: DNA, vectors: pd.DataFrame):
        super().__init__(sample_name, ref_name, first, last, ref_seq)
        if (vectors.shape[1] != self.length
                or (vectors.columns != self.positions).any()):
            raise ValueError("Columns of vectors do not match positions.")
        self._vectors = vectors
        self._cover_count = None
        self._mm_count = None
        self._del_count = None

    @staticmethod
    def series(f):
        def wrapper(self):
            return pd.Series(f(self), index=self.positions)
        return wrapper

    @property
    def vectors(self):
        return self._vectors

    @property
    def num_vectors(self):
        return self.vectors.shape[0]

    @property
    def mismatch_frac(self):
        return self.mismatch_count / self.coverage_count

    @property
    @series
    def deletion_count(self):
        if self._del_count is None:
            self._del_count = (self.vectors.values == DELET).sum(axis=0)
        return self._del_count

    @property
    def deletion_frac(self):
        return self.deletion_count / self.coverage_count

    @property
    def mutation_count(self):
        return self.mismatch_count + self.deletion_count

    @property
    def mutation_frac(self):
        return self.mutation_count / self.coverage_count

    @property
    def colors(self):
        return [self.color_dict[chr(base)] for base in self.region_seq]

    def plot_muts(self):
        fig, ax = plt.subplots()
        data = self.mutation_frac
        ax.bar(data.index, data, color=self.colors)
        plt.show()
'''
