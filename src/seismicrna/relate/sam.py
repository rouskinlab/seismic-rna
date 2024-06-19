from functools import cached_property, wraps
from logging import getLogger
from pathlib import Path
from typing import Callable, TextIO

from ..core import path
from ..core.ngs import (SAM_DELIM,
                        count_total_reads,
                        run_view_xam,
                        run_flagstat,
                        xam_paired)

logger = getLogger(__name__)


def read_name(line: str):
    """ Get the name of the read in the current line of a SAM file. """
    return line.split(SAM_DELIM, 1)[0]


def _range_of_records(get_records_func: Callable):
    @wraps(get_records_func)
    def wrapper(sam_file: TextIO, start: int, stop: int):
        logger.debug(
            f"Reading records from {sam_file.name} from {start} to {stop}"
        )
        sam_file.seek(start)
        records = get_records_func(sam_file)
        n_records = 0
        while True:
            if sam_file.tell() < stop:
                n_records += 1
                yield next(records)
            else:
                if sam_file.tell() > stop:
                    raise ValueError(f"Stopped at position {sam_file.tell()} "
                                     f"of {sam_file.name}, but expected {stop}")
                break
        logger.debug(f"Read {n_records} records in {sam_file.name} "
                     f"from {start} to {stop}")

    return wrapper


@_range_of_records
def _iter_records_single(sam_file: TextIO):
    """ Yield the line for every read in the file. """
    while line := sam_file.readline():
        yield line, ""


@_range_of_records
def _iter_records_paired(sam_file: TextIO):
    """ Yield every pair of lines in the file. """
    prev_line: str = ""
    prev_name: str = ""
    while line := sam_file.readline():
        if prev_line:
            # Read name is the first field of the line.
            name = read_name(line)
            # The previous read has not yet been yielded.
            if prev_name == name:
                # The current read is the mate of the previous read.
                yield prev_line, line
                prev_line = ""
                prev_name = ""
            else:
                # The previous read is paired, but its mate is not in
                # the SAM file. This situation can occur if Bowtie2 is
                # run in mixed alignment mode and two paired mates fail
                # to align as a pair but one mate aligns individually.
                yield prev_line, ""
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_line = line
                prev_name = name
        else:
            # Save the current read so that if its mate is the next
            # read, it will be returned as a pair.
            prev_line = line
            prev_name = read_name(line)
    if prev_line:
        # In case the last read has not yet been yielded, do so.
        yield prev_line, ""


class XamViewer(object):

    def __init__(self, xam_input: Path, tmp_dir: Path, batch_size: int):
        self.xam_input = xam_input
        self.tmp_dir = tmp_dir
        self.batch_size = batch_size

    @cached_property
    def _sample_ref(self):
        fields = path.parse(self.xam_input, *path.XAM_SEGS)
        return fields[path.SAMP], fields[path.REF]

    @property
    def sample(self):
        sample, _ = self._sample_ref
        return sample

    @property
    def ref(self):
        _, ref = self._sample_ref
        return ref

    @cached_property
    def flagstats(self) -> dict:
        return run_flagstat(self.xam_input)

    @cached_property
    def paired(self):
        if (paired := xam_paired(self.flagstats)) is None:
            logger.warning(f"Failed to determine whether {self.xam_input} has "
                           "single- or paired-end reads. Most likely, the file "
                           "contains no reads. It could also be corrupted.")
        return bool(paired)

    @cached_property
    def n_reads(self):
        return count_total_reads(self.flagstats)

    @cached_property
    def tmp_sam_path(self):
        """ Get the path to the temporary SAM file. """
        return path.build(*path.XAM_STAGE_SEGS,
                          top=self.tmp_dir,
                          sample=self.sample,
                          cmd=path.CMD_REL_DIR,
                          stage=path.STAGE_REL_SAMS,
                          ref=self.ref,
                          ext=path.SAM_EXT)

    def create_tmp_sam(self):
        """ Create the temporary SAM file. """
        run_view_xam(self.xam_input, self.tmp_sam_path)

    def delete_tmp_sam(self):
        """ Delete the temporary SAM file. """
        self.tmp_sam_path.unlink(missing_ok=True)

    def open_tmp_sam(self):
        """ Open the temporary SAM file as a file object. """
        if not self.tmp_sam_path.is_file():
            # Create the temporary SAM file if it does not yet exist.
            self.create_tmp_sam()
        return open(self.tmp_sam_path)

    def _iter_batch_indexes(self):
        """ Yield the start and end positions of every batch in the SAM
        file, where each batch should have about `records_per_batch`
        records. Assume that for nearly all records in paired-end SAM
        files, both mates are present. In the extreme case that only one
        mate is present for every paired-end record, there can be up to
        `2 * records_per_batch` records in a batch. """
        if self.batch_size <= 0:
            raise ValueError(f"batch_size must be a positive integer, "
                             f"but got {self.batch_size}")
        logger.info(f"Began computing batch indexes for {self}, aiming for "
                    f"{self.batch_size} reads per batch")
        # Number of lines to skip between batches: the number of records
        # per batch minus one (to account for the one line that is read
        # at the beginning of each batch, which ensures that every batch
        # has at least one line) times the number of mates per record.
        n_skip = (self.batch_size - 1) * (self.paired + 1)
        # Count the batches.
        batch = 0
        with self.open_tmp_sam() as sam_file:
            # Current position in the SAM file.
            position = sam_file.tell()
            while line := sam_file.readline():
                # The current batch starts at the current position.
                batch_start = position
                # Skip either the prescribed number of lines or to the
                # end of the file, whichever limit is reached first.
                i_skip = 0
                while i_skip < n_skip and (line := sam_file.readline()):
                    i_skip += 1
                # Update the position to the end of the current line.
                position = sam_file.tell()
                # Files of paired-end reads require an extra step.
                if self.paired and line:
                    # Check if the current and next lines are mates.
                    if read_name(line) == read_name(sam_file.readline()):
                        # If the read names match, then the lines are
                        # mates and should be in the same batch. Advance
                        # the variable position to point to the end of
                        # the next read.
                        position = sam_file.tell()
                    else:
                        # Otherwise, they are not mates. Backtrack to
                        # the end of the current read, to which the
                        # variable position points.
                        sam_file.seek(position)
                # Yield the number and positions of the batch.
                logger.debug(f"Batch {batch}: {batch_start} - {position}")
                yield batch, batch_start, position
                # Increment the batch number.
                batch += 1

    @cached_property
    def indexes(self):
        return {batch: (start, stop)
                for batch, start, stop in self._iter_batch_indexes()}

    def iter_records(self, batch: int):
        """ Iterate through the records of the batch. """
        start, stop = self.indexes[batch]
        with self.open_tmp_sam() as sam_file:
            if self.paired:
                yield from _iter_records_paired(sam_file, start, stop)
            else:
                yield from _iter_records_single(sam_file, start, stop)

    def __str__(self):
        return f"alignment map {self.xam_input}"

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
