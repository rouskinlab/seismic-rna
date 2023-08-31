from __future__ import annotations
from functools import cache, wraps
from logging import getLogger
from typing import Callable, TextIO

from ..core.xam import SAM_DELIM, SAM_HEADER, FLAG_PAIRED

logger = getLogger(__name__)


def _reset_seek(func: Callable):
    """ Decorator to reset the position in the SAM file after the
    decorated function returns. """

    @wraps(func)
    def wrapper(sam_file: TextIO, *args, **kwargs):
        previous_position = sam_file.tell()
        try:
            return func(sam_file, *args, **kwargs)
        finally:
            sam_file.seek(previous_position)

    return wrapper


def _range_of_records(get_records_func: Callable):
    @wraps(get_records_func)
    def wrapper(sam_file: TextIO, start: int, stop: int):
        logger.debug(
            f"Reading records from {sam_file.name} from {start} to {stop}")
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
    """ Yield the read name and line for every read in the file. """
    while line := sam_file.readline():
        yield line.split(SAM_DELIM, 1)[0], line, ""


@_range_of_records
def _iter_records_paired(sam_file: TextIO):
    prev_line: str = ""
    prev_name: str = ""
    while line := sam_file.readline():
        if prev_line:
            # Read name is the first field of the line.
            name = line.split(SAM_DELIM, 1)[0]
            # The previous read has not yet been yielded.
            if prev_name == name:
                # The current read is the mate of the previous read.
                yield prev_name, prev_line, line
                prev_line = ""
                prev_name = ""
            else:
                # The previous read is paired, but its mate is not in
                # the SAM file. This situation can occur if Bowtie2 is
                # run in mixed alignment mode and two paired mates fail
                # to align as a pair but one mate aligns individually.
                yield prev_name, prev_line, ""
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_line = line
                prev_name = name
        else:
            # Save the current read so that if its mate is the next
            # read, it will be returned as a pair.
            prev_line = line
            prev_name = line.split(SAM_DELIM, 1)[0]
    if prev_line:
        # In case the last read has not yet been yielded, do so.
        yield prev_name, prev_line, ""


@cache
@_reset_seek
def _find_first_record(sam_file: TextIO):
    """ Return the position of the first record in the SAM file. """
    # Seek to the beginning of the file.
    sam_file.seek(0)
    # Get the first position in the file.
    position = sam_file.tell()
    # Read the file until finding a line that is not part of the header.
    while sam_file.readline().startswith(SAM_HEADER):
        # If the line is part of the header, then get the position at
        # the end of the line.
        position = sam_file.tell()
    # The position at which the first record starts is either the first
    # position in the file (if there are no header lines) or the end of
    # the last header line that was read.
    return position


def _seek_first_record(sam_file: TextIO):
    """ Seek to the beginning of the first record in the SAM file. """
    sam_file.seek(_find_first_record(sam_file))
    return sam_file.tell()


@cache
@_reset_seek
def is_paired(sam_file: TextIO):
    """ Return whether the reads in the SAM file are paired-end. """
    _seek_first_record(sam_file)
    first_line = sam_file.readline()
    try:
        # Try to use the flag of the first read to determine whether the
        # SAM file has paired-end or single-end reads.
        flag = first_line.split(SAM_DELIM, 2)[1]
        paired = bool(int(flag) & FLAG_PAIRED)
    except (IndexError, ValueError):
        # If that failed, default to single-end and issue a warning.
        paired = False
        logger.warning(f"Could not determine whether {sam_file.name} has "
                       f"single- or paired-end reads. Most likely, the file "
                       f"contains no reads. It might also be corrupted.")
    else:
        logger.debug(f"SAM file {sam_file.name} has "
                     f"{'paired' if paired else 'single'}-ended reads")
    return paired


def iter_records(sam_file: TextIO, start: int, stop: int):
    """ Return an iterator of records between positions start and stop
    in the SAM file. """
    return (_iter_records_paired(sam_file, start, stop) if is_paired(sam_file)
            else _iter_records_single(sam_file, start, stop))


def read_name(line: str):
    """ Get the name of the read in the current line of a SAM file. """
    return line.split(SAM_DELIM, 1)[0]


def iter_batch_indexes(sam_file: TextIO, records_per_batch: int):
    """ Yield the start and end positions of every batch in the SAM
    file, where each batch should have about `records_per_batch`
    records. Assume that for nearly all records in paired-end SAM
    files, both mates are present. In the extreme case that only one
    mate is present for every paired-end record, there can be up to
    `2 * records_per_batch` records in a batch. """
    paired = is_paired(sam_file)
    if records_per_batch <= 0:
        raise ValueError(f"records_per_batch must be a positive integer, "
                         f"but got {records_per_batch}")
    logger.info(f"Began computing batch indexes for {sam_file} with "
                f"{'paired' if paired else 'single'}-end reads, aiming for "
                f"{records_per_batch} records per batch")
    # Number of lines to skip between batches: the number of records
    # per batch minus one (to account for the one line that is read
    # at the beginning of each batch, which ensures that every batch
    # has at least one line) times the number of mates per record.
    n_skip = (records_per_batch - 1) * (paired + 1)
    # Start at the beginning of the first record.
    position = _seek_first_record(sam_file)
    # Yield batches until the SAM file is exhausted. If there are no
    # records in the file, then this loop will exit immediately.
    batch = 0
    while line := sam_file.readline():
        # The current batch starts at the current position.
        batch_start = position
        # Skip either the prescribed number of lines or to the end
        # of the file, whichever limit is reached first.
        i_skip = 0
        while i_skip < n_skip and (line := sam_file.readline()):
            i_skip += 1
        # Update the position to the end of the current line.
        position = sam_file.tell()
        # One extra step is required for files of paired-end reads.
        if paired and line:
            # Check if the current line is the mate of the next line.
            if read_name(line) == read_name(sam_file.readline()):
                # If the read names match, then the current line is the
                # mate of the next line: they should be part of the same
                # batch. End the batch after the next read, which is now
                # the current position in the file.
                position = sam_file.tell()
            else:
                # Otherwise, they are not mates. End the batch after the
                # current line.
                sam_file.seek(position)
        # Yield the number and the start and end positions of the batch.
        logger.debug(f"Batch {batch}: {batch_start} - {position}")
        yield batch, batch_start, position
        # Increment the batch number.
        batch += 1
