from __future__ import annotations
from functools import cache, wraps
from logging import getLogger
from typing import BinaryIO, Callable

from ..core.xam import SAM_DELIMITER, SAM_HEADER, FLAG_PAIRED


logger = getLogger(__name__)


def _reset_seek(func: Callable):
    """ Decorator to reset the position in the SAM file after the
    decorated function returns. """
    @wraps(func)
    def wrapper(sam_file: BinaryIO, *args, **kwargs):
        prev_pos = sam_file.tell()
        try:
            return func(sam_file, *args, **kwargs)
        finally:
            sam_file.seek(prev_pos)
    return wrapper


def _range_of_records(get_records_func: Callable):
    @wraps(get_records_func)
    def wrapper(sam_file: BinaryIO, start: int, stop: int):
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
        logger.debug(f"Read {n_records} records from {sam_file.name} "
                     f"from {start} to {stop}")
    return wrapper


@_range_of_records
def _iter_records_single(sam_file: BinaryIO):
    """ Yield the read name and line for every read in the file. """
    while line := sam_file.readline():
        yield line.split(SAM_DELIMITER, 1)[0], line, b""


@_range_of_records
def _iter_records_paired(sam_file: BinaryIO):
    prev_line: bytes = b""
    prev_name: bytes = b""
    while line := sam_file.readline():
        if prev_line:
            # Read name is the first field of the line.
            name = line.split(SAM_DELIMITER, 1)[0]
            # The previous read has not yet been yielded.
            if prev_name == name:
                # The current read is the mate of the previous read.
                yield prev_name, prev_line, line
                prev_line = b""
                prev_name = b""
            else:
                # The previous read is paired, but its mate is not in
                # the SAM file. This situation can occur if Bowtie2 is
                # run in mixed alignment mode and two paired mates fail
                # to align as a pair but one mate aligns individually.
                yield prev_name, prev_line, b""
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_line = line
                prev_name = name
        else:
            # Save the current read so that if its mate is the next
            # read, it will be returned as a pair.
            prev_line = line
            prev_name = line.split(SAM_DELIMITER, 1)[0]
    if prev_line:
        # In case the last read has not yet been yielded, do so.
        yield prev_name, prev_line, b""


@cache
@_reset_seek
def _find_first_record(sam_file: BinaryIO):
    """ Return the position of the first record in the SAM file. """
    sam_file.seek(0)
    while (line := sam_file.readline()).startswith(SAM_HEADER):
        pass
    return sam_file.tell() - len(line)


@cache
@_reset_seek
def is_paired(sam_file: BinaryIO):
    """ Return whether the reads in the SAM file are paired-end. """
    sam_file.seek(_find_first_record(sam_file))
    first_line = sam_file.readline()
    try:
        flag = first_line.split(SAM_DELIMITER, 2)[1]
        paired = bool(int(flag) & FLAG_PAIRED)
    except (IndexError, ValueError):
        logger.critical(f"Failed to determine whether {sam_file.name} has "
                        f"single- or paired-end reads. Most likely, the file "
                        f"contains no reads. It might also be corrupted.")
        return False
    logger.debug(f"SAM file {sam_file.name} has "
                 f"{'paired' if paired else 'single'}-ended reads")
    return paired


def iter_records(sam_file: BinaryIO, start: int, stop: int):
    """ Return an iterator of records between positions start and stop
    in the SAM file. """
    return (_iter_records_paired(sam_file, start, stop) if is_paired(sam_file)
            else _iter_records_single(sam_file, start, stop))


def iter_batch_indexes(sam_file: BinaryIO, records_per_batch: int):
    """ Yield the start and end positions of every batch in the SAM
    file, where each batch should have about `records_per_batch`
    records. Assume that for nearly all records in paired-end SAM
    files, both mates are present. In the extreme case that only one
    mate is present for every paired-end record, there can be up to
    `2 * records_per_batch` records in a batch. """
    paired = is_paired(sam_file)
    logger.debug(f"Computing batch indexes for {sam_file} with "
                 f"{'paired' if paired else 'single'}-end reads, aiming for "
                 f"{records_per_batch} records per batch")
    if records_per_batch <= 0:
        raise ValueError(f"records_per_batch must be a positive integer, "
                         f"but got {records_per_batch}")
    # Number of lines to skip between batches: the number of records
    # per batch minus one (to account for the one line that is read
    # at the beginning of each batch, which ensures that every batch
    # has at least one line) times the number of mates per record.
    n_skip = (records_per_batch - 1) * (paired + 1)
    # Start at the beginning of the first record.
    sam_file.seek(_find_first_record(sam_file))
    # Yield batches until the SAM file is exhausted. If there are no
    # records in the file, then this loop will exit immediately.
    batch = 0
    while line := sam_file.readline():
        # The start position of the batch is the beginning of the
        # line that was just read. Since the position in the file is
        # now the end of that line, the line's length is subtracted.
        start_of_batch = sam_file.tell() - len(line)
        # Read either the prescribed number of lines or to the end
        # of the file, whichever limit is reached first.
        for _, line in zip(range(n_skip), sam_file, strict=False):
            pass
        if paired:
            # Get the end position of the current line.
            end_curr = sam_file.tell()
            # Try to get the name of the read in the current line.
            try:
                name_curr = line.split(SAM_DELIMITER, 1)[0]
            except IndexError:
                # Getting the read name failed, so it is impossible
                # to tell whether this read is the mate of the next
                # read. End this batch at the current read so that
                # the next iteration can check whether the next read
                # and the read after the next read are mates.
                logger.error(f"Cannot find read name in '{line.decode()}'")
                logger.debug(f"Batch {batch}: {start_of_batch} - {end_curr}")
                yield start_of_batch, end_curr
                continue
            # Try to read the next line.
            try:
                line = next(sam_file)
            except StopIteration:
                # The end of the file has been reached, so the end
                # of the current line is the end of the batch.
                logger.debug(f"Batch {batch}: {start_of_batch} - {end_curr}")
                yield start_of_batch, end_curr
                continue
            # Get the end position of the next line.
            end_next = sam_file.tell()
            # Try to get the name of the read in the next line.
            try:
                name_next = line.split(SAM_DELIMITER, 1)[0]
            except IndexError:
                # Getting the read name failed, so it is impossible
                # to tell whether this next read is the mate of the
                # current read. End this batch at the next read so
                # that this batch will have an even number of lines,
                # which, when in doubt, is best for paired reads.
                logger.error(f"Cannot find read name in '{line.decode()}'")
                logger.debug(f"Batch {batch}: {start_of_batch} - {end_next}")
                yield start_of_batch, end_next
                continue
            # Determine the end of the batch based on the names of
            # the current and the next reads.
            if name_curr == name_next:
                # If the read names match, then they are mates and
                # should be part of the same batch. So end the batch
                # at the end of the next read.
                logger.debug(f"Batch {batch}: {start_of_batch} - {end_next}")
                yield start_of_batch, end_next
            else:
                # Otherwise, the reads are not mates. If the batch
                # were to end at the end of the next read, then it
                # would be impossible for the next read to be in the
                # same batch as the read after the next read, even
                # if those two reads were mates. To allow those two
                # reads to be compared and potentially be placed in
                # the same batch, the current batch must end after
                # the current read.
                logger.debug(f"Batch {batch}: {start_of_batch} - {end_curr}")
                yield start_of_batch, end_curr
                # Because each iteration of the loop must start at
                # the position at which the last iteration ended,
                # move back to that position in the file.
                sam_file.seek(end_curr)
        else:
            # If the file does not contain paired-end reads, then
            # the end of the batch is the end of the current read.
            end_curr = sam_file.tell()
            logger.debug(f"Batch {batch}: {start_of_batch} - {end_curr}")
            yield start_of_batch, end_curr
        # Increment the batch number.
        batch += 1
