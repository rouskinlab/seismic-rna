from functools import cached_property
from pathlib import Path

from ..core import path
from ..core.extern import cmds_to_pipe
from ..core.logs import logger
from ..core.ngs import (SAM_DELIM,
                        FLAG_PAIRED,
                        FLAG_PROPER,
                        FLAG_UNMAP,
                        FLAG_SECONDARY,
                        FLAG_QCFAIL,
                        FLAG_DUPLICATE,
                        FLAG_SUPPLEMENTARY,
                        ShellCommand,
                        collate_xam_cmd,
                        count_total_reads,
                        run_flagstat,
                        view_xam_cmd,
                        xam_paired)


def get_line_attrs(line: str) -> tuple[str, bool, bool]:
    """ Read attributes from a line in a SAM file. """
    name, flag_str, _ = line.split(SAM_DELIM, 2)
    flag = int(flag_str)
    paired = bool(flag & FLAG_PAIRED)
    proper = bool(flag & FLAG_PROPER)
    return name, paired, proper


def _iter_batch_indexes(sam_file: Path, batch_size: int, paired: bool):
    logger.routine(f"Began calculating batch indexes in {sam_file}")
    if batch_size <= 0:
        raise ValueError(f"batch_size must be a positive integer, "
                         f"but got {batch_size}")
    # Number of lines to skip between batches: the number of records
    # per batch minus one (to account for the one line that is read
    # at the beginning of each batch, which ensures that every batch
    # has at least one line) times the number of mates per record.
    n_skip = (batch_size - 1) * (paired + 1)
    batch = 0
    with open(sam_file) as f:
        # Current position in the SAM file.
        position = f.tell()
        while line := f.readline():
            # The current batch starts at the current position.
            batch_start = position
            # Skip either the prescribed number of lines or to the
            # end of the file, whichever limit is reached first.
            i_skip = 0
            while i_skip < n_skip and (line := f.readline()):
                i_skip += 1
            # Update the position to the end of the current line.
            position = f.tell()
            # Files of paired-end reads require an extra step.
            if paired and line:
                # Check if the current and next lines are mates.
                read_name, read_paired, read_proper = get_line_attrs(line)
                if line := f.readline():
                    # The next line exists.
                    next_name, next_paired, next_proper = get_line_attrs(line)
                    if not (read_paired and next_paired):
                        raise ValueError(f"{sam_file} does not have only "
                                         "paired-end reads")
                    names_match = read_name == next_name
                    if names_match and read_proper != next_proper:
                        raise ValueError(
                            f"Read {repr(read_name)} has only one properly "
                            f"paired mate, which indicates a bug"
                        )
                    if names_match and read_proper:
                        # If the read names match, then the lines are
                        # mates and should be in the same batch. Advance
                        # the variable position to point to the end of
                        # the next read.
                        position = f.tell()
                    else:
                        # Otherwise, they are not mates. Backtrack to
                        # the end of the current read, to which the
                        # variable position points.
                        f.seek(position)
            # Yield the number and positions of the batch.
            logger.detail(
                f"Batch {batch} of {sam_file}: {batch_start} - {position}"
            )
            yield batch, batch_start, position
            # Increment the batch number.
            batch += 1
    logger.routine(f"Ended calculating batch indexes in {sam_file}")


def _iter_records_single(sam_file: Path, start: int, stop: int):
    """ Yield every single-end read in the file. """
    logger.routine(f"Began iterating through single-end records in {sam_file}")
    with open(sam_file) as f:
        f.seek(start)
        while f.tell() < stop:
            line = f.readline()
            name, paired, _ = get_line_attrs(line)
            if paired:
                raise ValueError(
                    f"Read {repr(name)} in {sam_file} is not single-end"
                )
            yield name, line, ""
    logger.routine(f"Ended iterating through single-end records in {sam_file}")


def _iter_records_paired(sam_file: Path, start: int, stop: int):
    """ Yield every paired-end read in the file. """
    logger.routine(f"Began iterating through paired-end records in {sam_file}")
    with open(sam_file) as f:
        f.seek(start)
        prev_line = ""
        prev_name = ""
        prev_proper = False
        while f.tell() < stop:
            line = f.readline()
            name, paired, proper = get_line_attrs(line)
            if not paired:
                raise ValueError(f"Read {repr(name)} in {sam_file} "
                                 "is not paired-end")
            if prev_line:
                # The previous read has not yet been yielded.
                if prev_name == name:
                    # This read and the previous read are mates.
                    if proper != prev_proper:
                        raise ValueError(f"Read {repr(name)} in {sam_file} "
                                         "has only one properly paired mate, "
                                         "which indicates a bug")
                    if proper:
                        # The current read is properly paired with its
                        # mate: yield them together.
                        yield name, prev_line, line
                    else:
                        # The current read is not properly paired with
                        # its mate: yield them separately.
                        yield name, prev_line, prev_line
                        yield name, line, line
                    # Reset the attributes of the previous read.
                    prev_line = ""
                    prev_name = ""
                    prev_proper = False
                else:
                    # The previous read is paired, but its mate is not
                    # in the SAM file.
                    if prev_proper:
                        raise ValueError(f"Read {repr(prev_name)} "
                                         f"in {sam_file} is properly paired "
                                         "but has no mate, "
                                         "which indicates a bug")
                    yield prev_name, prev_line, prev_line
                    # Save the current read so that if its mate is the
                    # next read, it can be returned as a pair.
                    prev_line = line
                    prev_name = name
                    prev_proper = proper
            else:
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_line = line
                prev_name = name
                prev_proper = proper
        if prev_line:
            # In case the last read has not yet been yielded, do so.
            if prev_proper:
                raise ValueError(f"Read {repr(prev_name)} in {sam_file} is "
                                 "properly paired but has no mate, which "
                                 "indicates a bug")
            yield prev_name, prev_line, prev_line
    logger.routine(f"Ended iterating through paired-end records in {sam_file}")


def tmp_xam_cmd(xam_in: Path, xam_out: Path, paired: bool, n_procs: int = 1):
    """ Collate and create a temporary XAM file. """
    flags_req = 0
    flags_exc = (FLAG_UNMAP
                 | FLAG_SECONDARY
                 | FLAG_QCFAIL
                 | FLAG_DUPLICATE
                 | FLAG_SUPPLEMENTARY)
    if paired:
        flags_req |= FLAG_PAIRED
        # Collate the XAM file so paired mates are adjacent.
        collate_step = collate_xam_cmd(xam_in,
                                       None,
                                       tmp_pfx=xam_out.with_suffix(""),
                                       fast=True,
                                       n_procs=max(n_procs - 1, 1))
        # Remove the header.
        view_step = view_xam_cmd(None,
                                 xam_out,
                                 flags_req=flags_req,
                                 flags_exc=flags_exc)
        return cmds_to_pipe([collate_step, view_step])
    # Remove the header.
    flags_exc |= FLAG_PAIRED
    return view_xam_cmd(xam_in,
                        xam_out,
                        flags_req=flags_req,
                        flags_exc=flags_exc)


run_tmp_xam = ShellCommand("collating reads and creating temporary SAM file",
                           tmp_xam_cmd)


class XamViewer(object):

    def __init__(self,
                 xam_input: Path,
                 tmp_dir: Path,
                 branch: str,
                 batch_size: int,
                 n_procs: int = 1):
        self.xam_input = xam_input
        self.tmp_dir = tmp_dir
        self.branch = branch
        self.batch_size = batch_size
        self.n_procs = n_procs

    @cached_property
    def _sample_ref_ancestors_flat(self):
        fields = path.parse(self.xam_input, path.XAM_SEGS)
        return fields[path.SAMPLE], fields[path.REF], fields[path.BRANCHES]

    @property
    def sample(self):
        sample, ref, ancestors_flat = self._sample_ref_ancestors_flat
        return sample

    @property
    def ref(self):
        sample, ref, ancestors_flat = self._sample_ref_ancestors_flat
        return ref

    @cached_property
    def ancestors(self):
        sample, ref, ancestors_flat = self._sample_ref_ancestors_flat
        if len(ancestors_flat) > 1:
            raise ValueError(f"{self} must have ≤ 1 ancestor branches, "
                             f"but got {ancestors_flat}")
        ancestor = ancestors_flat[0] if ancestors_flat else ""
        return {path.ALIGN_STEP: ancestor}

    @cached_property
    def branches(self):
        return path.add_branch(path.RELATE_STEP, self.branch, self.ancestors)

    @cached_property
    def flagstats(self):
        return run_flagstat(self.xam_input, n_procs=self.n_procs)

    @cached_property
    def paired(self):
        """ Whether the reads are paired. """
        return xam_paired(self.flagstats)

    @cached_property
    def n_reads(self):
        """ Total number of reads. """
        return count_total_reads(self.flagstats)

    @cached_property
    def tmp_sam_path(self):
        """ Get the path to the temporary SAM file. """
        return path.build(path.XAM_STAGE_SEGS,
                          {path.TOP: self.tmp_dir,
                           path.SAMPLE: self.sample,
                           path.STEP: path.RELATE_STEP,
                           path.BRANCHES: self.branches,
                           path.STAGE: path.STAGE_REL_SAMS,
                           path.REF: self.ref,
                           path.EXT: path.SAM_EXT})

    def create_tmp_sam(self):
        """ Create the temporary SAM file. """
        if not self.tmp_sam_path.is_file():
            run_tmp_xam(self.xam_input,
                        self.tmp_sam_path,
                        paired=self.paired,
                        n_procs=self.n_procs)

    def delete_tmp_sam(self):
        """ Delete the temporary SAM file. """
        self.tmp_sam_path.unlink(missing_ok=True)

    def _iter_batch_indexes(self):
        """ Yield the start and end positions of every batch in the SAM
        file, where each batch should have about `records_per_batch`
        records. Assume that for nearly all records in paired-end SAM
        files, both mates are present. In the extreme case that only one
        mate is present for every paired-end record, there can be up to
        `2 * records_per_batch` records in a batch. """
        logger.routine(f"Began computing batch indexes for {self}")
        self.create_tmp_sam()
        yield from _iter_batch_indexes(self.tmp_sam_path,
                                       self.batch_size,
                                       self.paired)
        logger.routine(f"Ended computing batch indexes for {self}")

    @cached_property
    def indexes(self):
        return {batch: (start, stop)
                for batch, start, stop in self._iter_batch_indexes()}

    def iter_records(self, batch: int):
        """ Iterate through the records of the batch. """
        logger.routine(f"Began iterating records for {self} batch {batch}")
        start, stop = self.indexes[batch]
        if self.paired:
            yield from _iter_records_paired(self.tmp_sam_path, start, stop)
        else:
            yield from _iter_records_single(self.tmp_sam_path, start, stop)
        logger.routine(f"Ended iterating records for {self} batch {batch}")

    def __str__(self):
        return f"Alignment map {self.xam_input}"
