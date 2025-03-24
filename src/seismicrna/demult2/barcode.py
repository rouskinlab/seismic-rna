import ahocorasick
from tqdm import tqdm
from functools import cached_property

from collections import namedtuple, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
import numpy as np

from ..core.seq.xna import DNA
from ..core.seq.refs import RefSeqs

from .neighbor import get_neighbors

from ..core.logs import logger

# Fields in the regions file.
FIELD_REF = "Reference"
FIELD_BARCODE5 = "Barcode5"
FIELD_BARCODE3 = "Barcode3"
FIELD_NAME = "Name"
FIELD_BARCODE = "Barcode"
FIELD_READ_POS = "Read Position"

RegionTuple = namedtuple("RegionTuple", ("pos5", "pos3"))


def get_ref_barcodes(ref_meta_file: Path):
    """
    Parse a file defining each barcode by the name of its reference and
    either its 5' and 3' coordinates or its sequence and name.
    Return one map from each reference to a 5'/3' coordinate
    pair and read position and another from each name to sequence
    and read position.

    Parameters
    ----------
    ref_meta_file: Path
        CSV file of a table that defines the basrcodes. The table must
        have columns labeled "Reference", "Barcode5", "Barcode3",
        "Name", "Barcode", and "Read Position". Others are ignored.

    Returns
    -------
    dict[str, tuple[int, int, int]]
        A mapping from ref to coords
    dict[str, tuple[DNA, int]]
        A mapping from name to barcode
    """

    # Initialize dictionaries mapping references to barcodes
    coords: dict[str, tuple[int, int, int]] = dict()
    bcs: dict[str, tuple[DNA, int]] = dict()

    def map_bc(mapping: dict[str, tuple], key: DNA, value: str):
        """ Add one region to the map if not already present. """
        # Check whether the mapping already contains the key.
        try:
            prev = mapping[key]
        except KeyError:
            # The mapping does not already contain the key: add it.
            mapping[key] = value
        else:
            # Check whether the value already mapped by the key matches
            # the value given currently.
            if prev == value:
                # If so, then warn about it.
                logger.warning(
                    f"Key {key} mapped to {repr(value)} multiple times"
                )
            else:
                # If not, then raise an error because it is ambiguous
                # which value to use.
                raise ValueError(
                    f"Key {key} mapped to {repr(prev)} and to {repr(value)}"
                )

    # Read every row of the regions file.
    refs_meta = pd.read_csv(ref_meta_file)
    lines = zip(refs_meta[FIELD_REF],
                refs_meta[FIELD_BARCODE5],
                refs_meta[FIELD_BARCODE3],
                refs_meta[FIELD_NAME],
                refs_meta[FIELD_BARCODE],
                refs_meta[FIELD_READ_POS])
    refs = set()
    for i, (ref, end5, end3, name, bc, read_pos) in enumerate(lines, start=1):
        end5 -= 1
        if not pd.isnull(ref):
            if ref in refs:
                raise ValueError(f"{ref} may only appear once as a reference or name")
        if not pd.isnull(name):
            if name in refs:
                raise ValueError(f"{name} may only appear once as a reference or name")

        try:
            if not pd.isnull(bc) and pd.isnull(read_pos):
                raise ValueError(f"To directly specify barcodes, {FIELD_BARCODE} and {FIELD_READ_POS} columns must be populated")

            ref_opts = [not pd.isnull(ref), not pd.isnull(end5), not pd.isnull(end3)]
            if any(ref_opts) and not all(ref_opts):
                raise ValueError(f"To extract barcodes from references, {FIELD_REF}, {FIELD_BARCODE5} and {FIELD_BARCODE3} columns must be populated")

            if not pd.isnull(ref):
                ref = str(ref)
            elif not pd.isnull(name):
                ref = str(name)
            else:
                raise ValueError(f"Either {FIELD_REF} or {FIELD_NAME} column must be filled")

            if pd.isnull(read_pos):
                read_pos = end5 # TODO: Fix this

            # Check whether coordinates or primers were given.
            has_coords = not (pd.isnull(end5) or pd.isnull(end3))
            has_bc = not (pd.isnull(bc) or pd.isnull(read_pos))
            if has_coords and has_bc:
                raise ValueError(f"Got both coordinates ({end5}, {end3}) "
                                 f"and barcode ({bc}, staring at {read_pos})")
            elif has_coords:
                # Map the reference and coordinates to the region.
                map_bc(coords, ref, (int(end5), int(end3), int(read_pos)))
            elif has_bc:
                # Map the reference and primers to the region.
                map_bc(bcs, ref, (DNA(bc), read_pos))
            else:
                raise ValueError("Got neither coordinates nor barcodes")
        except Exception as error:
            logger.error(error)
    return coords, bcs


def get_coords_by_name(coords: Iterable[tuple[str, int | DNA, int, int | None]]):
    coords_map: dict[str, set[tuple[int | DNA, int, int | None]]] = defaultdict(set)
    for coord in coords:
        if len(coord) == 4:
            ref, end5, end3, start_pos = coord
            coord_val = end5-1, end3
        elif len(coord) == 3:
            ref, bc, start_pos = coord
            coord_val = bc, start_pos
        else:
            raise ValueError(f"Unexpected number of values for coord {coord}")
        if coord_val in coords_map[ref]:
            logger.warning(f"Skipping duplicate coordinates: {coord}")
        else:
            coords_map[ref].add(coord_val)
    return coords_map


def coords_to_seq(seq: DNA, end5: int, end3: int):
    """ Extract the sequence between inclusive, 1-indexed coordinates """
    return seq[end5:end3]

def expand_by_tolerance(input: tuple[int, ...], tolerance: int) -> set[int]:
    """
    Get the full set of values in input within tolerance
    """
    return {inp + delta for inp in input for delta in range(-tolerance, tolerance + 1)}


class RefBarcodes(object):
    """ A collection of barcodes, mapped to references/names. """

    def __init__(self,
                 ref_seqs: Iterable[tuple[str, DNA]], *,
                 refs_meta_file: Path | None = None,
                 coords: Iterable[tuple[str, int, int, int]] = (),
                 bcs: Iterable[tuple[str, DNA, int]] = (),
                 mismatches: int = 0,
                 index_tolerance: int = 0,
                 allow_n: bool = False):
        ref_seqs = RefSeqs(ref_seqs)
        # Group coordinates and primers by reference.
        cli_coords = get_coords_by_name(coords)
        cli_bcs = get_coords_by_name(bcs)
        # Get the names of the regions from the regions file, if any.
        meta_coords = dict()
        meta_bcs = dict()
        if refs_meta_file is not None:
            try:
                meta_coords, meta_bcs = get_ref_barcodes(refs_meta_file)
                # Combines the coordinates from the regs_file and from the
                # coord parameter.
            except Exception as error:
                logger.error(error)
        all_names  = [key for d in (meta_coords, meta_bcs, cli_coords, cli_bcs) for key in d]
        all_unique = len(all_names) == len(set(all_names))
        if not all_unique:
            raise ValueError(f"Duplicate demultiplexing names/references in {all_names}")

        coords = meta_coords | cli_coords

        bcs_from_coords = dict()
        for ref, coord in coords.items():
            end5, end3, read_pos = coord
            bcs_from_coords[ref] = (coords_to_seq(ref_seqs.get(ref), end5, end3), read_pos)

        self.ref_lengths = {ref: len(seq) for ref, seq in ref_seqs.iter()}
        self._bcs = cli_bcs | meta_bcs | bcs_from_coords
        self.mismatches = mismatches
        self.index_tolerance = index_tolerance
        self.allow_n = allow_n
        self.rc = True # TODO Handle optional RC for barcodes

        self.num_refs = len(self.ref_lengths)

        self.automaton_index = 0

        self.automaton
        if self.rc:
            self.rc_automaton

    @property
    def names(self):
        """ Reference names. """
        return list(self._bcs.keys())

    @property
    def uniq_names(self):
        return set(self.names)

    @property
    def as_dict(self):
        """ Get a dict of name:barcode pairs. """
        return self._bcs

    @property
    def barcodes(self):
        """ List all barcodes. """
        return [bc for bc, _ in self._bcs.values()]

    @property
    def read_positions(self):
        """ List all read positions. """
        return [read_pos for _, read_pos in self._bcs.values()]

    @property
    def max_barcode_len(self):
        return max(map(len, self.barcodes))

    @property
    def read_pos_range(self):
        """ List the range in which a barcode can fall in a read. """
        min_pos = max(min(self.read_positions) - self.index_tolerance, 0)
        max_pos = max(self.read_positions) + self.max_barcode_len + self.index_tolerance # TODO: Max barcode len is not with max read_position
        return (min_pos, max_pos)

    @property
    def slice_position(self):
        slice_positions = list()
        read_range_start, read_range_end = self.read_pos_range
        range_size = read_range_end - read_range_start
        for read_position, barcode in zip(self.read_positions, self.barcodes):
            slice_position = (read_position - read_range_start + len(barcode)) - 1
            slice_positions.append((slice_position, slice_position + range_size + 1)) # +1 for the space character inserted between reads.
        return slice_positions

    @property
    def rc_read_positions(self):
        """ List all reverse complement barcode positions. """
        return [self.ref_lengths.get(name, 0)-(read_pos+len(bc)) for name, (bc, read_pos) in self._bcs.items()] #TODO: Fix rc indexing of non-ref barcodes


    @property
    def rc_read_pos_range(self):
        """ List the range in which an rc barcode can fall in a read. """
        min_pos =  max(min(self.rc_read_positions) - self.index_tolerance, 0)
        max_pos = max(self.rc_read_positions) + self.max_barcode_len + self.index_tolerance
        return (min_pos, max_pos)


    @property
    def rc_slice_position(self):
        slice_positions = list()
        read_range_start, read_range_end = self.rc_read_pos_range
        range_size = read_range_end - read_range_start
        for read_position, rc_barcode in zip(self.rc_read_positions, self.rc_barcodes): # TODO: Make more efficient. Also combine with normal barcodes to avoide code duplication.
            slice_position = (read_position - read_range_start + len(rc_barcode)) - 1
            slice_positions.append((slice_position, slice_position + range_size + 1)) # +1 for the space character inserted between reads.
        return slice_positions


    @property
    def count(self):
        """ Total number of barcodes. """
        return len(self._bcs)

    @property
    def rc_barcodes(self):
        """ Reverse complement of the barcodes. """
        return [bc.rc for bc in self.barcodes]

    @property
    def pairs(self):
        return [pair for pair in zip(self.barcodes, self.read_positions)]

    @property
    def by_pos(self):
        by_pos_dict = defaultdict(list)
        for barcode, pos in self.pairs:
            by_pos_dict[pos].append(barcode)
        return by_pos_dict

    @property
    def rc_pairs(self):
        return [pair for pair in zip(self.rc_barcodes, self.rc_read_positions)]

    @property
    def rc_by_pos(self):
        by_pos_dict = defaultdict(list)
        for barcode, pos in self.rc_pairs:
            by_pos_dict[pos].append(barcode)
        return by_pos_dict


    def get_automaton(self, barcodes: list[tuple[str, DNA]], start: int = 0, check: list[ahocorasick.Automaton] = []): # TODO: May be able to remove the check argument depending on edge cases.
        automaton = ahocorasick.Automaton()
        collisions = dict()
        sets = set()
        orig_set = set()
        for barcode_idx, (name, barcode) in tqdm(enumerate(barcodes)):
            orig_set.add(str(barcode))
            if (name, barcode) in sets:
                raise ValueError(f"Already encountered {(name, barcode)}")
            if self.mismatches:
                barcode_set = get_neighbors(barcode, max_mismatches=self.mismatches, allow_n=self.allow_n)
            else:
                barcode_set = set([barcode])
            for barcode in barcode_set:
                if (name, barcode) in sets:
                    raise ValueError(f"Already encountered {(name, barcode)}")
                sets.add((name, barcode))
                barcode_set = set([str(barcode)])
                for idx, barcode in enumerate(barcode_set):
                    for check_automaton in check:
                        if check_automaton.exists(barcode):
                            assert (name, barcode) not in collisions, f"{(name, barcode)} already collided. {collisions[(name, barcode)]} {check_automaton.get(barcode)}"
                            collisions[(name, barcode)] = check_automaton.get(barcode)
                    automaton.add_word(str(barcode), barcode_idx + start)
        if collisions:
            raise ValueError(f"The following barcodes collide with --mismatch-tolerance {self.mismatches}: "
                             f"{[f'{key} and {val}' for key, val in collisions.items()]}")
        automaton.make_automaton()
        self.automaton_index += len(barcodes)
        return automaton

    @cached_property
    def name_map(self):
        if self.rc:
            names = self.names + self.names
            num_refs = self.num_refs * 2
        else:
            names = self.names
            num_refs = self.num_refs
        name_map = {ref_idx: name for ref_idx, name in zip(np.arange(num_refs), names)}
        return name_map

    @cached_property
    def valid_positions(self):
        if self.rc:
            slice_positions = self.slice_position + self.rc_slice_position
            num_refs = self.num_refs * 2
        else:
            slice_positions = self.slice_position
            num_refs = self.num_refs

        if self.index_tolerance > 0:
            valid_positions = {ref_idx: expand_by_tolerance(read_pos, self.index_tolerance) for ref_idx, read_pos in zip(np.arange(num_refs), slice_positions)}
        else:
            valid_positions = {ref_idx: read_pos for ref_idx, read_pos in zip(np.arange(num_refs), slice_positions)}
        return valid_positions


    @cached_property
    def automaton(self):
        return self.get_automaton([(name, barcode) for name, barcode in zip(self.names, self.barcodes)], start=self.automaton_index)

    @cached_property
    def rc_automaton(self):
        return self.get_automaton([(name, rc_barcode) for name, rc_barcode in zip(self.names, self.rc_barcodes)], start=self.automaton_index)

    def __str__(self):
        return f"{type(self).__name__} ({self.count}): {self._bcs}"
