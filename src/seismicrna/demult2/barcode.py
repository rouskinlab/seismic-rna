import ahocorasick

from collections import namedtuple, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..core.seq.xna import DNA
from ..core.seq.refs import RefSeqs

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
            coord_val = end5, end3
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
    """ Extract the sequence between inclusive coordinates """
    return seq[end5-1:end3]


class RefBarcodes(object):
    """ A collection of barcodes, mapped to references/names. """

    def __init__(self,
                 ref_seqs: Iterable[tuple[str, DNA]], *,
                 refs_meta_file: Path | None = None,
                 coords: Iterable[tuple[str, int, int, int]] = (),
                 bcs: Iterable[tuple[str, DNA, int]] = ()):
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

    @property
    def names(self):
        """ Reference names. """
        return set(self._bcs)

    @property
    def as_dict(self):
        """ Get a dict of name:barcode pairs. """
        return {name: bc for name, (bc, _) in self._bcs.items()}

    @property
    def barcodes(self):
        """ List all barcodes. """
        return [bc for bc, _ in self._bcs.values()]

    @property
    def read_positions(self):
        """ List all read positions. """
        return [read_pos for _, read_pos in self._bcs.values()]

    @property
    def rc_read_positions(self):
        """ List all reverse complement barcode positions. """
        return [self.ref_lengths.get(name, 0)+1-(read_pos+len(bc)) for name, (bc, read_pos) in self._bcs.items()] #TODO: Fix rc indexing of non-ref barcodes

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

    @property
    def automaton(self):
        automaton = ahocorasick.Automaton()
        pairs = self.pairs
        pairs.extend(self.rc_pairs)
        as_dict = self.as_dict
        for name, barcode in as_dict.items():
            automaton.add_word(str(barcode), (name, str(barcode)))
            automaton.add_word(str(barcode.rc), (name, str(barcode.rc)))
        automaton.make_automaton()
        return automaton

    def __str__(self):
        return f"{type(self).__name__} ({self.count}): {self._bcs}"
