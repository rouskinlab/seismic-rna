
from typing import Iterable

from .xna import CompressedSeq, XNA


class RefSeqs(object):
    """ Store reference sequences. """

    def __init__(self, seqs: Iterable[tuple[str, XNA]] = ()):
        self._data: dict[str, CompressedSeq] = dict()
        for name, seq in seqs:
            self.add(name, seq)

    def add(self, name: str, seq: XNA):
        """ Add a sequence to the collection via its name. """
        compressed = seq.compress()
        try:
            # Check whether this name was already used for a sequence.
            other = self._data[name]
        except KeyError:
            # If not, then compress and store the sequence.
            self._data[name] = compressed
        else:
            # If so, then confirm all sequences with this name match.
            if compressed != other:
                raise ValueError(f"Got multiple sequences for {repr(name)}: "
                                 f"{seq} ≠ {other.decompress()}")

    def get(self, name: str):
        """ Get a sequence from the collection via its name. """
        return self._data[name].decompress()

    def iter(self):
        """ Yield every sequence and its name. """
        for name, seq in self._data.items():
            yield name, seq.decompress()

    def __iter__(self):
        return self.iter()

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
