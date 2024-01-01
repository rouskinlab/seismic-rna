
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
