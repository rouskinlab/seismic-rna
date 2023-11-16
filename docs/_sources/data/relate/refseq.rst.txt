
Reference Sequence
------------------------------------------------------------------------

The reference sequence is saved as a ``RefseqIO`` object.

Reference sequence: Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``RefseqIO`` object stores a reference sequence as a space-efficient
``CompressedSeq`` object in the private attribute ``_s``.
The original sequence is available via the cached property ``refseq``.

In a standard sequence, each base is represented by one character, and
hence occupies 1 byte (8 bits).
But because each base has only 4 possibilities, each base requires only
log\ :sub:`2`\ (4) = 2 bits, shown in the following table:

==== ===========
Base Binary Code
==== ===========
A    ``00``
C    ``01``
G    ``10``
T    ``11``
==== ===========

A ``CompressedSeq`` object uses this code to compress 4 bases
(= 8 bits/byte ÷ 2 bits/nt) into each byte.
Each has the following attributes:

- ``b`` (``bytes``): sequence of bytes, each byte encoding 4 nucleotides
- ``s`` (``int``): length of the sequence (number of nucleotides)
- ``n`` (``tuple[int, ...]``): 0-indexed position of each ``N``
- ``r`` (``bool``): ``True`` if the sequence is RNA, ``False`` if DNA

Reference sequence: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose that a DNA sequence is ``CAGNTTCGAN``.
This sequence would be compressed into ::

    b = b'!\x9f\x00'
    s = 10
    n = (3, 9)
    r = False

For details on the algorithm, see :doc:`../../algos/compseq`.
