
Reference Sequence
------------------------------------------------------------------------

The reference sequence is saved as a ``RefseqIO`` object.

Reference sequence: Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``RefseqIO`` object stores a reference sequence as a space-efficient
``CompressedSeq`` object in the private attribute ``_s``.
(The original sequence is available via the cached property ``refseq``.)
A ``CompressedSeq`` object has the following attributes:

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

For details on the algorithm, see :doc:`../../algos/compress`.
