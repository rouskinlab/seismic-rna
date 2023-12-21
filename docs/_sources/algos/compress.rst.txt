
Algorithm for Sequence Compression/Decompression
========================================================================

Background on Sequence Compression/Decompression
------------------------------------------------------------------------

In a nucleic acid sequence, each base is represented by one character in
`ASCII code`_ and hence occupies 1 byte (8 bits).
But because each base has only 4 possibilities, each base requires only
log\ :sub:`2`\ (4) = 2 bits, shown in the following table:

==== =================
Base 2-Bit Binary Code
==== =================
A    ``00``
C    ``01``
G    ``10``
T/U  ``11``
==== =================

For efficient storage of (long) nucleic acid sequences, the following
algorithm compresses 4 bases (= 8 bits/byte ÷ 2 bits/nt) in each byte,
and then restores the original sequence when needed.

Algorithm for Sequence Compression
------------------------------------------------------------------------

Algorithm for Sequence Compression: Procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, attributes of the nucleic acid sequence are recorded:

- rna (``bool``): whether the sequence is RNA
- length (``int``): number of bases in the sequence
- ns (``tuple[int, ...]``): positions of ``N`` bases, if any

The type and the length are scalar values with constant sizes regardless
of the length of the sequence.
The 0-indexed positions of ambiguous bases is an array that at worst (if
all bases are ``N``) scales linearly with the length of the sequence and
at best (if no bases are ``N``) is a constant (small) size.
Because most reference sequences contain no or very few bases that are
``N``, recording the ``N`` positions generally requires little space.

Then, each non-overlapping segment of four bases is encoded as one byte
by concatenating the 2-bit codes (above) of the bases in the segment, in
reverse order.
Ambiguous (``N``) bases are arbitrarily encoded as ``00``.
Because this step requires the length of the sequence to be a multiple
of 4, the sequence is padded on its 3' side with ``A`` (an arbitrary
choice) until its length is a multiple of 4.

Algorithm for Sequence Compression: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given the DNA sequence ``CAGNTTCGAN``, the attributes are extracted:

- rna: ``False``
- length: ``10``
- ns: ``(3, 9)`` (note: these positions are 0-indexed)

The sequence is then padded with ``A`` at the 3' end until its length
becomes a multiple of 4 (in this case, length 12):

``CAGNTTCGANAA``

Then, each 4-base segment is transformed into one byte by encoding each
base and concatenating the codes in order from code 4 to code 1:

====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========
Number Sequence Base 1 Base 2 Base 3 Base 4 Code 4 Code 3 Code 2 Code 1 Byte
====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========
1      CAGN     C      A      G      N      00     10     00     01     00100001
2      TTCG     T      T      C      G      10     01     11     11     10011111
3      ANAA     A      N      A      A      00     00     00     00     00000000
====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========

Thus, the compressed byte string is ``[00100001, 10011111, 00000000]``.
Note that this string is only 3 bytes, compared to 10 for the original.
As a ``bytes`` object, the representation is ``b'!\x9f\x00'``.

Algorithm for Sequence Decompression
------------------------------------------------------------------------

Algorithm for Sequence Decompression: Procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Beginning with a compressed sequence of bytes, each byte is transformed
into four bases by decoding each 2-bit chunk to a base using the above
table, then reversing the order of the four bases.
The code ``11`` is decoded to ``U`` if the `rna` attribute is ``True``
and to ``T`` if ``False``.
The sequence at this point must have a number of bases divisible by 4.
It is cut to the correct number of bases using the `length` attribute.
Finally, every position in the attribute `ns` is masked to ``N``.

Algorithm for Sequence Decompression: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose that a compressed sequence has the following attributes:

- compressed byte string: ``[00100001, 10011111, 00000000]``
- rna: ``False``
- length: ``10``
- ns: ``(3, 9)``

To decompress the sequence, each byte is split into four 2-bit segments,
decoded, reversed, and reassembled:

====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========
Number Byte     Code 1 Code 2 Code 3 Code 4 Base 4 Base 3 Base 2 Base 1 Sequence
====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========
1      00100001 00     10     00     01     C      A      G      A      CAGA
2      10011111 10     01     11     11     T      T      C      G      TTCG
3      00000000 00     00     00     00     A      A      A      A      AAAA
====== ======== ====== ====== ====== ====== ====== ====== ====== ====== ========

The resulting sequence, ``CAGATTCGAAAA``, is trimmed to 10 nt: ``CAGATTCGAA``.
Finally, (0-indexed) positions 3 and 9 are replaced with ``N``: ``CAGNTTCGAN``.

.. _`ASCII code`: https://en.wikipedia.org/wiki/ASCII
